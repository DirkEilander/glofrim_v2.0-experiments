# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 15:26:19 2018

@authors: Willem van Verseveld (verseve) & Arjen Haag (haag)
"""

# required modules
import rasterio
from rasterio.mask import mask
from shapely.geometry import box
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
import geopandas as gpd
import os
import shutil
import sys
import configparser as cp
import numpy as np
import xarray as xr
import gdal
from scipy import ndimage as nd
import pandas as pd
import glob
import json
from scipy.optimize import curve_fit

# waterbodies
import pcraster as pcr
from functools import partial
from shapely.ops import transform
import pyproj
import waterbodies as setup_waterbody_maps
import wflow_lake_intbl as setup_lake_intbl
import wflow_reservoir_intbl as setup_reservoir_intbl

# riverwidths
import derive_river_widths as setup_river_widths

# AnnualDischarge parameter
import catchment_FLO1K as flo1k

# upscaled slope
from upscaled_slope import get_slope

# logging
import setup_logging
from datetime import datetime

# modelbuilder
import geojson
import subprocess

# MERIT / pyflwdir
import pyflwdir

from merit.merit_model_data import get_merit_basin_bbox, upscale_merit_basin, network_merit_basin, resample_merit_basin
from merit.wflow_topomaps import wflow_topomaps


def fill(data, invalid=None):
    """
    Replace the value of invalid 'data' cells (indicated by 'invalid') by the value of the nearest valid data cell.
    
    Input:
        data:    numpy array of any dimension
        invalid: binary array of same shape as 'data'. True cells set where data value should be replaced.
                 [if None (default), use: invalid = np.isnan(data)]
    
    Output: filled array
    """
    if invalid is None: invalid = np.isnan(data)
    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]


def transform_coordinates(coords):
    """ 
    Transform coordinates from geodetic to cartesian.
    
    Input:
        coords: a set of lan/lon coordinates (e.g. a tuple or an array of tuples)
    
    Output: same set of coords, converted to cartesian coordinates
    """
    # WGS 84 reference coordinate system parameters
    A = 6378.137 # major axis [km]   
    E2 = 6.69437999014e-3 # eccentricity squared    
    
    coords = np.asarray(coords).astype(np.float)

    # is coords a tuple? Convert it to an one-element array of tuples
    if coords.ndim == 1:
        coords = np.array([coords])
    
    # convert to radiants
    lat_rad = np.radians(coords[:,0])
    lon_rad = np.radians(coords[:,1]) 
    
    # convert to cartesian coordinates
    r_n = A / (np.sqrt(1 - E2 * (np.sin(lat_rad) ** 2)))
    x = r_n * np.cos(lat_rad) * np.cos(lon_rad)
    y = r_n * np.cos(lat_rad) * np.sin(lon_rad)
    z = r_n * (1 - E2) * np.sin(lat_rad)
    
    return np.column_stack((x, y, z))


def getFeatures(gdf):
    """Function to parse features from GeoDataFrame to rasterio"""
    return [json.loads(gdf.to_json())['features'][0]['geometry']]


def reproject_raster(src, src_profile, dst_profile, resampling, threads=1):
    """Reprojects a raster/grid using rasterio.warp.reproject"""
    arr = np.empty((dst_profile['height'],dst_profile['width'])).astype(np.float32)
    rasterio.warp.reproject(
        source = src,
        destination = arr,
        src_crs = src_profile['crs'],
        src_nodata = src_profile['nodata'],
        src_transform = src_profile['transform'],
        dst_transform = dst_profile['transform'],
        dst_crs = dst_profile['crs'],
        resampling = resampling,
        num_threads=threads)
    return arr


def rasterio_mask_shapes(clone, dst_res, src_profile, clone_profile):
    bnds = clone.bounds
    bnd0 = bnds[0]-2*dst_res
    bnd1 = bnds[1]-2*dst_res
    bnd2 = bnds[2]+2*dst_res
    bnd3 = bnds[3]+2*dst_res
    if src_profile['crs'] == clone_profile['crs']:
        bbox = box(bnd0, bnd1, bnd2, bnd3)
        geo = gpd.GeoDataFrame({'geometry': bbox}, index=[0])
    else:
        d = 10
        x = np.concatenate((np.repeat(bnd0, d), np.linspace(bnd0, bnd2,d), np.repeat(bnd2, d), np.linspace(bnd2, bnd0,d)))
        y = np.concatenate((np.linspace(bnd1, bnd3, d), np.repeat(bnd3, d), np.linspace(bnd3, bnd1, d), np.repeat(bnd1, d)))
        bbox = Polygon(zip(x, y))
        geo = gpd.GeoDataFrame({'geometry': bbox}, index=[0])
        geo.crs = clone_profile['crs']
        geo = geo.to_crs(src_profile['crs'])
    coords = getFeatures(geo)
    return coords


def create_M(ksat, theta, zi, method, M_minmax, directory_out, logger):
    """
    Creates M and related map(s).
    
    Parameters:
        ksat : [dict] Ksat maps
        theta : [dict] theta maps
        zi : [array] soilgrid depths 
        method : [int] method to create M, 0 = numpy linalg, 1 = scipy optimize, 2 = both (default)
        M_minmax : [int] value used to constrain M
        directory_out : [string] output path
        logger : [logging object] instance of Python logging module
    
    Output: M and related map(s)
    """
    
    # helper functions
    def func(x, b):
        return np.exp(-b * x)
    
    def contrain_M(M_, popt_0, M_minmax):
        M_[(M_ > 0) & (popt_0 == 0)] = M_minmax
        M_[ M_ > 100000] = M_minmax
        M_[M_ < 0] = M_minmax
        return M_
    
    def do_linalg(z_i_, ks, row, col):
        idx = ~np.isinf(np.log(ks[:,row,col]))
        return np.linalg.lstsq(z_i_[idx,np.newaxis], np.log(ks[idx,row,col]), rcond=None)
    
    def do_curve_fit(z_i_, ks, row, col, p0):
        idx = ~np.isinf(np.log(ks[:,row,col]))
        return curve_fit(func, z_i_[idx], ks[idx,row,col] ,p0=p0)
    
    def write_file(data, path, ext):
        with rasterio.open(path + ext + '.tif', 'w', **dst_profile) as dst:
            dst.write(np.float32(data),1)
    
    def do_translate(path, ext, options='-of PCRaster'):
        gdal.Translate(path + ext + '.map', path + ext + '.tif', options=options)
    
    def write_maps(M_, file_ext, popt_0):
        # write M to file before contrains are applied
        write_file(M_, path_M, '_original' + file_ext)
        do_translate(path_M, '_original' + file_ext)
        # contrain M to specified value(s)
        M_ = contrain_M(M_, popt_0, M_minmax)
        ks_0_new = ks_[0]
        # write M and Ks to files after contrains are applied on M
        write_file(M_, path_M, file_ext)
        write_file(ks_0_new, path_Ks, file_ext)
        do_translate(path_M, file_ext)
        do_translate(path_Ks, file_ext)
    
    # check method
    if method not in [0,1,2]:
        logger.warning('Unrecognized input for method of M parameter calculation! Using default value instead.')
        method = 2
    
    # start M parameter calculation 
    ks_ = []
    for z in zi:
        ks_.append(ksat['KsatVer_' + str(z) + 'cm'].read(1))
    
    ts = theta['thetaS'].read(1)
    tr = theta['thetaR'].read(1)
    
    dst_profile = theta['thetaR'].profile
    
    rows = ks_[0].shape[0]
    cols = ks_[0].shape[1]
    
    ks = (np.array( ks_ )/ks_[0])
    
    z_i_ = zi * 10.0
    
    popt_0 = np.zeros(ks_[0].shape)
    
    path_M = os.path.join(directory_out, 'M')
    path_Ks = os.path.join(directory_out, 'KsatVer')
    
    if (method == 0 or method ==2):
        
        logger.info('fit zi - log(Ksat) with numpy linalg regression (y = b*x) -> M_')
        for row in range(0, rows):
            print(str(round((float(row)/float(rows))*100,2)) + "% completed of curve fit")
            for col in range(0, cols):
                d = do_linalg(z_i_, ks, row, col)
                popt_0[row,col] = d[0]
        M_ = (ts - tr)/(-popt_0)
        write_maps(M_, '_', popt_0)
    
    if (method == 1 or method ==2):
        
        logger.info('fit zi - Ksat with curve_fit (scipy.optimize) -> M')
        for row in range(0, rows):
            print(str(round((float(row)/float(rows))*100,2)) + "% completed of curve fit")
            for col in range(0, cols):
                # try curve fitting with certain p0
                try:
                    popt, pcov = do_curve_fit(z_i_, ks, row, col, (1e-3 ))
                except RuntimeError:
                    # try curve fitting with lower p0
                    try:
                        popt, pcov = do_curve_fit(z_i_, ks, row, col, (1e-4 ))
                    except RuntimeError:
                        # do linalg  regression instead (method 0)
                        popt = np.array(do_linalg(z_i_, ks, row, col))
                        popt[0] = popt[0] * -1.0
                
                popt_0[row,col] = popt[0]
        M_ = (ts - tr)/(popt_0)
        write_maps(M_, '', popt_0)


def make_clone(dst_res, settings, interp_soilthick, M_method, M_minmax, directory_in, directory_out, clonefile, logger):
    """
    Creates maps for wflow model (staticmaps).
    
    Parameters:
        dst_res : [float] resolution of output maps
        settings : [string] path to settings file 
        interp_soilthick : [boolean] control for interpolation/filling of zeros in soil thickness map
        M_method : [int] method to create M, 0 = numpy linalg, 1 = scipy optimize, 2 = both (default)
        M_minmax : [int] value used to constrain M
        directory_in : [string] input path
        directory_out : [string] output path
        clonefile : [string] filename of (PCRaster) clone map file (by default wflow_dem.map)
        logger : [logging object] instance of Python logging module
    
    Output: map files
    """
    
    settings = pd.read_csv(settings)
    
    clone_path = os.path.join(directory_out, clonefile)
    clone = rasterio.open(clone_path)
    
    clone_profile = clone.profile
    
    soilgrids_depths = np.array([0,5,15,30,60,100,200])
    
    for index, row in settings.iterrows():
        
        files = []
        files.extend((glob.glob(os.path.join(directory_in, row.folder_in, row.files))))
        
        for index, filepath in enumerate(files):
            logger.info('Reading ' + filepath)
            src = rasterio.open(filepath)
            src_profile = src.profile
            
            if clone.crs == None:
                logger.warning('*** clone file ' + clone_path + ' without CRS, CRS set to EPSG:4326 ***')
                clone_profile['crs']  = rasterio.crs.CRS.from_epsg(4326)
            
            shapes = rasterio_mask_shapes(clone, dst_res, src_profile, clone_profile)
            
            logger.info('trim file ' + str(filepath))
            out_grid, out_transform = mask(src, shapes, crop=True, all_touched=True)
            
            nx, ny = out_grid[0].shape[1], out_grid[0].shape[0]
            
            grid = out_grid[0] * float(row.mult_factor)
            
            if row.parameter == 'LAI':
                grid[np.where(grid == src_profile['nodata'] * float(row.mult_factor))] = 0.0
            
            if src_profile['nodata'] != None:
                grid[np.where(grid == src_profile['nodata'] * float(row.mult_factor))] = np.nan
                logger.info('fill nodata for parameter ' + row.parameter)
                grid = fill(grid)
            
            scr_file = os.path.basename(files[index])
            dst_tiff_file = scr_file
            
            if scr_file.startswith('KsatVer'):
                dst_map_file = '_'.join(scr_file.split('_')[0:2]) + '.map'
            elif scr_file.startswith('lambda'):
                dst_map_file = ('_'.join(scr_file.split('_')[0:2])).replace('lambda', 'c') + '.map'
            elif scr_file.startswith('LAI'):
                dst_map_file = scr_file.split('_')[0].ljust(8, '0') + '.' +  scr_file.replace('.tif','').split('_')[-1].zfill(3)
            elif scr_file.startswith('GLOBCOVER'):
                dst_map_file = 'wflow_landuse.map'
            elif scr_file.startswith('RootingDepth'):
                dst_map_file = scr_file.replace('tif','') + 'map'
            else:
                dst_map_file = scr_file.split('_')[0] + '.map'
            
            # update the relevant parts of the profiles
            dst_profile = src.meta.copy()
            dst_profile.update({
                'transform': clone_profile['transform'],
                'crs': clone_profile['crs'],
                'dtype' : np.float32,
                'width': clone_profile['width'],
                'height': clone_profile['height']
            })
            
            src_profile.update({
                'transform' : out_transform,
                'width': nx,
                'height': ny
            })
            
            if row.scale_method == 'average':
                resample_method = rasterio.warp.Resampling.average
            if row.scale_method == 'mode':
                resample_method = rasterio.warp.Resampling.mode
            
            if row.conversion == 'log':
                grid = np.log(grid)
            
            logger.info('Resample '+ row.parameter + ' to resolution ' + str(dst_res))
            out = reproject_raster(grid, src_profile, dst_profile, resample_method, threads=4)
            
            if row.conversion == 'log':
                out = np.exp(out)
            
            if (row.parameter == 'soilthickness') and (interp_soilthick):
                logger.info('Interpolating/filling zeros for parameter ' + row.parameter)
                out = fill(out, out==0)
            
            # KsatHorFrac
            if row.parameter == 'KsatHorFrac':
                KsatVer = out_grid[0]
                if src_profile['nodata'] != None:
                    KsatVer[np.where(KsatVer == src_profile['nodata'])] = np.nan
                    KsatVer = fill(KsatVer)
                if row.conversion == 'log':
                    out = out/np.exp(reproject_raster(np.log(KsatVer), src_profile, dst_profile, resample_method, threads=4))
                else:
                    out = out/reproject_raster(KsatVer, src_profile, dst_profile, resample_method, threads=4)
                dst_tiff_file = 'KsatHorFrac.tif'
                dst_map_file = 'KsatHorFrac.map'
            
            if row.parameter == 'lambda':
                logger.info('Convert '+ row.parameter + ' to parameter c')
                out = (3. + (2./out))
                
            path_tif =  os.path.join(directory_out, dst_tiff_file)
            logger.info('write resampled '+ row.parameter + ' to file ' + dst_tiff_file)
            with rasterio.open(path_tif, 'w', **dst_profile) as dst:
                dst.write(out,1)
            
            path_map = os.path.join(directory_out, dst_map_file)
            logger.info('convert ' +  dst_tiff_file  + ' to PCRaster file ' +  dst_map_file)
            gdal.Translate(path_map, path_tif, options = '-of PCRaster')
    
    logger.info('calculating parameter M...')
    files = []
    
    files.extend((glob.glob(os.path.join(directory_out,"KsatVer*.tif"))))
    files.extend((glob.glob(os.path.join(directory_out,"theta*.tif"))))
    
    input_ksat = {}
    input_theta = {}
    for index, filepath in enumerate(files):
        inputfile = os.path.basename(filepath)
        logger.info('read file ' + inputfile )
        if inputfile.startswith('KsatVer'):
            input_ksat['_'.join(inputfile.split('_')[0:2])] = rasterio.open(filepath)
        elif inputfile.startswith('theta'):
            input_theta[inputfile.split('_')[0]] = rasterio.open(filepath)
    
    create_M(input_ksat,input_theta,soilgrids_depths,M_method,M_minmax,directory_out,logger)
    
    del input_ksat, input_theta, clone, src
    
    for f in glob.glob(os.path.join(directory_out, "*.tif")):
        try:
            os.remove(f)
        except:
            logger.error('Could not remove ' + f)
    
    for f in glob.glob(os.path.join(directory_out,"*.aux.xml")):
        try:
            os.remove(f)
        except:
            logger.error('Could not remove ' + f)
    
    c_parameter_files = ['c_5cm.map','c_15cm.map','c_60cm.map','c_200cm.map']
    
    for i,f in enumerate(c_parameter_files):
        try:
            shutil.copy(os.path.join(directory_out,f),os.path.join(directory_out,'c_'+ str(i) + '.map'))
        except:
            logger.error('Could not copy ' + f)
    
    clim_dir = os.path.join(directory_out, "clim")   
    if not os.path.exists(clim_dir):
        os.mkdir(clim_dir)
    
    LAI_files = glob.glob(os.path.join(directory_out, 'LAI*'))     
    for f in LAI_files:
        try:
            shutil.move(f, os.path.join(clim_dir,os.path.basename(f)))
        except:
            logger.error('Could not move ' + f)
    
    logger.info('Creating SoilMinThickness.map by copying and renaming SoilThickness.map')
    try:
        shutil.copy(os.path.join(directory_out,'SoilThickness.map'), os.path.join(directory_out,'SoilMinThickness.map'))
    except:
        logger.error('Could not copy and rename SoilThickness.map')
    
    logger.info('Creating RootingDepth.map by copying and renaming RootingDepth_d75_300x300m.map')
    try:
        shutil.copy(os.path.join(directory_out,'RootingDepth_d75_300x300m.map'), os.path.join(directory_out,'RootingDepth.map'))
    except:
        logger.error('Could not copy and rename RootingDepth.map')


def check_key(dictionary, key, default_value):
    """
    Returns the value assigned to the 'key' in the ini file.
    
    Parameters:
        dictionary : [dict] 
        key : [string|int] 
        default_value : [string] 
    
    Output: Value assigned to the 'key' into the 'dictionary' (or default value if not found)
    """
    if key in dictionary.keys():
        return dictionary[key]
    else:
        return default_value


def replaceLinesIniFile(keys, new, path, logger):
    """
    Replaces a specific line in existing ini file with known contents.
    Uses regular Python text parser instead of ConfigParser to circumvent issues with headers, comments and duplicate sections.
    
    Parameters:
        keys: [list] key(s) / line(s) where to replace value(s)
        new : [list] value(s) to replace
        path : [string] path to ini file
        logger : [logging object] instance of Python logging module
    
    Output: ini file adjusted with certain lines added or removed
    """
    logger.info('Rewriting ini file at ' + path)

    # open .ini file
    with open(path, 'r') as f:
        lines = f.readlines()
    
    # contruct replacement content based on keys
    if 'reservoirs' in keys:
        # text specifically for reservoirs
        txt_for_res = "# Reservoirs\n" \
            "ReserVoirSimpleLocs=staticmaps/wflow_reservoirlocs.map,staticmap,0.0,0\n" \
            "ResTargetFullFrac=intbl/ResTargetFullFrac.tbl,tbl,0.8,0,staticmaps/wflow_reservoirlocs.map\n" \
            "ResTargetMinFrac=intbl/ResTargetMinFrac.tbl,tbl,0.4,0,staticmaps/wflow_reservoirlocs.map\n" \
            "ResMaxVolume=intbl/ResMaxVolume.tbl,tbl,0.0,0,staticmaps/wflow_reservoirlocs.map\n" \
            "ResMaxRelease=intbl/ResMaxRelease.tbl,tbl,1.0,0,staticmaps/wflow_reservoirlocs.map\n" \
            "ResDemand=intbl/ResDemand.tbl,tbl,1.0,0,staticmaps/wflow_reservoirlocs.map\n" \
            "ReservoirSimpleAreas=staticmaps/wflow_reservoirareas.map,staticmap,0.0,0\n" \
            "ResSimpleArea = intbl/ResSimpleArea.tbl,tbl,0,0,staticmaps/wflow_reservoirlocs.map\n"
        # check if present in ini file
        try:
            for txt_to_add_line in txt_for_res.split('\n')[0:-1]:
                index = lines.index(txt_to_add_line+'\n')
            res_in_ini = True
        except:
            res_in_ini = False
        if keys == 'reservoirs_add':
            if not res_in_ini:
                line_below_which_to_add = "LAI=staticmaps/clim/LAI,monthlyclim,1.0,1\n"
                index = lines.index(line_below_which_to_add)
                lines[index] = line_below_which_to_add + '\n' + txt_for_res
            else:
                logger.info('Reservoir parameters already included in model ini file. Skipping rewrite.')
        elif keys == 'reservoirs_rem':
            if res_in_ini:
                for txt_to_rem_line in txt_for_res.split('\n')[0:-1]:
                    index = lines.index(txt_to_rem_line+'\n')
                    lines[index] = ''
                lines[index+1] = '' # to prevent double blank lines between segments
            else:
                logger.info('Reservoir parameters not found in model ini file. Cannot remove.')
    elif 'lakes' in keys:
        # text specifically for lakes
        txt_for_lakes_1 = "# Lakes\n" \
            "LakeLocs=staticmaps/wflow_lakelocs.map,staticmap,0.0,0\n" \
            "LakeAreasMap=staticmaps/wflow_lakeareas.map,staticmap,0.0,0\n" \
            "LinkedLakeLocs=intbl/LinkedLakeLocs.tbl,tbl,0,0,staticmaps/wflow_lakelocs.map\n" \
            "LakeStorFunc=intbl/LakeStorFunc.tbl,tbl,1,0,staticmaps/wflow_lakelocs.map\n" \
            "LakeOutflowFunc=intbl/LakeOutflowFunc.tbl,tbl,3,0,staticmaps/wflow_lakelocs.map\n" \
            "LakeArea=intbl/LakeArea.tbl,tbl,1,0,staticmaps/wflow_lakelocs.map\n" \
            "LakeAvgLevel=intbl/LakeAvgLevel.tbl,tbl,1,0,staticmaps/wflow_lakelocs.map\n" \
            "LakeAvgOut=intbl/LakeAvgOut.tbl,tbl,1,0,staticmaps/wflow_lakelocs.map\n" \
            "LakeThreshold=intbl/LakeThreshold.tbl,tbl,0,0,staticmaps/wflow_lakelocs.map\n" \
            "Lake_b=intbl/Lake_b.tbl,tbl,50,0,staticmaps/wflow_lakelocs.map\n" \
            "Lake_e=intbl/Lake_e.tbl,tbl,2.0,0,staticmaps/wflow_lakelocs.map\n"
        txt_for_lakes_2 = "estimatelakethresh=0\n"
        # check if present in ini file
        try:
            for txt_to_add_line in txt_for_lakes_1.split('\n')[0:-1]:
                index = lines.index(txt_to_add_line+'\n')
            for txt_to_add_line in txt_for_lakes_2.split('\n')[0:-1]:
                index = lines.index(txt_to_add_line+'\n')
            lakes_in_ini = True
        except:
            lakes_in_ini = False
        if keys == 'lakes_add':
            if not lakes_in_ini:
                line_below_which_to_add = "LAI=staticmaps/clim/LAI,monthlyclim,1.0,1\n"
                index = lines.index(line_below_which_to_add)
                lines[index] = line_below_which_to_add + '\n' + txt_for_lakes_1
                line_below_which_to_add_2 = "UStoreLayerThickness = 100,300,800\n"
                index_2 = lines.index(line_below_which_to_add_2)
                lines[index_2] = line_below_which_to_add_2 + txt_for_lakes_2
            else:
                logger.info('Lake parameters already included in model ini file. Skipping rewrite.')
        elif keys == 'lakes_rem':
            if lakes_in_ini:
                for txt_to_rem_line in txt_for_lakes_1.split('\n')[0:-1]:
                    index = lines.index(txt_to_rem_line+'\n')
                    lines[index] = ''
                lines[index+1] = '' # to prevent double blank lines between segments
                for txt_to_rem_line in txt_for_lakes_2.split('\n')[0:-1]:
                    index = lines.index(txt_to_rem_line+'\n')
                    lines[index] = ''
                lines[index+1] = '' # to prevent double blank lines between segments
            else:
                logger.info('Lake parameters not found in model ini file. Cannot remove.')
    else:
        # keys / default lines in .ini file
        key_lines = {
            'AnnualDischarge': 'AnnualDischarge\t= 2290\n',
            'Alpha': 'Alpha\t\t= 120\n'
        }
        # helper function: rewrite key/line
        def rewriteKey(key, value, logger=None):
            logger.info('Rewriting ' + key + ' in ini file...')
            to_replace = key_lines[key].split(' ')[-1]
            try:
                index = lines.index(key_lines[key])
                lines[index] = key_lines[key].replace(to_replace,  str(value) + '\n')
            except ValueError:
                try:
                    index = lines.index(key_lines[key].split(to_replace)[0] + str(value) + '\n')
                    logger.info('Specified replacement did already take place, skipping this rewrite.')
                except ValueError as e:
                    logger.error(str(e))
        # replace specified keys/lines
        for i in range(len(keys)):
            if keys[i] in key_lines.keys():
                rewriteKey(keys[i], new[i], logger)
            else:
                logger.error("Specified key not recognized: '" + key[i] + "'. Cannot replace value for this key in ini file!")
    
    # rewrite .ini file
    with open(path, 'w') as new_f:
        for line in lines:
            new_f.write(line)


def copyFiles(dst, src_files, script_root, src_scripts, logger):
    """
    Copies files from various locations to the specified directory.
    
    Parameters:
        dst : [string] destination path
        src_files : [list] list of source paths (files)
        src_scripts : [list] list of source paths (scripts)
        logger : [logging object] instance of Python logging module
    
    Output: ReadMe.txt file created in specified path
    """
    for src_file in src_files:
        try:
            shutil.copy(src_file, os.path.join(dst, os.path.basename(src_file)))
        except:
            logger.warning('Could not copy ' + src_file)
    for src_script in src_scripts:
        try:
            shutil.copy(os.path.join(script_root, src_script), os.path.join(dst, src_script))
        except:
            logger.warning('Could not copy ' + os.path.join(script_root, src_script))


def createReadMe(path_readme, setup_folder='setup'):
    """
    Create a ReadMe file in each catchment folder, explaining the automated setup.
    
    Parameters:
        path_readme : [string] full path including file extension where to create file
        setup_folder : [string] name of setup folder within the same directory
    
    Output: ReadMe.txt file created in specified path
    """
    lines = [
        "This model was automatically set up using a combination of different Python scripts and data sources. \n\n",
        "Several folders and files have been updated, including 'intbl', 'staticmaps' and the model .ini file. \n",
        "This means that other folders and files, specifically those in the 'data\parameters' folder, might contain outdated information. \n\n",
        "Setup and debug information is stored in the '" + setup_folder + "' folder. This contains, amongst others:\n",
        "- log file, with the following levels:\n",
        "  - INFO:     information deemed relevant during and/or after execution of setup\n"
        "  - WARNING:  potential issue was encountered for which a workaround was in place, should probably be checked\n"
        "  - ERROR:    potential issue was encountered for which no workaround was in place, should definitely be checked\n"
        "  - CRITICAL: could not execute code that is vital for success, no setup/processing could be carried out\n"
        "- ini file\n",
        "- settings file\n",
        "- used scripts\n",
        "- figures related to the set up of lakes and/or reservoir intbl's and the Annual Discharge parameter\n",
    ]
    with open(os.path.join(path_readme, 'ReadMe.txt'), 'w') as read_me:
        for line in lines:
            read_me.write(line)


def changeGeoJSONgeomType(path, file_old='catchments_original.geojson', file_new='catchments_v2.geojson', geomType_old='Polygon', geomType_new='MultiLineString', remove_old=True, logger=None):
    """
    Changes the geometry type in a (Geo)JSON file.
    
    Parameters:
        path : [string] path to directory where old/temporary file is located
        file_old : [string] filename of old/temporary file
        file_new : [string] filename of new/to-be-created file
        geomType_old : [string] geometry type in old/temporary file
        geomType_new : [string] geometry type for new/to-be-created file
        remove_old : [boolean] delete/remove old/temporary file
        logger : [logging object] instance of Python logging module
    
    Output: ReadMe.txt file created in specified path
    """
    logger.info('Changing geometry type in catchments file from ' + geomType_old + ' to ' + geomType_new)
    with open(os.path.join(path, file_old), 'r') as f:
        data = json.load(f)
    for feature in data['features']:
        if feature['geometry']['type'] == geomType_old:
            feature['geometry']['type'] = geomType_new
        else:
            logger.warning('Feature in GeoJSON does not have geometry type ' + geomType_old + '!')
    with open(os.path.join(path, file_new), 'w') as f:
        json.dump(data, f)#, indent=2)
    if remove_old:
        os.remove(os.path.join(path, file_old))


def runModelbuilder(folder, catchment_id, resolution, modelbuilder_path, catchments_path, rivers_path, logger):
    """
    Runs the modelbuilder to obtain new topographic base data.
    
    Parameters:
        folder : [string] path to current catchment/case
        catchment_id : [string] id/name of current catchment/case
        resolution : [float] resolution for model
        modelbuilder_path : [string] path to modelbuilder scripts
        rivers_path : [string] path to existing rivers geojson (potential modelbuilder input)
        logger : [logging object] instance of Python logging module
    
   Output: various files created with modelbuilder, including logs in newly created 'modelbuilder' folder within catchment/case
    """
    # update modelbuilder settings file
    modelbuilder_settings_path = os.path.join(modelbuilder_path, 'settings.json')
    logger.info('Updating modelbuilder settings file at ' + modelbuilder_settings_path)
    logger.info('Using catchment geometry from ' + catchments_path)
    with open(modelbuilder_settings_path) as f:
        modelbuilder_settings = geojson.load(f)
    with open(catchments_path) as c:
        catchment_json = geojson.load(c)
    geom_type = catchment_json['features'][0]['geometry']['type']
    modelbuilder_settings['features'][0]['properties'] = resolution
    modelbuilder_settings['features'][0]['geometry']['type'] = geom_type
    modelbuilder_settings['features'][0]['geometry']['coordinates'] = catchment_json['features'][0]['geometry']['coordinates']
    with open(modelbuilder_settings_path, 'w') as f:
        geojson.dump(modelbuilder_settings, f)
    logger.info('Modelbuilder settings file updated.')
    # clearing modelbuilder log files
    modelbuilder_log_1      = 'wtools_create_grid.log'
    modelbuilder_log_2      = 'wtools_static_maps.log'
    modelbuilder_log_1_path = os.path.join(modelbuilder_path, modelbuilder_log_1)
    modelbuilder_log_2_path = os.path.join(modelbuilder_path, modelbuilder_log_2)
    logger.info('Clearing modelbuilder log files...')
    if os.path.exists(modelbuilder_log_1_path):
        os.remove(modelbuilder_log_1_path)
    if os.path.exists(modelbuilder_log_2_path):
        os.remove(modelbuilder_log_2_path)
    # run modelbuilder from command line
    logger.info('Running modelbuilder...')
    curr_work_dir = os.getcwd()
    try:
        os.chdir(modelbuilder_path)
        logger.info('Changed working directory to ' + modelbuilder_path)
    except OSError:
        logger.fatal('Could not change working directory to ' + modelbuilder_path)
        sys.exit()
    if rivers_path != None:
        if geom_type != 'Point':
            modelbuilder_command = "python modelbuilder.py --geojson-path settings.json --name " + catchment_id + " --cellsize " + str(resolution) + " --river-path " + rivers_path + " --region-filter region"
        else:
            modelbuilder_command = "python modelbuilder.py --geojson-path settings.json --name " + catchment_id + " --cellsize " + str(resolution) + " --river-path " + rivers_path
    else:
        if geom_type != 'Point':
            modelbuilder_command = "python modelbuilder.py --geojson-path settings.json --name " + catchment_id + " --cellsize " + str(resolution) + " --region-filter region"
        else:
            modelbuilder_command = "python modelbuilder.py --geojson-path settings.json --name " + catchment_id + " --cellsize " + str(resolution)
    logger.info('Running the following in command line: ' + modelbuilder_command)
    try:
        subprocess.run(modelbuilder_command)
    except:
        logger.fatal('Modelbuilder command line operation failed!')
        sys.exit()
    # copy/overwrite all modelbuilder data to current folder
    logger.info('Replacing outdated directories of '+ folder + ' with modelbuilder output...')
    to_be_replaced_dirs = os.listdir(os.path.join(modelbuilder_path, catchment_id))
    for temp_dir in to_be_replaced_dirs:
        if os.path.isdir(os.path.join(modelbuilder_path, catchment_id, temp_dir)):
            if os.path.exists(os.path.join(folder, temp_dir)):
                logger.info("Deleting outdated '" + temp_dir + "' folder")
                shutil.rmtree(os.path.join(folder, temp_dir))
            logger.info("Copying modelbuilder '" + temp_dir + "' folder")
            shutil.copytree(os.path.join(modelbuilder_path, catchment_id, temp_dir), os.path.join(folder, temp_dir))
        else:
            logger.info("Overwriting outdated '" + temp_dir + "' file")
            shutil.copy(os.path.join(modelbuilder_path, catchment_id, temp_dir), os.path.join(folder, temp_dir))
    logger.info('Deleting modelbuilder output...')
    shutil.rmtree(os.path.join(modelbuilder_path, catchment_id))
    # copy modelbuilder files to new folder
    modelbuilder_folder = os.path.join(folder, 'modelbuilder')
    if not os.path.exists(modelbuilder_folder):
        os.makedirs(modelbuilder_folder, exist_ok=True)
    else:
        # clear folder if it already existed
        try:
            for temp_file in os.listdir(modelbuilder_folder):
                os.remove(os.path.join(modelbuilder_folder, temp_file))
        except OSError:
            pass
    logger.info('Copying modelbuilder settings and log files to ' + modelbuilder_folder)
    shutil.copy(modelbuilder_settings_path, os.path.join(modelbuilder_folder, 'settings.json'))
    shutil.copy(modelbuilder_log_1_path, os.path.join(modelbuilder_folder, modelbuilder_log_1))
    shutil.copy(modelbuilder_log_2_path, os.path.join(modelbuilder_folder, modelbuilder_log_2))
    # change working directory back to what is was before modelbuilder run
    try:
        os.chdir(curr_work_dir)
        logger.info('Changed working directory back to ' + curr_work_dir)
    except OSError:
        logger.fatal('Could not change working directory back to ' + curr_work_dir)
        sys.exit()
    logger.info('Modelbuilder process finished. Continuing with setup procedure...')


def getWaterbodiesCatchment(catchment, waterbodies_path):
    """
    Obtain waterbodies (i.e. lakes or reservoirs) within catchment.
    
    Parameters:
        catchment : [geopandas object] catchment bounds
        waterbodies_path : [string] path to dataset for lakes or reservoirs
    
    Output: waterbodies within catchment
    """
    #get catchment boundaries
    c = catchment.geometry.bounds
    #load data
    waterbodies = gpd.read_file(waterbodies_path)
    #intersect catchment bounding box with reservoirs
    s_idx = waterbodies.sindex
    inds = list(s_idx.intersection((c.minx[0],c.miny[0],c.maxx[0],c.maxy[0])))
    if len(inds) > 0:
        #intersect actual catchment bounds with reservoirs
        match_intersect = waterbodies.iloc[inds]
        catch_intersect = gpd.overlay(catchment, match_intersect, how='intersection')
        if catch_intersect.geometry.size == 0:
            catch_intersect = gpd.GeoDataFrame([])
    else:
        catch_intersect = gpd.GeoDataFrame([])
    return catch_intersect


def checkWaterbodies(type, catchment, waterbodies_path, min_area=0, logger=None):
    """
    Checks if waterbodies (i.e. lakes or reservoirs) of sufficient size are present within catchment.
    
    Parameters:
        type: [string] 'lake' or 'reservoir'
        catchment : [geopandas object] catchment bounds
        waterbodies_path : [string] path to dataset for lakes or reservoirs
        min_area : [float) minimum area of waterbodies to include
        logger : [logging object] instance of Python logging module
    
    Output: dictionary with waterbodies within catchment and boolean whether to execute processing
    """
    #initialize output variables
    catch_intersect = None
    do_process = False
    #get catchment boundaries
    c = catchment.geometry.bounds
    #load data
    waterbodies = gpd.read_file(waterbodies_path)
    #intersect catchment bounding box with water bodies
    s_idx = waterbodies.sindex
    inds = list(s_idx.intersection((c.minx[0],c.miny[0],c.maxx[0],c.maxy[0])))
    if len(inds) > 0:
        #intersect actual catchment bounds with water bodies
        match_intersect = waterbodies.iloc[inds]
        catch_intersect = gpd.overlay(catchment, match_intersect, how='intersection')
        if catch_intersect.geometry.size > 0:
            if min_area > 0:
                #filter on area
                catch_intersect = catch_intersect[catch_intersect['Lake_area']>=min_area]
                if catch_intersect.geometry.size > 0:
                    do_process = True
                else:
                    logger.warning('No ' + type + 's of sufficient size found within catchment! Skipping ' + type + ' procedures!')
            else:
                do_process = True
            #update area value (because only parts of water bodies might fall within catchment)
            if do_process:
                proj = partial(pyproj.transform, pyproj.Proj(init='epsg:4326'), pyproj.Proj(init='epsg:3857'))
                def updateArea(polygon, proj=proj):
                    return transform(proj, polygon).area / 1e6 # km2
                catch_intersect['Lake_area'] = catch_intersect.apply(lambda row: updateArea(row['geometry']), axis=1)
        else:
            logger.warning('No ' + type + 's found within catchment! Skipping ' + type + ' procedures!')
    else:
        logger.warning('No ' + type + 's found within catchment! Skipping ' + type + ' procedures!')
    return {'data': catch_intersect, 'do_process':do_process}


def processWaterbodies(type, do_process, waterbodies, out_intbl_dir, out_clone_dir, model_ini_path, res_range_min, res_range_max, res_method, debug_path, logger):
    """
    Processing for waterbodies (i.e. lakes or reservoirs).
    
    Parameters:
        type : [string] 'lake' or 'reservoir'
        do_process : [boolean] whether to process/add for specified type
        out_intbl_dir : [string] path to where newly created intbl files will be placed
        out_clone_dir : [string] path to where newly created map files will be placed
        model_ini_path : [string] path to model ini file
        res_range_min : [list|float] range (min,max) for reservoir minimum fractions
        res_range_max : [list|float] range (min,max) for reservoir full fractions
        res_method : [int] method of intbl calculation for reservoirs (0, 1 or 2)
        debug_path : [string] path to where debug/setup information is stored
        logger : [logging object] instance of Python logging module
    
    Output: new intbl and map files (when do_process=True) as well as an adjusted model ini file
    """
    if do_process:
        logger.info('Setting up ' + type + ' intbls in ' + out_intbl_dir)
        if type == 'reservoir':
            setup_reservoir_intbl.make_tbls(waterbodies, out_intbl_dir, range_min=res_range_min, range_max=res_range_max, method=res_method, debug_path=debug_path, logger=logger)
        elif type == 'lake':
            setup_lake_intbl.make_tbls(waterbodies, out_intbl_dir, debug_path=debug_path, logger=logger)
        logger.info('Finished ' + type + ' intbls setup')
        logger.info('Setting up ' + type + ' maps in ' + out_clone_dir)
        if type == 'reservoir':
            setup_waterbody_maps.make_maps(waterbodies, False, out_clone_dir, "wflow_reservoirareas", "wflow_reservoirlocs", logger=logger)
        elif type == 'lake':
            setup_waterbody_maps.make_maps(waterbodies, True, out_clone_dir, "wflow_lakeareas", "wflow_lakelocs", logger=logger)
        logger.info('Finished ' + type + ' maps setup')
        logger.info('Adding ' + type + ' parameters to model ini file...')
        replaceLinesIniFile(type+'s_add', None, model_ini_path, logger=logger)
    else:
        logger.info('Removing ' + type + ' parameters from model ini file...')
        replaceLinesIniFile(type+'s_rem', None, model_ini_path, logger=logger)


def createMapFile(data, outdir, mapfile, map_meta, logger):
    """
    Helper function to create PCRaster .map file by first creating a temporary GeoTIFF file and then converting that to desired output.
    """
    #create temporary GeoTIFF file
    logger.info('Staring process to write ' + os.path.join(outdir, mapfile + '.map'))
    logger.info('Writing temporary GeoTIFF file...')
    map_meta.update(compress='lzw') # add compression
    map_meta['driver']='GTiff' # make sure driver is set to GeoTIFF
    with rasterio.open(os.path.join(outdir, mapfile + '.tif'), 'w', **map_meta) as out:
        out.write(data, 1)
    out.close()
    #convert temporary GeoTIFF file to PCRaster .map file
    logger.info('Converting temporary GeoTIFF file to PCRaster map file...')
    gdal.Translate(os.path.join(outdir, mapfile + '.map'), os.path.join(outdir, mapfile + '.tif'), options = '-of PCRaster')
    #delete temporary GeoTIFF file (and .map.aux.xml file created during conversion)
    logger.info('Deleting temporary GeoTIFF file and .map.aux.xml file that is created during conversion..')
    if os.path.exists(os.path.join(outdir, mapfile + '.tif')):
        os.remove(os.path.join(outdir, mapfile + '.tif'))
    else:
        logger.warning('Could not find temporary GeoTIFF file! Skipping removal.')
    if os.path.exists(os.path.join(outdir, mapfile + '.map.aux.xml')):
        os.remove(os.path.join(outdir, mapfile + '.map.aux.xml'))
    else:
        logger.warning('Could not find .map.aux.xml file created during conversion! Skipping removal.')


def updateRiverWidths(path, reservoirs, lakes, flag=-2, logger=None):
    """
    Removes large river widths from river width map at locations of lakes and reservoirs.
    
    Parameters:
        path : [string] path to current catchment/case maps folder (i.e. directory of river width map)
        reservoirs : [GeoDataFrame] reservoirs within catchment
        lakse : [GeoDataFrame] lakes within catchment
        flag : [int] value to assign to waterbody pixels in riverwidth map
        logger : [logging object] instance of Python logging module
    
    Output: updated wflow_riverwidth map.
    """
    logger.info('Updating wflow_riverwidth.map by replacing lake/reservoir pixels with a value of ' + str(flag))
    #read in current riverwidths map and obtain map metadata
    riverwidth_map = rasterio.open(os.path.join(path, 'wflow_riverwidth.map'), dtype=np.uint)
    riverwidths    = riverwidth_map.read(1)
    map_meta = riverwidth_map.meta.copy()
    if map_meta['crs'] == None:
        map_meta['crs'] = rasterio.crs.CRS.from_epsg(4326)
    if not lakes.empty:
        #rasterize lakes
        lake_ids     = lakes.Hylak_id
        out_arr      = pcr.pcr2numpy(pcr.readmap(os.path.join(path, 'wflow_subcatch.map')), map_meta['nodata'])
        out_arr      = (out_arr/out_arr)-1 #make sure default array contains zeros only
        lake_shapes  = ((geom,value) for geom, value in zip(lakes.geometry, lake_ids))
        lakes_raster = rasterio.features.rasterize(shapes=lake_shapes, fill=0, out=out_arr, transform=map_meta['transform'], all_touched=True)
        #flag riverwidths at lakes
        riverwidths = np.where(lakes_raster > 0, flag, riverwidths)
    if not reservoirs.empty:
        #rasterize reservoirs
        reservoir_ids     = reservoirs.ID
        out_arr           = pcr.pcr2numpy(pcr.readmap(os.path.join(path, 'wflow_subcatch.map')), map_meta['nodata'])
        out_arr           = (out_arr/out_arr)-1 #make sure default array contains zeros only
        reservoir_shapes  = ((geom,value) for geom, value in zip(reservoirs.geometry, reservoir_ids))
        reservoirs_raster = rasterio.features.rasterize(shapes=reservoir_shapes, fill=0, out=out_arr, transform=map_meta['transform'], all_touched=True)
        #flag riverwidths at reservoirs
        riverwidths = np.where(reservoirs_raster> 0 , flag, riverwidths)
    #copy original map
    logger.info('Copying old map to wflow_riverwidth_original.map')
    shutil.copy(os.path.join(path, 'wflow_riverwidth.map'), os.path.join(path, 'wflow_riverwidth_original.map'))
    #save new map
    createMapFile(riverwidths, path, 'wflow_riverwidth_new', map_meta, logger)
    #overwrite original map with updated values
    logger.info('Saving new map as wflow_riverwidth.map')
    riverwidth_map.close()
    os.remove(os.path.join(path, 'wflow_riverwidth.map'))
    shutil.copy(os.path.join(path, 'wflow_riverwidth_new.map'), os.path.join(path, 'wflow_riverwidth.map'))
    os.remove(os.path.join(path, 'wflow_riverwidth_new.map'))


def checkStaticmaps(staticmaps_check, out_clone_dir, logger):
    """
    Checks if all staticmaps specified in setup ini file were created.
    
    Parameters:
        staticmaps_check : [list|string] files that should be present in staticmaps folder
        out_clone_dir : [string] path to staticmaps folder
        logger : [logging object] instance of Python logging module
    
    Output: info in log (and can potentially delete some files that are deemed irrelevant)
    """
    if staticmaps_check != None:
        for temp_file in staticmaps_check:
            if not os.path.exists(os.path.join(out_clone_dir, temp_file + '.map')):
                logger.warning(temp_file + ' map from setup ini file was not created!')
        # remove staticmaps that are no longer in use / not relevant
        logger.info('Cleaning up staticmaps folder...')
        for temp_file in os.listdir(out_clone_dir):
            if temp_file != 'clim' and temp_file.split('.map')[0] not in staticmaps_check:
                temp_path = os.path.join(out_clone_dir, temp_file)
                logger.warning("Deleting unexpected/outdated file " + temp_path)
                os.remove(temp_path)


def checkIntbl(intbls_check, out_intbl_dir, logger):
    """
    Checks if all intbl's specified in setup ini file were created.
    
    Parameters:
        intbls_check : [list|string] files that should be present in intbl folder
        out_intbl_dir : [string] path to intbl folder
        logger : [logging object] instance of Python logging module
    
    Output: info in log
    """
    if intbls_check != None:
        for temp_file in intbls_check:
            if not os.path.exists(os.path.join(out_intbl_dir, temp_file + '.tbl')):
                logger.error(temp_file + ' intbl from setup ini file not found as template intbl! Could not copy!')


def relocateStation(lat, lon, area, uparea_dataset, nodata=-9999, point_buffer_val=0.04, point_buffer_style=3, area_margin=0.5, digits_new_coords=4):
    """
    Relocates gauging station based on upstream area.
    
    Parameters:
        lat : [float] latitude of gauging station location
        lon : [float] longitude of gauging station location
        uparea_dataset : [array] 2D array of upstream area (e.g. file loaded with rasterio)
        nodata : [int|float] no data / missing value in upstream area dataset
        point_buffer_style : [int] type of buffer around point (1=circle [default in function], 3=square [default here])
        point_buffer_val : [float] search radius (buffer around lat/lon point)
        area_margin : [float] margin for matching upstream area (0.5 = between 50% lower or higher)
        digits_new_coords : [int] number of digits to store for output, i.e. new lat/lon coordinate
    
    Output: relocated gauging station (new lat/lon coordinates)
    """
    # get data around current location
    data = rasterio.mask.mask(uparea_dataset, [Point(lon, lat).buffer(point_buffer_val, cap_style=point_buffer_style)], crop=True, all_touched=True)
    # get data mask
    mask_missings = data[0][0]==nodata#map_meta['nodata']
    mask_margin   = (data[0][0]>=area*(1-area_margin)) & (data[0][0]<=area*(1+area_margin))
    mask_final    = mask_missings | ~mask_margin
    # check if location can be relocated
    if not mask_final.all():
        # mask data
        data_masked = np.ma.masked_where(mask_final, data[0][0])
        # get distances from current location
        x_dist        = np.arange(data_masked.shape[1])-data_masked.shape[1]/2+np.mod(data_masked.shape[1],2)/2
        y_dist        = data_masked.shape[0]/2-np.mod(data_masked.shape[0],2)/2-np.arange(data_masked.shape[0])
        x_dist,y_dist = np.meshgrid(x_dist,y_dist)
        dists         = np.rint(np.sqrt(x_dist**2+y_dist**2)).astype(np.int)
        # find location that fullfills criteria of exceedence and has minimum distance to source
        ind_lat, ind_lon = np.unravel_index(np.ma.masked_where(mask_final, dists).argmin(), data[0][0].shape)
        # construct coordinates of new location from indices
        new_lat = np.round(data[1][5] + (ind_lat + 0.5) * data[1][4], digits_new_coords)
        new_lon = np.round(data[1][2] + (ind_lon + 0.5) * data[1][0], digits_new_coords)
        return Point(new_lon, new_lat)
    else:
        return np.NaN


def getabspath(path, root):
    "return absolute path if relative path given"
    if not os.path.isabs(path):
        path = os.path.normpath(os.path.join(root, path))
    return path


def makeSingleIniClones(ini_file):
    """
    Does all setup processing for a single ini file. This can be for a single catchment, or for all catchments in a certain folder.
    
    Parameters:
        ini_file : [string] ini file with relevant configuration
    
    Output: A lot of files created, deleted and/or changed. Returns logging statistics to be used in top level log file.
    """
    
    # read the ini-file
    root = os.path.dirname(os.path.abspath(ini_file))
    script_root = os.path.dirname(os.path.realpath(__file__))
    config = cp.ConfigParser()
    try:
        config.read(ini_file)
    except:
        sys.exit("ERROR: Not possible to open 'ini'- file.")
    
    if config.has_section("STRUCTURE"):
        directory_topo     = getabspath(check_key(config["STRUCTURE"],"input_topo", r"p:/wflow_global/static_data/base/hydro_merit"), root)
        directory_in       = getabspath(check_key(config["STRUCTURE"],"input_other", r"p:/wflow_global"), root)
        directory_out      = getabspath(check_key(config["STRUCTURE"],"output", "output"), root)
        paramsfolder       = check_key(config["STRUCTURE"],"parameters", "static_data")
        settingsfile       = check_key(config["STRUCTURE"],"file", "settings_scaling.csv")
        intblfolder        = getabspath(check_key(config["STRUCTURE"],"intbl", r"p:/wflow_global/static_data/wflow_sbm_parameters/intbl_template"), root)
        reservoirs_path    = check_key(config["STRUCTURE"],"reservoirs", r"p:/wflow_global/static_data/base/waterbodies/reservoir-db.gpkg")
        lakes_path         = check_key(config["STRUCTURE"],"lakes", r"p:/wflow_global/static_data/base/waterbodies/lake-db.gpkg")
        catchmentsfolder   = check_key(config["STRUCTURE"],"catchments", r"data/catchments/catchments.geojson")
        riversfolder       = check_key(config["STRUCTURE"],"rivers", r"data/rivers/rivers.geojson")
        modelbuilder_path  = check_key(config["STRUCTURE"],"modelbuilder", "modelbuilder")
        path               = check_key(config["STRUCTURE"],"path", "staticmaps")
        clonefile          = check_key(config["STRUCTURE"],"clone", "wflow_dem.map")
        clonefolder        = check_key(config["STRUCTURE"],"clonefolder", "*/") 
        discharges_path    = getabspath(check_key(config["STRUCTURE"],"discharges", r"p:/wflow_global/static_data/mean_discharge_1k/FLO1K.ts.1960.2015.qav.nc"), root)
        path_grdc_stations = check_key(config["STRUCTURE"],"path_grdc_stations", r"p:/wflow_global/static_data/gauging_stations/grdc_stations.xlsx")
        model_ini          = check_key(config["STRUCTURE"],"model_ini", "wflow_sbm")
        setup_folder       = check_key(config["STRUCTURE"],"setup_info", "setup")
    if config.has_section("CONTROLS"):
        do_modelbuilder     = bool(int(check_key(config["CONTROLS"],"do_modelbuilder", 0)))
        use_current_rivers  = bool(int(check_key(config["CONTROLS"],"use_current_rivers", 0)))
        use_merit_derived   = bool(int(check_key(config["CONTROLS"],"use_merit_derived", 1)))
        use_pyflwdir_point  = bool(int(check_key(config["CONTROLS"],"use_pyflwdir_point", 0)))
        upstream_from_point = bool(int(check_key(config["CONTROLS"],"upstream_from_point", 0)))
        min_stream_order    = int(int(check_key(config["CONTROLS"],"min_stream_order", 6)))
        get_custom_widths   = bool(int(check_key(config["CONTROLS"],"get_custom_widths", 1)))
        do_lakes            = bool(int(check_key(config["CONTROLS"],"do_lakes", 0)))
        do_reservoirs       = bool(int(check_key(config["CONTROLS"],"do_reservoirs", 0)))
        debug_discharge     = bool(int(check_key(config["CONTROLS"],"debug_discharge", 1)))
        template_ini        = bool(int(check_key(config["CONTROLS"],"template_ini", 1)))
        save_pyflwdir       = bool(int(check_key(config["CONTROLS"],"save_pyflwdir", 0)))
        get_pyflwdir_riv    = bool(int(check_key(config["CONTROLS"],"get_pyflwdir_riv", 0)))
        interp_soilthick    = bool(int(check_key(config["CONTROLS"],"interp_soilthick", 1)))
        get_grdc_gauges     = bool(int(check_key(config["CONTROLS"],"grdc_gauges", 1)))
    if config.has_section("PARS"):
        resolution       = float(check_key(config["PARS"],"resolution", 0.008333333333333333))
        alpha            = check_key(config["PARS"],"alpha", 60)
        M_method         = int(check_key(config["PARS"],"M_method", 2))
        M_minmax         = int(check_key(config["PARS"],"M_minmax", 100000))
        riv_upa          = float(check_key(config["PARS"],"pyflwdir_riv_upa", 30.))
        smooth_len       = float(check_key(config["PARS"],"pyflwdir_smooth_len", 1e4))
        ucat_ratio       = int(check_key(config["PARS"],"pyflwdir_ucat_ratio", 10))
        res_min_area     = float(check_key(config["PARS"],"res_min_area_km2", 0))
        lake_min_area    = float(check_key(config["PARS"],"lake_min_area_km2", 3))
        res_intbl_method = int(check_key(config["PARS"],"res_intbl_method", 1))
        res_minfrac_min  = float(check_key(config["PARS"],"res_minfrac_min", 0.0))
        res_minfrac_max  = float(check_key(config["PARS"],"res_minfrac_max", 0.9))
        res_fullfrac_min = float(check_key(config["PARS"],"res_fullfrac_min", 0.1))
        res_fullfrac_max = float(check_key(config["PARS"],"res_fullfrac_max", 1.0))
    if config.has_section("FILES"):
        intbls_check     = check_key(config["FILES"],"tbls", None)
        intbls_check     = intbls_check.split('\n')
        staticmaps_check = check_key(config["FILES"],"maps", None)
        staticmaps_check = staticmaps_check.split('\n')
    if config.has_section("TESTING"):
        tests = check_key(config["TESTING"],"do_test", None)
    else:
        tests = None

    # tests
    if tests != None:
        tests = tests.split(',')
    else:
        tests = [tests]
    
    # paths
    clones = glob.glob(os.path.join(directory_out, clonefolder))
    parameters_path = os.path.join(directory_in, paramsfolder)
    #settings_path = os.path.join(directory_out, settingsfile)
    settings_path = os.path.join(script_root, settingsfile)

    # variable for storing information on each individual catchment
    log_stats = []
    
    # quick checks before going into loop
    if use_merit_derived and do_modelbuilder:
        sys.exit("ERROR: Running modelbuilder while using MERIT derived data! This will cause MERIT derived data to be overwritten with modelbuilder results, which makes it impossible to use MERIT derived data! Please choose one of the two, and note that it is advised to use MERIT derived data when possible.")
    if get_custom_widths and not use_merit_derived:
        sys.exit("ERROR: Custom river widths can only be obtained when using MERIT derived data! Please adjust setup ini file!")
    if get_grdc_gauges and not use_merit_derived:
        sys.exit("ERROR: GRDC stations for gauges map can only be used when also using MERIT derived data! Please adjust setup ini file!")
    if use_merit_derived:
        if not os.path.exists(directory_topo):
            sys.exit("ERROR: Path to topographic base data does not exist! Check setup ini file and/or topographic input directory!")
    if not os.path.exists(parameters_path):
        sys.exit("ERROR: Path to parameter base data does not exist! Check setup ini file and/or input directory!")
    if len(clones) == 0:
        sys.exit("ERROR: Folder(s) where model should be set up do not exist! Check setup ini file and/or output directory!")
    
    # loop over each catchment
    for folder in clones:
        
        # catchment
        basin_id = -1
        catchment_id    = os.path.basename(os.path.normpath(folder))
        catchments_path = os.path.join(folder, catchmentsfolder)
        if not os.path.exists(os.path.split(catchments_path)[0]):
            os.makedirs(os.path.split(catchments_path)[0], exist_ok=True)
        model_ini_path  = os.path.join(folder, model_ini + '.ini')
        
        # create setup/debug folder
        setup_path = os.path.join(folder, setup_folder)
        if not os.path.exists(setup_path):
            os.makedirs(setup_path, exist_ok=True)
        else:
            # clear folder if it already existed
            try:
                for temp_file in os.listdir(setup_path):
                    os.remove(os.path.join(setup_path, temp_file))
            except OSError:
                pass
        
        # get logger
        logger = setup_logging.setupLogging(setup_path)
        logger.info('Starting setup procedure of catchment/case ' + catchment_id + ' (' + ini_file + ')')
        logger.info('Running modelbuilder:           ' + str(do_modelbuilder))
        logger.info('Using MERIT derived data:       ' + str(use_merit_derived))
        logger.info('Deriving custom river widths:   ' + str(get_custom_widths))
        logger.info('Using GRDC stations for gauges: ' + str(get_grdc_gauges))
        logger.info('Adding reservoirs:              ' + str(do_reservoirs))
        logger.info('Adding lakes:                   ' + str(do_lakes))
        
        if not None in tests:
            logger.warning("  <<<  RUNNING SCRIPTS IN TEST MODE!  >>>  ")
            if len(tests) > 1:
                logger.info('Test parameters:')
                for test_subject in tests:
                    logger.info(' - ' + test_subject)
            else:
                logger.info('Test parameter: ' + tests[0])
        
        # modelbuilder
        if do_modelbuilder:
            # check potential input/paths
            logger.info('Starting modelbuilder process...')
            if not os.path.exists(modelbuilder_path):
                sys.exit('ERROR: Path to modelbuilder could not be found! Check path in ini file!')
            if use_current_rivers:
                rivers_path = os.path.join(folder, riversfolder)
                logger.info('Using existing rivers from ' + rivers_path)
                if not os.path.exists(rivers_path):
                    logger.fatal('Could not find existing rivers!')
                    sys.exit()
            else:
                rivers_path = None
            # run modelbuilder
            runModelbuilder(folder, modelbuilder_path, catchments_path, rivers_path, logger)
        
        # check on paths which are vital for setup
        do_proceed = True
        out_clone_dir = os.path.join(folder, path)
        if not os.path.exists(out_clone_dir):
            try:
                os.mkdir(out_clone_dir)
            except:
                logger.fatal('Folder with/for maps does not exist and could not be created! Cannot set up maps! (' + out_clone_dir + ')')
                do_proceed = False
        if not os.path.exists(intblfolder):
            logger.fatal('intbl template folder does not exists! Cannot copy files! (' + intblfolder + ')')
            do_proceed = False
        
        # get hydro MERIT basin ID and bounding box (if this is to be used)
        if use_merit_derived:
            if not use_pyflwdir_point:
                xy2 = None
                # get basin ID
                if basin_id == -1:
                    if catchment_id.isnumeric():
                        logger.info('Reading hydro MERIT catchment ID from model folder name: ' + catchment_id)
                        basin_id = int(catchment_id)  # model folder is basin id
                    else:
                        try:
                            check_files = glob.glob(os.path.join(folder, 'data', 'catchments', '*')) + glob.glob(os.path.join(folder, '*.basinid'))
                            for temp_file in check_files:
                                basin_id_str = os.path.basename(temp_file).split('.')[0]
                                if basin_id_str.isnumeric():
                                    logger.info(f'Reading hydro MERIT catchment ID from file in catchments folder: {basin_id_str}')
                                    basin_id = int(basin_id_str)  # file in catchments folder is basin id
                                    break
                            
                        except:
                            pass
                # check if basin ID can be found in CSV files
                if basin_id != -1:
                    try:
                        # find relevant CSV file
                        pfaf_id            = basin_id // 10**7
                        fn_outlets         = os.path.join(directory_topo, f'pfaf{pfaf_id:02d}_outlets.csv')
                        merit_basin_lookup = pd.read_csv(fn_outlets, index_col='pfaf3')
                        # get basin bounding box listed in CSV
                        bbox               = merit_basin_lookup.loc[basin_id, ['xmin', 'ymin', 'xmax', 'ymax']].values
                    except:
                        basin_id = -1
                if basin_id == -1:
                    logger.fatal('Could not find correct hydro MERIT catchment ID! Cannot execute required setup for this catchment/case!')
                    do_proceed = False
            else:
                # get xy (lon,lat) point
                try:
                    for temp_file in os.listdir(os.path.join(folder, 'data', 'catchments')) + glob.glob(os.path.join(folder, '*.xy')):
                        if temp_file != 'catchments.geojson':
                            # check if file contains xy coords tuple
                            fn_temp_file = getabspath(temp_file, os.path.join(folder, 'data', 'catchments'))
                            with open(fn_temp_file, 'r') as f:
                                xy = tuple(float(s) for s in f.readline().strip().split(','))
                            logger.info(f'Getting hydro MERIT catchment from point (x: {xy[0]}, y:{xy[1]})')
                            # get basin bounding box and ID from point
                            bbox, basin_id, xy2 = get_merit_basin_bbox(xy, directory_topo, upstream_from_point=upstream_from_point, min_sto=min_stream_order)
                            if not upstream_from_point:
                                logger.info(f'Hydro-MERIT catchment ID at point (x: {xy[0]}, y:{xy[1]}): {basin_id}')
                                xy2 = None # this is important to not add a pit at the point location later
                            else:
                                logger.info(f'Point snapped to (x: {xy2[0]:.5f}, y:{xy2[1]:5f}) based on a minimum required stream order of {min_stream_order}')
                            break
                except:
                    logger.fatal('Could not find correct hydro MERIT catchment from point location! Cannot execute required setup for this catchment/case!')
                    do_proceed = False
        
        # reservoir and lakes ini; moved forward incase do_proceed == False
        res_catch_intersect     = None
        do_reservoirs_catchment = False
        lakes_catch_intersect   = None
        do_lakes_catchment      = False

        # start of actual processing
        if do_proceed:
            
            if use_merit_derived:
            
                # get topographic base data from hydro MERIT
                logger.info('Starting pyflwdir processing...')
                
                # get scale ratio from MERIT native resolution and desired model resolution
                res = 1/1200. # hydro MERIT data has 3 arcsec resolution
                scale_ratio = resolution / res
                if scale_ratio.is_integer():
                    scale_ratio = int(scale_ratio)
                    logger.info('Scale ratio (model/MERIT) is ' + str(scale_ratio))
                else:
                    logger.warning('Scale ratio (model/MERIT) is not an integer! (' + str(scale_ratio) + ')')
                    scale_ratio = round(scale_ratio)
                    logger.warning('Rounding scale ratio to ' + str(scale_ratio))
                    resolution = res * scale_ratio
                    logger.warning('Model resolution changed to ' + str(resolution))
                
                # upscale flow direction and get subgrid basins (and update model resolution if required)
                logger.info('Upscaling hydro MERIT data...')
                upscale_merit_basin(scale_ratio, bbox, basin_id, directory_topo, folder, xy=xy2, logger=logger)
                
                # use flow direction network to derive relevant topographic maps
                logger.info('Obtaining topographic maps from hydro MERIT data...')
                network_merit_basin(folder, smooth_len=smooth_len, ucat_ratio=ucat_ratio, riv_shape=get_pyflwdir_riv, basin_shape=True, logger=logger)
                
                # perform simple resampling to model resolution of slope and elevation maps
                logger.info('Resampling slope and elevation maps to model resolution...')
                resample_merit_basin(directory_topo, folder)
                
                # write wflow topographic static maps based on just created GTiff maps
                logger.info('Writing wflow topographic staticmaps...')
                wflow_topomaps(folder, riv_upa=riv_upa, logger=logger)
                
                # move catchments geojson to its expected locations
                logger.info('Copying catchments.geojson derived from hydro MERIT to ' + catchments_path)
                shutil.copy(os.path.join(folder, 'catchments.geojson'), catchments_path)
                os.remove(os.path.join(folder, 'catchments.geojson'))
                
                # move other data created with pyflwdir to dedicated folder (or remove if specified in setup ini)
                try:
                    shutil.rmtree(os.path.join(folder, 'pyflwdir'))
                    logger.info('Clearing any previously stored pyflwdir files...')
                except:
                    logger.warning('Could not clear previously stored pyflwdir files!')
                if save_pyflwdir:
                    logger.info('Moving data created with pyflwdir to ' + os.path.join(folder, 'pyflwdir'))
                    try:
                        os.mkdir(os.path.join(folder, 'pyflwdir'))
                    except:
                        pass
                else:
                    logger.info('Removing all other data now created with pyflwdir...')
                for file_or_dir in os.listdir(folder):
                    if os.path.isfile(os.path.join(folder, file_or_dir)):
                        if not file_or_dir.split('.')[-1] in ['txt', 'ini']:
                            if save_pyflwdir:
                                try:
                                    shutil.copy(os.path.join(folder, file_or_dir), os.path.join(folder, 'pyflwdir', file_or_dir))
                                except:
                                    logger.warning('Could not copy ' + file_or_dir)
                            try:
                                os.remove(os.path.join(folder, file_or_dir))
                            except:
                                logger.warning('Could not delete ' + file_or_dir)
                    else:
                        if file_or_dir in ['flwdir', 'river', 'upadff_1perc']:
                            if save_pyflwdir:
                                try:
                                    os.mkdir(os.path.join(folder, 'pyflwdir', file_or_dir))
                                except:
                                    pass
                                for file_or_dir_2 in os.listdir(os.path.join(folder, file_or_dir)):
                                    if os.path.isfile(os.path.join(folder, file_or_dir, file_or_dir_2)):
                                        try:
                                            shutil.copy(os.path.join(folder, file_or_dir, file_or_dir_2), os.path.join(folder, 'pyflwdir', file_or_dir, file_or_dir_2))
                                        except:
                                            logger.warning('Could not copy ' + file_or_dir_2)
                            try:
                                shutil.rmtree(os.path.join(folder, file_or_dir))
                            except:
                                logger.warning('Could not delete ' + file_or_dir)
                logger.info('Finished pyflwdir processing.')
                
                if get_custom_widths:
                    logger.info('Overwriting MERIT Hydro riverwidths with custom riverwidths...')
                    setup_river_widths.createWidthsMap(out_clone_dir, path_grdc_stations, os.path.join(parameters_path, 'riverwidth'), debug_path=setup_path, logger=logger)
            
            # read catchment geojson file, check if valid, fix if not
            catchment = gpd.read_file(catchments_path)
            if catchment.geometry.is_valid.all():
                pass
            else:
                logger.warning('Catchment is not valid, probably due to self-intersection point(s)! Fixing...')
                logger.info('Copying current catchment file to ' + setup_path)
                shutil.copy(catchments_path, os.path.join(setup_path, 'catchments_original.geojson'))
                # changing geometry type, saving to new file
                changeGeoJSONgeomType(setup_path, logger=logger)
                # read in new catchments geojson file and check if valid
                catchment = gpd.read_file(os.path.join(setup_path, 'catchments_v2.geojson'))
                if catchment.geometry.is_valid.all():
                    logger.info('Catchment is now valid. Converting from MultiLineStrings back to MultiPolygon(s).')
                    for catch_i in range(len(catchment.geometry)):
                        temp_polys = []
                        for catch_j in range(len(catchment.geometry[catch_i])):
                            temp_coords = catchment.geometry[catch_i][catch_j].coords
                            if temp_coords[0] == temp_coords[-1]:
                                temp_polys.append(Polygon(temp_coords))
                            else:
                                logger.error('First and last coordinate are not the same for feature '+ str((catch_i+1)*(catch_j+1)) + '! Cannot create valid closed polygon!')
                        temp_geom = MultiPolygon(temp_polys)
                        catchment.geometry[catch_i] = temp_geom
                else:
                    logger.error('Catchment still not valid! Proceeding with setup but good results cannot be guaranteerd! Check lakes and reservoirs when finished!')
            
            # get reservoirs of sufficient size within catchment
            if do_reservoirs:
                reservoirs_check        = checkWaterbodies('reservoir', catchment, reservoirs_path, min_area=res_min_area, logger=logger)
                res_catch_intersect     = reservoirs_check['data']
                do_reservoirs_catchment = reservoirs_check['do_process']
            
            # get lakes of sufficient size within catchment
            if do_lakes:
                lakes_check           = checkWaterbodies('lake', catchment, lakes_path, min_area=lake_min_area, logger=logger)
                lakes_catch_intersect = lakes_check['data']
                do_lakes_catchment    = lakes_check['do_process']
            
            # gauging stations
            if get_grdc_gauges:
                logger.info('Using GRDC stations to construct wflow_gauges.map file...')
                if use_merit_derived:
                    uparea_src = rasterio.open(os.path.join(out_clone_dir, 'wflow_uparea.map'), 'r')
                    #uparea_data = uparea_src.read(1)
                    uparea_meta = uparea_src.meta
                    # read GRDC stations as DataFrame from Excel file
                    logger.info('Reading GRDC stations from ' + path_grdc_stations)
                    gauging_stations = pd.read_excel(path_grdc_stations, na_values='n.a.')
                    logger.info("Removing stations with no valid 'area' value...")
                    gauging_stations = gauging_stations.loc[~gauging_stations['area'].isnull()]
                    gauging_stations = gauging_stations.loc[gauging_stations['area']>0]
                    # convert to GeoDataFrame
                    gauging_stations = gpd.GeoDataFrame(gauging_stations, geometry=[Point(x, y) for x, y in zip(gauging_stations.long, gauging_stations.lat)])
                    # get only stations within catchment
                    gauging_stations = gauging_stations.loc[~gpd.tools.sjoin(gauging_stations, catchment, how='left')['index_right'].isnull()]
                    logger.info('Number of potentially valid GRDC stations found within catchment: ' + str(len(gauging_stations)))
                    # check if there are any stations to proceed
                    if not gauging_stations.empty:
                        # relocate stations
                        logger.info('Relocating GRDC stations so that they are located on model river cells...')
                        #relocateStation(lat, lon, area, uparea_dataset, nodata=-9999, point_buffer_val=0.04, point_buffer_style=3, area_margin=0.5, digits_new_coords=4)
                        gauging_stations['geometry'] = gauging_stations.apply(lambda row, uparea_src=uparea_src, uparea_meta=uparea_meta: relocateStation(row['lat'], row['long'], row['area'], uparea_src, uparea_meta['nodata']), axis=1)
                        # remove stations that could not be relocated (i.e. no upstream area within allowed margin found within window around base location)
                        gauging_stations = gauging_stations.loc[~gauging_stations['geometry'].isnull()]
                        logger.info('Number of GRDC stations left after relocating: ' + str(len(gauging_stations)))
                        # remove stations located in waterbodies
                        logger.warning('Step to remove GRDC stations located within lakes or reservoirs has not been implemented yet! Check if wflow_gauges.map is correct!')
                        # convert to map
                        out_arr = pcr.pcr2numpy(pcr.readmap(os.path.join(out_clone_dir, 'wflow_subcatch.map')), uparea_meta['nodata'])
                        out_arr = (out_arr/out_arr)-1 #make sure default array contains zeros only
                        shapes = ((geom,value) for geom, value in zip(gauging_stations.geometry, gauging_stations['grdc_no']))
                        gauges_arr = rasterio.features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=uparea_meta['transform'], all_touched=False)
                        gauges_arr[gauges_arr==0] = uparea_meta['nodata'] # prevent zeros in output map
                        createMapFile(gauges_arr, out_clone_dir, 'wflow_gauges', uparea_meta, logger)
                        logger.info('New wflow_gauges.map created from GRDC stations.')
                    else:
                        logger.warning('No valid GRDC stations found within catchment, skipping this procedure!')
                else:
                    logger.warning('Creation of wflow_gauges.map from GRDC stations only implemented with hydroMERIT! Cannot be executed with current setup!')
            
            # catchment/model ini file (template)
            if template_ini:
                if model_ini.endswith('.ini'):
                    template_ini_path = getabspath(model_ini, root)
                else:
                    template_ini_path = os.path.join(script_root, model_ini+'.ini')
                logger.info('Overwriting model ini file with template from ' + template_ini_path)
                if os.path.exists(template_ini_path):
                    if os.path.exists(model_ini_path):
                        os.remove(model_ini_path)
                    shutil.copy(template_ini_path, model_ini_path)
                else:
                    logger.error('Could not find specified template ini file!')
            
            # staticmaps
            if not 'no_generic' in tests:
                logger.info('Setting up maps in ' + out_clone_dir)
                make_clone(resolution, settings_path, interp_soilthick, M_method, M_minmax, parameters_path, out_clone_dir, clonefile, logger)
            else:
                logger.info('TEST MODE: skipping creation of new generic map files.')
            if not use_merit_derived:
                if not 'no_slope' in tests:
                    get_slope(path_geojson=os.path.join(modelbuilder_path, 'settings.json'), path_model_abs=folder, path_DEM_rel=os.path.join('data', 'dem'), dst_res=0.005, region_filter="region", do_delete_temp=True, logger=logger)
                else:
                    logger.info('TEST MODE: skipping creation of new upscaled slope map.')
            # check if all staticmaps specified in setup ini file were created
            if not 'no_generic' in tests:
                checkStaticmaps(staticmaps_check, out_clone_dir, logger)
                logger.info('Finished maps setup')
            
            # update riverwidth map using all lakes and reservoirs within catchment (no restrictions on size)
            if use_merit_derived and not get_custom_widths and not 'no_rivwidth' in tests:
                reservoirs_all = getWaterbodiesCatchment(catchment, reservoirs_path)
                lakes_all      = getWaterbodiesCatchment(catchment, lakes_path)
                updateRiverWidths(out_clone_dir, reservoirs_all, lakes_all, flag=-2, logger=logger)
            
            # intbl's
            out_intbl_dir = os.path.join(folder, 'intbl')
            if not os.path.exists(out_intbl_dir):
                os.makedirs(out_intbl_dir, exist_ok=True)
            if not 'no_intbl' in tests:
                logger.info("Clearing contents of intbl folder...")
                if os.path.exists(out_intbl_dir):
                    for temp_file in os.listdir(out_intbl_dir):
                        os.remove(os.path.join(out_intbl_dir, temp_file))
                logger.info("Copying template intbl's to intbl folder...")
                for temp_file in os.listdir(intblfolder):
                    if temp_file.split('.tbl')[0] not in intbls_check:
                        logger.warning('Template intbl copied to intbl folder was not listed in setup ini file: ' + temp_file.split('.tbl')[0])
                    shutil.copy(os.path.join(intblfolder, temp_file), os.path.join(out_intbl_dir, temp_file))
                # check if all intbl specified in setup ini file were created
                checkIntbl(intbls_check, out_intbl_dir, logger)
                logger.info('Finished intbl setup')
            else:
                logger.info("TEST MODE: skipping creation of new intbl's.")
            
            # reservoirs/lakes
            processWaterbodies('reservoir', do_reservoirs_catchment, res_catch_intersect, out_intbl_dir, out_clone_dir, model_ini_path, [res_minfrac_min,res_minfrac_max], [res_fullfrac_min,res_fullfrac_max], res_intbl_method, setup_path, logger)
            processWaterbodies('lake', do_lakes_catchment, lakes_catch_intersect, out_intbl_dir, out_clone_dir, model_ini_path, None, None, None, setup_path, logger)
            
            # catchment/model ini file
            if not 'no_discharge' in tests:
                logger.info('Calculating AnnualDischarge...')
                if os.path.exists(discharges_path):
                    annualDischarge = flo1k.getAnnualDischarge(discharges_path, catchment, debug_discharge, debug_path=setup_path, logger=logger)
                else:
                    logger.error('Path to discharge dataset is invalid! Could not calculate AnnualDischarge!')
            else:
                annualDischarge = 2290 # default value for testing, can be used to suppress slow calculation from data
                logger.info('TEST MODE: skipping calculation of AnnualDischarge, using pre-set default value instead.')
            logger.info('AnnualDischarge is ' + str(annualDischarge) + ' m3/s')
            replaceLinesIniFile(['AnnualDischarge', 'Alpha'], [annualDischarge, alpha], model_ini_path, logger=logger)
            
            # clear outdated folders
            logger.info('Clearing outdated folders...')
            outdated_folders = ['outstate']
            for outdated_folder in outdated_folders:
                outdated_path = os.path.join(folder, outdated_folder)
                if os.path.exists(outdated_path) and os.path.isdir(outdated_path):
                    logger.info('Clearing ' + outdated_path)
                    try:
                        for outdated_file in os.listdir(outdated_path):
                            os.remove(os.path.join(outdated_path, outdated_file))
                    except:
                        logger.error('Could not clear file(s) from ' + outdated_path)
        
        # copy relevant files to setup/debug folder
        logger.info('Copying files (settings and scripts) used in setup procedure to ' + setup_path)
        list_scripts = [sys.argv[0], ini_file, 'catchment_FLO1K.py', 'setup_logging.py']
        if not use_merit_derived:
            list_scripts.extend(['upscaled_slope.py'])
        if (do_proceed) and (do_reservoirs_catchment or do_lakes_catchment):
            list_scripts.extend(['waterbodies.py', 'wflow_lake_intbl.py', 'wflow_reservoir_intbl.py', 'reservoir_intbl_utils.py'])
        copyFiles(setup_path, [settings_path], script_root, list_scripts, logger)
        
        # ReadMe file
        logger.info('Writing ReadMe file in ' + folder)
        createReadMe(folder, setup_folder)
        
        # check log stats
        logger_stats = setup_logging.getStats(logger)
        
        # add log stats with formatting
        log_stats.append("Catchment ID: " + catchment_id + "\n")
        log_stats.append("Warnings:     " + str(logger_stats['warn']) + "\n")
        log_stats.append("Errors:       " + str(logger_stats['err']) + "\n")
        log_stats.append("Criticals:    " + str(logger_stats['crit']) + "\n\n")
        
        # end of setup for this catchment/case
        logger.info('Finished setup procedure of ' + catchment_id + ' (' + ini_file + ')')
        setup_logging.closeLogging()
        print('')
        
    return log_stats


def make_clones(ini_files):
    # initialize top level logging variable
    top_lvl_log = ["This file contains collected logging statistics of each individual catchment.\n",
                   "Detailed log files for each catchment can be found in their respective directories.\n\n",
                   "Warning  = Potential issue encountered, but a workaround was used instead. Good to check.\n",
                   "Error    = Potential issue encountered without any workarounds. Should be checked.\n",
                   "Critical = Critical issue(s) that prevented processing of this catchment.\n\n",
                   "Timestamp:    "+ datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n\n"]
    # run setup function for each ini file, returning top level logging information for each
    print('')
    if len(ini_files) == 1:
        top_lvl_log.append("Ini file:     " + ini_files[0] + "\n\n")
        top_lvl_log.append(makeSingleIniClones(ini_files[0]))
    elif len(ini_files) > 1:
        for ini_file in ini_files:
            top_lvl_log.append("Ini file:     " + ini_file + "\n\n")
            top_lvl_log.append(makeSingleIniClones(ini_file))
    else:
        sys.exit('ERROR: Please provide an ini file argument!')
    # write top level logging information to file
    log_path = os.path.join(os.getcwd(), 'LOG.txt')
    print('Writing top level log file to ' + log_path)
    with open(log_path, 'w') as log:
        for line in top_lvl_log:
            if isinstance(line, str):
                log.write(line)
            else:
                for line2 in line:
                    log.write(line2)

if __name__ == "__main__":
    make_clones(sys.argv[1:])