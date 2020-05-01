# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 13:12:00 2019

@author: haag

Includes parts adapted from script by Verseveld (2018) for project 11202412-007,
which in turn was based on https://gis.stackexchange.com/questions/289775/python-mask-netcdf-data-using-shapefile-xarray-geopandas

This script extract discharge information for a specific catchment from the FLO1K dataset.

The discharge information is used to set the AnnualDischarge parameter in a Wflow model's ini file.
The FLO1K dataset was created by Barbarossa et al. (2018): https://www.nature.com/articles/sdata201852
"""

import os
import math
import geojson
import xarray as xr
import numpy as np
import geopandas as gpd

from rasterio import features
from affine import Affine

import setup_logging

import matplotlib.pyplot as plt

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()


def transform_from_latlon(lat, lon):
    """
    Input 1D array of lat / lon and output an Affine transformation.
    Copied from Verseveld (2018).
    """
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    trans = Affine.translation(lon[0], lat[0])
    scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
    return trans * scale


def rasterize(shapes, coords, latitude='lat', longitude='lon', fill=np.nan, **kwargs):
    """
    Rasterize a list of (geometry, fill_value) tuples onto the given xray coordinates. This only works for 1d latitude and longitude arrays.
    Copied from Verseveld (2018).
    
    usage:
    -----
    1. read shapefile to geopandas.GeoDataFrame
          `states = gpd.read_file(shp_dir+shp_file)`
    2. encode the different shapefiles that capture those lat-lons as different
        numbers i.e. 0.0, 1.0 ... and otherwise np.nan
          `shapes = (zip(states.geometry, range(len(states))))`
    3. Assign this to a new coord in your original xarray.DataArray
          `ds['states'] = rasterize(shapes, ds.coords, longitude='X', latitude='Y')`
    
    arguments:
    ---------
    : **kwargs (dict): passed to `rasterio.rasterize` function
    
    attrs:
    -----
    :transform (affine.Affine): how to translate from latlon to ...?
    :raster (numpy.ndarray): use rasterio.features.rasterize fill the values
      outside the .shp file with np.nan
    :spatial_coords (dict): dictionary of {"X":xr.DataArray, "Y":xr.DataArray()}
      with "X", "Y" as keys, and xr.DataArray as values
    
    returns:
    -------
    :(xr.DataArray): DataArray with `values` of nan for points outside shapefile
      and coords `Y` = latitude, 'X' = longitude.
    """
    transform = transform_from_latlon(coords[latitude], coords[longitude])
    out_shape = (len(coords[latitude]), len(coords[longitude]))
    raster = features.rasterize(shapes, out_shape=out_shape,
                                fill=fill, transform=transform,
                                dtype=float, **kwargs)
    spatial_coords = {latitude: coords[latitude], longitude: coords[longitude]}
    
    return xr.DataArray(raster, coords=spatial_coords, dims=(latitude, longitude))


def add_shape_coord_from_data_array(xr_da, catchment, coord_name):
    """
    Create a new coord for the xr_da indicating whether or not it is inside the shapefile.
    Creates a new coord - "coord_name" which will have integer values used to subset xr_da for plotting/analysis.
    
    Adapted from Verseveld (2018).
    
    Usage:
    -----
    precip_da = add_shape_coord_from_data_array(precip_da, "awash.shp", "awash")
    awash_da = precip_da.where(precip_da.awash==0, other=np.nan) 
    """
    # create a list of tuples (shapely.geometry, id), this allows for many different polygons within a .shp file (e.g. States of US)
    shapes = [(shape, n) for n, shape in enumerate(catchment.geometry)]
    
    # create a new coord in the xr_da which will be set to the id in `shapes`
    xr_da[coord_name] = rasterize(shapes, xr_da.coords, longitude='lon', latitude='lat')
    
    return xr_da


def extract_catchment_bounds(ds, catchment):
    """
    Extracts a subset of an dataset using rectangle catchment bounds with lat/lon coordinates.
    Dataset is expected to be loaded with xarray, catchment with geopandas.
    """
    if len(catchment.geometry) > 1:
        sys.exit('Catchment contains multiple polygons! This is not properly handled yet in script!')
    catchment_bounds = catchment.geometry[0].bounds
    return ds.sel(lat=slice(catchment_bounds[1],catchment_bounds[3]), lon=slice(catchment_bounds[0],catchment_bounds[2]))


def haversine(origin, destination):
    """
    Calculates distance between latitude longitude points, using the Haversine formula (see https://en.wikipedia.org/wiki/Haversine_formula)
    
    Author: Wayne Dyck (https://gist.github.com/rochacbruno/2883505)
    """
    lat1, lon1 = origin
    lat2, lon2 = destination
    radius = 6371 # km
    
    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c
    
    return d


def get_distances(list_points):
    """
    Calculates the distance matrix between a set of points.
    Adapted from https://stackoverflow.com/questions/22081503/distance-matrix-creation-using-nparray-with-pdist-and-squareform
    """
    # create empty distance array
    N = len(list_points)
    distance_matrix = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            distance_matrix[i, j] = haversine(list_points[i], list_points[j])
            distance_matrix[j, i] = distance_matrix[i, j]
    return distance_matrix


def writeToGeoJSON(points, distances, fc_props, path, filename='MaximumDischargePoints'):
    """
    Writes lat/lon points identified as location of maximum discharge to a geojson file.
    This allows for checking/testing after model has been setup, and stores this for later use.
    
    Parameters:
        points:     [list] list of lat/lon points of maximum discharge
        distances:  [array] distance matrix, output of get_distances function
        path:       [string] path to directory to place geojson file
        filename:   [string] name of to-be-constructed geojson file (without extension)
    
    output: geojson file
    """
    # construct geojson
    features = []
    for i in range(len(points)):
        props = {'id': i, 'distances': distances[i].round(2).tolist()}
        point = points[i]
        features.append(geojson.Feature(geometry=geojson.Point(point), properties=props))
    # create full path to file
    full_path = os.path.join(path, filename+'.geojson')
    # check if file already exists, and if so, remove
    if os.path.exists(full_path):
        os.remove(full_path)
    # write new file
    with open(full_path, 'w') as f:
        geojson.dump(geojson.FeatureCollection(features, properties=fc_props), f, sort_keys=True, ensure_ascii=False)


def writeToFigureTimeSeries(data, path, filename='FLO1K_catchment_max_discharge_time_series', extension='png'):
    """
    Writes a xarray time series dataset to a figure.
    
    Parameters:
        data:       [xarray dataset] xarray 1D array with time coordinates
        path:       [string] path to directory to place figure
        filename:   [string] name of to-be-constructed figure (without extension)
        extension:  [string] extension/type of figure
    
    output: png file
    """
    data.plot()
    plt.title('Yearly maximum discharge values over catchment')
    # create full path to file
    full_path = os.path.join(path, filename+'.'+extension)
    # check if file already exists, and if so, remove
    if os.path.exists(full_path):
        os.remove(full_path)
    # write new file
    plt.savefig(full_path, bbox_inches='tight')
    plt.close()


def writeToFigureMapPoints(data, points, path, filename='FLO1K_catchment_discharge_with_max_points', extension='png'):
    """
    Writes a xarray dataset with lat/lon coords to a figure, including the maximum discharge points, if any.
    
    Parameters:
        data:       [xarray dataset] xarray 2D array with spatial coordinates
        points:     [list] list of points with spatial coordinates
        path:       [string] path to directory to place figure
        filename:   [string] name of to-be-constructed figure (without extension)
        extension:  [string] extension/type of figure
    
    output: png file
    """
    data.plot.imshow()
    for point in points:
        plt.scatter(point[0], point[1], color='r')
    plt.title('Overall max discharge with yearly max points')
    # create full path to file
    full_path = os.path.join(path, filename+'.'+extension)
    # check if file already exists, and if so, remove
    if os.path.exists(full_path):
        os.remove(full_path)
    # write new file
    plt.savefig(full_path, bbox_inches='tight')
    plt.close()


def getAnnualDischarge(path_data, catchment, do_check=False, tolerance_abs=None, tolerance_frac=0.1, debug_path=None, logger=None):
    """
    Calculates AnnualDischarge (mean annual discharge over full time period) for use in Wflow ini file.
    Since the exact location of the river cells might not match between the Wflow model and the FLO1K dataset,
    it is assumed that the maximum discharge within a catchment represents the discharge we need (often the discharge at/near the outflow point).
    This assumption can be tested in this function by checking the location of the maximum discharge for each time step.
    If the locations are close together, it is probably correct, but if not, the assumption might not be valid.
    
    Parameters:
        path_data :     [string] path to FLO1K dataset (netCDF4 format)
        catchment:      [geopandas dataframe] catchment boundaries (geojson format) already loaded with geopandas
        do_check :      [boolean] controls whether to perform a check on the above mentioned assumption
        tolerance_abs:  [number] maximum allowed distance, in km, between points of maximum discharge, if not specified tolerance_frac is used instead
        tolerance_frac: [float] maximum allowed distance fraction (catchment bounds diagonal / maximum distance between points of maximum discharge)
        debug_path:     [string] full path to directory where setup/logging/debug files are stored
        logger          [logging object] instance of Python logging module
    
    output: [int] value for AnnualDischarge
    """
    # get dataset
    setup_logging.showInfo(logger, 'Opening discharge dataset and calculating mean yearly discharge...')
    ds = xr.open_dataset(path_data)
    # extract catchment bounds
    ds_bounds = extract_catchment_bounds(ds, catchment)
    # add variable with catchment polygon
    ds_bounds = add_shape_coord_from_data_array(ds_bounds, catchment, 'catchment')
    # mask data outside catchment
    ds_catchment = ds_bounds.where(ds_bounds.catchment==0, other=np.nan)
    # calculate mean annual discharge over full time period within catchment
    max_discharges = ds_catchment.qav.max(['lat', 'lon'])
    mean_discharge = float(np.mean(max_discharges))
    setup_logging.showInfo(logger, 'Finished calculating mean yearly discharge.')
    # check if maximum discharges are always from approximately the same location
    if do_check:
        setup_logging.showInfo(logger, 'Checking maximum discharge locations...')
        max_coords = []
        for time in ds_catchment.qav.time.values:
            temp_data = ds_catchment.qav.sel(time=time)
            max_loc   = temp_data.where(temp_data==temp_data.max(['lat','lon']), drop=True)
            max_lat   = max_loc.lat.values
            max_lon   = max_loc.lon.values
            if len(max_lat == 1) or len(max_lon == 1):
                if len(max_lat == 1):
                    for p in range(len(max_lon)):
                        temp_coord = (max_lon[p], max_lat[0])
                        if temp_coord not in max_coords:
                            max_coords.append(temp_coord)
                else:
                    for p in range(len(max_lat)):
                        temp_coord = (max_lon[0], max_lat[p])
                        if temp_coord not in max_coords:
                            max_coords.append(temp_coord)
            else:
                setup_logging.showWarning(logger, 'Found multiple latitude and longitude values for time ' + str(time) + '! Cannot combine to proper lat/lon points! Check on maximum discharge locations is not valid!')
        # calculate maximum distance between these points
        distances  = get_distances(max_coords)
        max_dist   = np.max(distances)
        setup_logging.showInfo(logger, 'Number of different points of maximum discharge found: ' + str(len(max_coords)))
        setup_logging.showInfo(logger, 'Maximum distance between points of maximum discharge is ' + str(np.round(max_dist,2)) + ' km')
        # check if maximum distance falls within allowed tolerance
        if tolerance_abs != None:
            setup_logging.showInfo(logger, 'Absolute tolerance for maximum distance is ' + str(tolerance_abs) + ' km')
            if max_dist > tolerance_abs:
                dist_diagonal = None
                dist_frac = None
                setup_logging.showWarning(logger, 'Maximum distance between points of maximum discharge exceeds allowed tolerance! Cannot guarantee that assumption used for calculation of AnnualDischarge is valid! Please check relevant files in setup/debug folder!')
            else:
                setup_logging.showInfo(logger, 'Maximum distance between points of maximum discharge falls within allowed tolerance. Assumption used for calculation of AnnualDischarge is deemed to be valid.')
        else:
            # calculate distance of diagonal line in catchment bounds
            catchment_bounds = catchment.geometry[0].bounds
            dist_diagonal = haversine((catchment_bounds[0],catchment_bounds[1]), (catchment_bounds[2],catchment_bounds[3]))
            dist_frac = max_dist/dist_diagonal
            setup_logging.showInfo(logger, 'Catchment bounds diagonal distance is ' + str(np.round(dist_diagonal,2)) + ' km')
            setup_logging.showInfo(logger, 'Fractional tolerance for maximum distance is ' + str(tolerance_frac))
            setup_logging.showInfo(logger, 'Distance fraction is ' + str(np.round(dist_frac,4)))
            if dist_frac > tolerance_frac:
                setup_logging.showWarning(logger, 'Maximum distance between points of maximum discharge exceeds allowed tolerance! Cannot guarantee that assumption used for calculation of AnnualDischarge is valid! Please check relevant files in setup/debug folder!')
            else:
                setup_logging.showInfo(logger, 'Maximum distance between points of maximum discharge falls within allowed tolerance. Assumption used for calculation of AnnualDischarge is deemed to be valid.')
        # write to geojson file
        if debug_path != None:
            try:
                props = {'max_dist':np.round(max_dist,2), 'dist_frac':np.round(dist_frac,4), 'dist_diagonal':np.round(dist_diagonal,2)}
                writeToGeoJSON(max_coords, distances, props, debug_path)
                setup_logging.showInfo(logger, 'GeoJSON of maximum discharge points created at ' + debug_path)
            except:
                setup_logging.showWarning(logger, 'Could not write points of maximum discharge to GeoJSON file!')
        # create and save figures
        if debug_path != None:
            try:
                writeToFigureTimeSeries(max_discharges, debug_path) 
                writeToFigureMapPoints(ds_catchment.qav.max('time'), max_coords, debug_path)
                setup_logging.showInfo(logger, 'Figures of catchment discharges and maximum discharge points created at ' + debug_path)
            except:
                setup_logging.showWarning(logger, 'Could not create figures of catchment discharges and maximum discharge points!')
    #return float(mean_discharge)
    return int(mean_discharge)

