# -*- coding: utf-8 -*-
"""
Created on Thu Aug 8 10:12:00 2019

@author: haag

Creates waterbody (reservoirs and/or natural lakes) maps for use in wflow models. Part of the wflow model setup scripts.

Adapted from script 'reservoirposition.py' from Marnix van der Vat, 180606, version 0.0.
"""

import geopandas as gpd
import fiona as fiona
import rasterio as rasterio
from rasterio import features
import gdal
import numpy as np
import math as math
import os
import pcraster as pcr

import setup_logging

def make_maps(waterbodies, do_lakes, outdir, area_filename, outflow_filename, limsnp=50, mxsnp=0, logger=None):
	"""
	Creates map files related to waterbodies (either reservoirs or natural lakes) for use in a wflow model.
	
	Parameters:
		waterbodies : [geopandas object] waterbodies within catchment
		do_lakes : [boolean] run for lakes (otherwise runs for reservoirs)
		outdir : [string] path where to create the new waterbody map files
		limsnp : [int] limit in flow accumulation value below which no snapping takes place
		mxsnp: [int] max distance in cell numbers for snapping of waterbody to river in flow accumulation (cells)
		logger : [logging object] instance of Python logging module
	
	output: waterbody maps in specified output directory
	"""
	
	#helper function to create PCRaster .map file by first creating a temporary GeoTIFF file and then converting that to desired output
	def createMapFile(data, outdir, mapfile):
		#create temporary GeoTIFF file
		setup_logging.showInfo(logger, 'Staring process to write ' + os.path.join(outdir, mapfile + '.map'))
		setup_logging.showInfo(logger, 'Writing temporary GeoTIFF file...')
		with rasterio.open(os.path.join(outdir, mapfile + '.tif'), 'w', **map_meta) as out:
			out.write(data, 1)
		out.close()
		#convert temporary GeoTIFF file to PCRaster .map file
		setup_logging.showInfo(logger, 'Converting temporary GeoTIFF file to PCRaster map file...')
		gdal.Translate(os.path.join(outdir, mapfile + '.map'), os.path.join(outdir, mapfile + '.tif'), options = '-of PCRaster')
		#delete temporary GeoTIFF file (and .map.aux.xml file created during conversion)
		setup_logging.showInfo(logger, 'Deleting temporary GeoTIFF file and .map.aux.xml file that is created during conversion..')
		if os.path.exists(os.path.join(outdir, mapfile + '.tif')):
			os.remove(os.path.join(outdir, mapfile + '.tif'))
		else:
			setup_logging.showWarning(logger, 'Could not find temporary GeoTIFF file! Skipping removal.')
		if os.path.exists(os.path.join(outdir, mapfile + '.map.aux.xml')):
			os.remove(os.path.join(outdir, mapfile + '.map.aux.xml'))
		else:
			setup_logging.showWarning(logger, 'Could not find .map.aux.xml file created during conversion! Skipping removal.')
	
	#helper function to convert PCRaster .map file to nominal data type
	def convertToNominal(mapdir, mapfile):
		mapfile += '.map'
		setup_logging.showInfo(logger, 'Converting to nominal type: ' + mapfile)
		try:
			temp_file = pcr.readmap(os.path.join(mapdir, mapfile))
			pcr.report(pcr.nominal(temp_file), os.path.join(mapdir, mapfile))
		except:
			setup_logging.showWarning(logger, 'Could not convert ' + mapfile + '!')
	
	#1. Define and read dummy input that will come from within the MapBuilder and from the ini file
	if not os.path.exists(os.path.join(outdir, 'wflow_acc.map')):
		#calculate flow accumulation (PCRaster Python)
		pcr._pcraster.setclone(os.path.join(outdir, 'wflow_dem.map'))
		setup_logging.showInfo(logger, 'Calculating flow accumulation...')
		ldd = pcr.readmap(os.path.join(outdir, 'wflow_ldd.map'))
		flow_acc = pcr.accuflux(ldd, 1)
		setup_logging.showInfo(logger, 'Writing flow accumulation to wflow_acc.map')
		pcr.report(flow_acc, os.path.join(outdir, 'wflow_acc.map'))
	
	#read the flow accumulation and obtain map metadata
	accummap=rasterio.open(os.path.join(outdir, 'wflow_acc.map'),dtype=np.uint)
	accum=accummap.read(1)
	map_meta = accummap.meta.copy()
	if map_meta['crs'] == None:
		setup_logging.showWarning(logger, 'No CRS found in map metadata, CRS set to EPSG:4326!')
		map_meta['crs'] = rasterio.crs.CRS.from_epsg(4326)
	map_meta.update(compress='lzw')
	map_meta['driver']='GTiff' #make sure driver is set to GeoTIFF (instead of PCRaster)
	
	#read the catchment
	with rasterio.open(os.path.join(outdir, "wflow_subcatch.map")) as src: 
		catchment=src.read(1)
	src.close()
	
	#define class for waterbody data
	class cls_Rsv:
		def __init__(self,id,name):
			self.id=id
			self.name=name
			self.maxaccum=-1
			self.Xmaxaccum=-999
			self.Ymaxaccum=-999
			self.Rmaxaccum=-999
			self.Cmaxaccum=-999
	
	if do_lakes:
		waterbody_ids   = waterbodies.Hylak_id
		waterbody_names = waterbodies.Lake_name
	else:
		waterbody_ids   = waterbodies.ID
		waterbody_names = waterbodies.GRanD_DAM_
	
	rsvs=list()
	IDs=list()
	for id,name in zip(waterbody_ids,waterbody_names): 
		rsvs.append(cls_Rsv(id,name))
		IDs.append(id)
	
	#create the rasterized map of waterbody IDs
	out_arr = pcr.pcr2numpy(pcr.readmap(os.path.join(outdir, 'wflow_acc.map')), map_meta['nodata'])
	out_arr = (out_arr/out_arr)-1 #make sure default array contains zeros only
	shapes = ((geom,value) for geom, value in zip(waterbodies.geometry, waterbody_ids))
	waterbody = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=map_meta['transform'], all_touched=True)
	waterbody[waterbody==0] = map_meta['nodata'] # prevent zeros in output map
	createMapFile(waterbody, outdir, area_filename)
	convertToNominal(outdir, area_filename)
	
	#2. determine location of maximum flow accumulation in each waterbody
	nrow=waterbody.shape[0]
	ncol=waterbody.shape[1]
	for irow in range(nrow):
		for icol in range(ncol):
			if catchment[irow,icol] > 0:
				id = waterbody[irow,icol]
				if id > 0:
					i=IDs.index(id)
					ac = accum[irow,icol]
					if ac > rsvs[i].maxaccum:
						rsvs[i].maxaccum=ac
						rsvs[i].Rmaxaccum=irow
						rsvs[i].Cmaxaccum=icol
	
	#3. snap waterbody location to a higher accum value with the highest gradient, if present in selected window
	resmap = np.zeros((nrow,ncol),dtype=np.float32)
	for rsv in rsvs:
		#snap only if maxaccum below the limit and above 0
		ac = rsv.maxaccum
		if ac > 0:
			py=rsv.Rmaxaccum
			px=rsv.Cmaxaccum
			mx=px
			my=py
			if ac < limsnp:
				gradmx = 0.
				for ny in range(py-mxsnp, py+mxsnp,1):
					for nx in range(px-mxsnp,px+mxsnp,1):
						dist= math.sqrt((px-nx)**2+(py-ny)**2)
						if dist>0: grad=(accum[ny,nx]-ac)/dist
						else: grad=-1.
						if grad > gradmx:
							mx=nx
							my=ny
							gradmx=grad
			resmap[my,mx] = float(rsv.id)
	resmap[resmap==0] = map_meta['nodata'] # prevent zeros in output map
	
	#4. Write map file with the outflow point of each waterbody indentified by its ID
	createMapFile(resmap, outdir, outflow_filename)
	convertToNominal(outdir, outflow_filename)
	
	#5. Check if rasterized waterbodies were created properly (continuous area, no 'splits')
	unique_input = len(np.unique(waterbody)) - 1
	continuous_waterbodies = pcr.clump(pcr.numpy2pcr(pcr.Nominal, waterbody, map_meta['nodata']))
	unique_output = len(np.unique(pcr.pcr2numpy(continuous_waterbodies, map_meta['nodata']))) - 1
	if unique_input != unique_output:
		setup_logging.showWarning(logger, 'Number of unique output waterbodies does not match input! This can result in duplicate IDs for different areas! Check ' + area_filename + '.map')
