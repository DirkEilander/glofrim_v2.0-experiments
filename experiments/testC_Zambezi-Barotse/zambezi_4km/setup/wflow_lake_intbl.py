# -*- coding: utf-8 -*-
"""
Created on Thu Aug 8 11:06:00 2019

@author: haag

Writes .tbl files required for lake module in wflow model. Also allows for the creation of figures and logging for a quick overview of output data and easier debugging.

Partially based on 'v2_wflow_reservoir_intbl.py' from Verseveld.
"""

import geopandas as gpd
import numpy as np
import os

import setup_logging

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def writeSingleArrayToFigure(array, path, title, filename, ylabel, extension='png'):
	"""
	Writes a single array to figure.
	
	Parameters:
		array:     [array] 1D array of factors
		path:      [string] path to directory to place figure
		title:     [string] title of figure
		filename:  [string] name of to-be-constructed figure (without extension)
		ylabel:    [string] text for label on y-axis
		extension: [string] extension/type of figure
	
	output: png file
	"""
	# plot data
	plt.plot(array)
	plt.gca().xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
	if title != None:
		plt.title(title)
	if ylabel != None:
		plt.ylabel(ylabel)
	# create full path
	full_path = os.path.join(path, filename+'.'+extension)
	# check if file already exists, and if so, remove
	if os.path.exists(full_path):
		os.remove(full_path)
	# write figure to file
	plt.savefig(full_path, bbox_inches='tight')
	plt.close()


def make_tbls(lakes, outdir, debug_path=None, logger=None):
	"""
	Creates tbl files related to lakes for use in a wflow model.
	
	Parameters:
		lakes : [geopandas object] lakes within catchment
		outdir : [string] path where to create the new lake tbl files
		debug_path: [string] path to directory where setup/logging/debug files are stored
		logger : [logging object] instance of Python logging module
	
	output: lake tbl files
	"""
	
	# tbl files for wflow lake function
	f_Lake_b       = open(os.path.join(outdir, 'Lake_b.tbl'), 'w')
	f_LakeAvgLevel = open(os.path.join(outdir, 'LakeAvgLevel.tbl'), 'w')
	f_LakeAvgOut   = open(os.path.join(outdir, 'LakeAvgOut.tbl'), 'w')
	f_LakeArea     = open(os.path.join(outdir, 'LakeArea.tbl'), 'w')
	
	# initialization of variables for storing setup/debug information (used to create relevant files at the end)
	if (debug_path != None) and (len(lakes) > 1):
		lake_bs = np.zeros(len(lakes))
		avglvls = np.zeros(len(lakes))
		avgouts = np.zeros(len(lakes))
		areas   = np.zeros(len(lakes))
	
	setup_logging.showInfo(logger, 'Number of lakes in catchment: ' + str(len(lakes)) + '. Writing intbls per lake...')
	
	for i, v in lakes.iterrows():
		
		#get index/count of loop (only used for logger/debug)
		if logger != None or debug_path != None:
			loop_index = np.where(lakes.index==i)[0][0]
		
		setup_logging.showInfo(logger, 'i: ' + str(loop_index) + ', lake index: ' + str(i) + ', HydroLAKE ID: ' + str(v.Hylak_id))
		
		#calculate relevant variables
		lake_b    = v.Dis_avg / v.Depth_avg**2 # alpha parameter from rating curve, simplified version for global setup [-]
		avglvl    = v.Depth_avg # [m]
		avgout    = v.Dis_avg # [m3/s]
		lake_area = v.Lake_area * 1e6 # [m2]
		
		#write wflow lake tbl files
		f_Lake_b.write(str(v.Hylak_id) + '\t' + str(lake_b) + '\n')
		f_LakeAvgLevel.write(str(v.Hylak_id) + '\t' + str(avglvl) + '\n')
		f_LakeAvgOut.write(str(v.Hylak_id) + '\t' + str(avgout) + '\n')
		f_LakeArea.write(str(v.Hylak_id) + '\t' + str(lake_area) + '\n')
		
		#store relevant values for creation of setup/debug files
		if (debug_path != None) and (len(lakes) > 1):
			lake_bs[loop_index] = lake_b
			avglvls[loop_index] = avglvl
			avgouts[loop_index] = avgout
			areas[loop_index]   = lake_area
	
	f_Lake_b.close()
	f_LakeAvgLevel.close()
	f_LakeAvgOut.close()
	f_LakeArea.close()
	
	if (debug_path != None) and (len(lakes) > 1):
		try:
			writeSingleArrayToFigure(lake_bs, debug_path, 'Lake Lake_b', 'Lakes_Lake_b', '')
			writeSingleArrayToFigure(avglvls, debug_path, 'Lake average water level', 'Lakes_AvgLevel', 'Water level (m)')
			writeSingleArrayToFigure(avgouts, debug_path, 'Lake average outflow', 'Lakes_AvgOut', 'Outflow (m3/s)')
			writeSingleArrayToFigure(areas/1e6, debug_path, 'Lake average area', 'Lakes_Area', 'Area (km2)')
			setup_logging.showInfo(logger, 'Figures related to lake intbls created at ' + debug_path)
		except:
			setup_logging.showWarning(logger, 'Could not create figures related to lake intbls!')
