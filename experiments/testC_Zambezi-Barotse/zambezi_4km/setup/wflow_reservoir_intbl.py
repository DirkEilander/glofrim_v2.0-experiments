# -*- coding: utf-8 -*-
"""
Created by Willem van Verseveld. Adapted by Arjen Haag, March/April 2019.

Writes .tbl files required for reservoir module in wflow model. Also allows for the creation of figures and logging for a quick overview of output data and easier debugging.
"""

import geopandas as gpd
import hydroengine as he
import numpy as np
import reservoir_intbl_utils as ut
import os

import setup_logging

import matplotlib.pyplot as plt


def writeMultiArrayToFigure(arrays, labels, path, filename, extension='png', title=None):
	"""
	Writes multiple arrays to figure.
	
	Parameters:
		arrays:    [list] list of 1D arrays
		labels:    [list] list of labels for legend
		path:      [string] path to directory to place figure
		filename:  [string] name of to-be-constructed figure (without extension)
		extension: [string] extension/type of figure
	
	output: png file
	"""
	# plot data
	for i in range(len(arrays)):
		plt.plot(arrays[i], label=labels[i])
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	if title == None:
		plt.title(filename.replace('_', ' '))
	# create full path
	full_path = os.path.join(path, filename+'.'+extension)
	# check if file already exists, and if so, remove
	if os.path.exists(full_path):
		os.remove(full_path)
	# write figure to file
	plt.savefig(full_path, bbox_inches='tight')
	plt.close()


def writeSingleArrayToFigure(array, path, title='Reservoir increase factors', filename='Reservoir_increase_factors', extension='png'):
	"""
	Writes a single array to figure.
	
	Parameters:
		array:     [array] 1D array of factors
		title:     [string] title of figure
		path:      [string] path to directory to place figure
		filename:  [string] name of to-be-constructed figure (without extension)
		extension: [string] extension/type of figure
	
	output: png file
	"""
	# plot data
	plt.plot(array)
	if title != None:
		plt.title(title)
	# create full path
	full_path = os.path.join(path, filename+'.'+extension)
	# check if file already exists, and if so, remove
	if os.path.exists(full_path):
		os.remove(full_path)
	# write figure to file
	plt.savefig(full_path, bbox_inches='tight')
	plt.close()


def writeHistogram(data, labels, path, filename='Reservoir_accuracy_levels_histogram', extension='png'):
	"""
	Writes histogram figure.
	
	Parameters:
		data:      [list|array] (list of) 1D array(s)
		labels:    [list|string] (list of) label(s) for legend
		path:      [string] path to directory to place figure
		filename:  [string] name of to-be-constructed figure (without extension)
		extension: [string] extension/type of figure
	
	output: png file
	"""
	plt.hist(data, np.arange(-4,5)-0.5, label=labels)
	plt.title('Histogram of reservoir accuracy levels')
	plt.legend(loc='upper right')
	# create full path
	full_path = os.path.join(path, filename+'.'+extension)
	# check if file already exists, and if so, remove
	if os.path.exists(full_path):
		os.remove(full_path)
	# write figure to file
	plt.savefig(full_path, bbox_inches='tight')
	plt.close()


def make_tbls(reservoirs, outdir, perc_norm=50, perc_min=20, mv=-99, range_min=[0.0,1.0], range_max=[0.0,1.0], method=1, debug_path=None, logger=None):
	"""
	Creates tbl files related to reservoirs for use in a wflow model.
	
	Parameters:
		reservoirs : [geopandas object] reservoirs within catchment
		outdir : [string] path where to create the new reservoir tbl files
		perc_norm : [int] percentile for normal (operational) surface area
		perc_min: [int] percentile for minimal (operational) surface area
		mv : [int] missing value
		method [int] controls which method to use for calculations (0 = original / with factor, 1 = without factor)
		debug_path: [string] path to directory where setup/logging/debug files are stored
		logger : [logging object] instance of Python logging module
	
	output: reservoir tbl files
	"""
	
	# check method
	if method not in [0,1]:
		setup_logging.showWarning(logger, 'Unrecognized input for method of reservoir intbl calculation! Using default value instead.')
		method = 1
	
	# tbl files for wflow simple reservoir function
	f_ResDemand = open(os.path.join(outdir, 'ResDemand.tbl'), 'w')
	f_ResMaxRelease = open(os.path.join(outdir, 'ResMaxRelease.tbl'), 'w')
	f_ResMaxVolume = open(os.path.join(outdir, 'ResMaxVolume.tbl'), 'w')
	f_ResTargetFullFrac = open(os.path.join(outdir, 'ResTargetFullFrac.tbl'), 'w')
	f_ResTargetMinFrac = open(os.path.join(outdir, 'ResTargetMinFrac.tbl'), 'w')
	f_ResSimpleArea = open(os.path.join(outdir, 'ResSimpleArea.tbl'), 'w')
	
	# initialization of variables for storing setup/debug information (used to create relevant files at the end)
	if (debug_path != None) and (len(reservoirs) > 1):
		min_fracs = np.zeros(len(reservoirs))
		max_fracs = np.zeros(len(reservoirs))
		factors   = np.zeros(len(reservoirs))
		acc_mins  = np.zeros(len(reservoirs))
		acc_norms = np.zeros(len(reservoirs))
		acc_dict  = {'1':0,'2':1,'3':2,'21':-1,'31':-2,'4':-3,'5':-4} # conversion of accuracy levels for better plotting
	
	setup_logging.showInfo(logger, 'Number of reservoirs in catchment: ' + str(len(reservoirs)) + '. Writing intbls per reservoir...')
	
	for i, v in reservoirs.iterrows():
		
		#get index/count of loop (only used for logger/debug)
		if logger != None or debug_path != None:
			loop_index = np.where(reservoirs.index==i)[0][0]
		
		setup_logging.showInfo(logger, 'i: ' + str(loop_index) + ', reservoir index: ' + str(i) + ', HydroLAKE ID: ' + str(v.Hylak_id))
		
		try:
			time_series = he.get_lake_time_series(v.Hylak_id,'water_area')
			area_series = np.array(time_series['water_area']) #[m2]
			area_series_nozeros =  area_series[area_series > 0]
			max_areaJRC = area_series_nozeros.max()/10000 #[ha]
			norm_area = np.percentile(area_series_nozeros, perc_norm, axis=0)/10000 #[ha]
			min_area = np.percentile(area_series_nozeros, perc_min, axis=0)/10000 #[ha]
		except Exception as e:
			if logger != None:
				logger.error(str(e))
				logger.warning('No HydroEngine time series availalable!')
			else:
				print(str(e))
				print('No HydroEngine time series availalable!')
			#TO DO: add lakes from HydroLAKES (not all lakes are part of GRanD database)
			max_areaJRC = mv
			norm_area = mv
			min_area = mv
		
		#calculate relevant variables
		Lake_area = v.Lake_area * 100
		max_area, min_cap, norm_cap, max_cap, factor, acc_min, acc_norm = \
			ut.reservoir_parameters(v.G_CAP_MIN,v.G_CAP_REP,v.Vol_total,min_area,norm_area,max_areaJRC,Lake_area,v.G_DAM_HGT_,v.Depth_avg,mv,method,logger)
		
		#estimate demand reservoir downstream (ResDemand) and maximum release (ResMaxRelease)
		#based on average annual discharge
		if v.Dis_avg == mv:
			ResDemand = mv
			ResMaxRelease = mv
		else:
			ResDemand = v.Dis_avg/2.0
			ResMaxRelease = v.Dis_avg*4.0
		
		#calculate fractions
		min_frac  = min_cap/max_cap
		full_frac = norm_cap/max_cap
		
		# limit fractions to specified range
		min_frac = min(max(min_frac, range_min[0]), range_min[1])
		full_frac = min(max(full_frac, range_max[0]), range_max[1])
		
		#write wflow reservoir tbl files
		f_ResDemand.write(str(v.ID) + '\t' + str(ResDemand) + '\n')
		f_ResMaxRelease.write(str(v.ID) + '\t' + str(ResMaxRelease) + '\n')
		f_ResMaxVolume.write(str(v.ID) + '\t' + str(max_cap * 1000000) + '\n') #[m3]
		f_ResTargetFullFrac.write(str(v.ID) + '\t' + str(full_frac) + '\n')
		f_ResTargetMinFrac.write(str(v.ID) + '\t' + str(min_frac) + '\n')
		f_ResSimpleArea.write(str(v.ID) + '\t' + str(max_area*10000) + '\n') #[m2]
		
		#store relevant values for creation of setup/debug files
		if (debug_path != None) and (len(reservoirs) > 1):
			min_fracs[loop_index] = min_frac
			max_fracs[loop_index] = full_frac
			factors[loop_index]   = factor
			acc_mins[loop_index]  = acc_dict[str(acc_min)]
			acc_norms[loop_index] = acc_dict[str(acc_norm)]
	
	f_ResDemand.close()
	f_ResDemand.close()
	f_ResMaxRelease.close()
	f_ResMaxVolume.close()
	f_ResTargetFullFrac.close()
	f_ResTargetMinFrac.close()
	f_ResSimpleArea.close()
	
	if (debug_path != None) and (len(reservoirs) > 1):
		try:
			writeMultiArrayToFigure([min_fracs, max_fracs], ['MinFrac', 'FullFrac'], debug_path, 'Reservoir_capacity_fractions')
			writeMultiArrayToFigure([acc_mins, acc_norms], ['min', 'norm'], debug_path, 'Reservoir_accuracy_levels')
			writeSingleArrayToFigure(factors, debug_path)
			writeHistogram([acc_mins, acc_norms], ['min', 'norm'], debug_path)
			setup_logging.showInfo(logger, 'Figures related to reservoir intbls created at ' + debug_path)
		except:
			setup_logging.showWarning(logger, 'Could not create figures related to reservoir intbls!')
