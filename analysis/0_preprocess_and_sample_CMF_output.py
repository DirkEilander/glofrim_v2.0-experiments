# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 13:47:00 2018

This is an auxiliary script to pre-process and sample CaMa-FLood (CMF) output from binary files.
Running this script may be required if CMF is not compiled with netCDF. Output from PCR-GLOBWB and 
Lisdflood-FP is either already stored as netCDF files or directly as txt-files and thus does not
require this pre-processing.
After specifying the output location, all binary files are read and simulated discharge is sampled.
For further processing, the values are stored and written to a txt-file.

@author: J.M. Hoch
@mail: j.m.hoch@uu.nl
"""

import rasterio as rs
import matplotlib.pyplot as plt
import numpy as np
import os


# CMF tif file required to obtain spatial properties of output maps

rs_fo = r'path/to/CMF/nextxy.tif'

# using rasterio to derive spatial properties
dataset = rs.open(rs_fo)
a = dataset.transform

# specify the location where output is sampled
# 1: Bahadurabad, 2: Hardinge Bridge, 3: location for comparison with Lisflood-FP
locationID = 3 

if locationID == 1:
    xy_ind = (66,28)
elif locationID == 2:
    xy_ind = (65,32)
elif locationID == 3:
    xy_ind = (70,35)

def get_cmf_outflw_binary(path_base, ncol, nrow, x_ind_obs, y_ind_obs):
    """
    This function reads all binary files in the output folders and appends the simulated discharge
    at the specified location.
    It also plots the study area and the specified location for a visual inspection as well as the
    time series of simulated discharge.
    Finally, it stores the sampled discharge to a txt-file for further processing.
    """
    
    # initiate array for appending discharge
    q = []
    
    # reading binary output files in output filder
    cmf_files = [os.path.join(path_base,f) for f in os.listdir(path_base) if os.path.isfile(os.path.join(path_base,f))]
        
    # for each binary file, the spatial extent has to be reconstructed
    # then, simulated discharge can be sampled for all time steps in the binary file
    for i in xrange(len(cmf_files)):
        fo = cmf_files[i]
        qs = np.fromfile(fo, 'f')
        qs = qs.reshape(len(qs)/(nrow*ncol), nrow, ncol)
        qsim = qs[:, y_ind_obs, x_ind_obs]
        q = np.append(q, qsim) 
    
    # plot discharge map with specified location
    plt.figure(1)
    plt.imshow(np.ma.masked_greater(qs[-1,:,:], 9999999))
    plt.scatter(x_ind_obs, y_ind_obs, c='r')
    
    # plot simulated discharge for all time steps and all binary files
    plt.figure(2)
    plt.plot(q/1000, label='CMF')
    plt.ylim(0,)
    
    return q

# specify folder where binary output is stored and apply function
CMFpath = r'path/to/folder/with/CMF/binOutput'
QCMF = get_cmf_outflw_binary(CMFpath, dataset.width, dataset.height, xy_ind[0], xy_ind[1])

# store output to txt-file per location
if locationID == 1:
    np.savetxt(r'path/to/sampled_CMF_output/QCMF_at_Bahadurabad.txt', QCMF)
elif locationID == 2:
    np.savetxt(r'path/to/sampled_CMF_output/QCMF_at_HardingeBridge.txt', QCMF)
elif locationID == 3:
    np.savetxt(r'path/to/sampled_CMF_output/QCMF_at_LFPlocation.txt', QCMF)

