# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 13:19:47 2018

This script plots the discharge at a common location in the Lisflood-FP (LFP)
model domain.
It requires the time series for all models: CMF, LFP, and PCR, and is used
to inspect the accuracy of each model output qualitatively.

@author: J.M. Hoch
@mail: j.m.hoch@uu.nl
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import datetime
import os

# DEFINE PLOTTING PERIOD
# since it does not need to correspond with entire modelling period
def date_range(start, end):
    r = (end+datetime.timedelta(days=1)-start).days
    return [start+datetime.timedelta(days=i) for i in range(r)]
start = datetime.date(2000,01,01)
end = datetime.date(2009,12,31)
end_plot = datetime.date(2010,01,01)
date = date_range(start, end)

# CMF discharge
# from script '0_preprocess_and_sample_CMF_output.py'
CMF_bench2LFP_fo = r'path/to/sampled_CMF_output/QCMF_at_LFPlocation.txt'
CMF_bench2LFP = np.loadtxt(CMF_bench2LFP_fo)

# PCR dischage
# using dump-files here obtained from ncview
PCR_bench2LFP_fo = r'path/to/sampled_PCR_output/QPCR_at_LFPlocation.dump'
PCR_bench2LFP = np.loadtxt(PCR_bench2LFP_fo, skiprows=4)[:,1]

# function reads the LFP discharge from file and plots it for a given observation station
def printLFPresults(fo , rows, leg, stationNr):
    do = np.loadtxt(fo,skiprows=rows)
    ts = len(do)
    ts = np.arange(0,ts,1)
        
    q = do[:,stationNr]

    return q

# sample LFP discharge from file for observation station
LFP_fo = r'path/to/sampled_LFP_output/QLFP.discharge'
QLFP = printLFPresults(LFP_fo, rows=3, leg='LFP', stationNr=-1)

# FIND COMMON TIME PERIOD
# in case time lengths differ
minLength = min(len(PCR_bench2LFP), len(CMF_bench2LFP))
lag1 = len(CMF_bench2LFP) - minLength

# LIMIT ARRAYS TO COMMON TIME PERIOD
CMF_bench2LFP_mL = CMF_bench2LFP[lag1:minLength+lag1]
PCR_bench2LFP_mL = PCR_bench2LFP[:minLength]
QLFP_mL = QLFP[:minLength]
date = date[:minLength]

# DEFINE TIME PERIOD FRO PLOT AND TIMESERIES OUTPUT
d1 = datetime.date(2004,01,01)
d2 = datetime.date(2009,01,01)
delta_d1 = (d1 - start).days
delta_d2 = (d2 - start).days

# PLOT OBSERVED AND SIMULATED DISCHAGE
fig = plt.figure(1, figsize=(40,20))
ax1 = fig.add_subplot(111)
l_pcr, = ax1.plot_date(x=date, y=PCR_bench2LFP_mL/1000, fmt='-', c='r', linewidth=8)
l_cmf, = ax1.plot_date(x=date, y=CMF_bench2LFP_mL/1000, fmt='-', c='b', linewidth=8)
l_lfp, = ax1.plot_date(x=date, y=QLFP_mL/1000, fmt='-', c='g', linewidth=8)
ax1.legend((l_pcr, l_cmf, l_lfp), ('PCR-DynRout', 'PCR->CMF', 'PCR->CMF->LFP'), loc=2, ncol=3, frameon=True)
ax1.set_yticks(np.arange(10,100.1,20))
ax1.set_ylim(0,100)
ax1.set_xlim(d1, d2)
ax1.xaxis.labelpad = 1
plt.setp(ax1.get_yticklabels())
plt.setp(ax1.get_xticklabels(), visible=True)
ax1.set_ylabel('discharge Q [$10^3 m^3 s^{-1}$]')

plt.tight_layout()