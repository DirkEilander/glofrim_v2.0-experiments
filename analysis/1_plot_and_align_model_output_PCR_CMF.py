# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 13:19:47 2018

This script loads (pre-processed) CMF and PCR output as well as observed discharge.
It subsequently alignes the lengths of the time series in case they differ.
Time series are plotted for both Hardinge Bridge and Bahadurabad.
For further analysis of for instance the KGE, the aligned time series are saved to file.

@author: J.M. Hoch
@mail: j.m.hoch@uu.nl
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import datetime

# DEFINE PLOTTING PERIOD
# since it does not need to correspond with entire modelling period
def date_range(start, end):
    r = (end+datetime.timedelta(days=1)-start).days
    return [start+datetime.timedelta(days=i) for i in range(r)]
start = datetime.date(2000,01,01)
end = datetime.date(2009,12,31)
end_plot = datetime.date(2010,01,01)
date = date_range(start, end)

# LOAD FILES
# observed discharge
obs_hab_fo = r'path/to/hab_iwm.txt' # at Hardinge Bridge
obs_bah_fo = r'path/to/bah_iwm.txt' # at Bahadurabad
obs_hab = np.loadtxt(obs_hab_fo)
obs_bah = np.loadtxt(obs_bah_fo)

# CMF discharge
# from script '0_preprocess_and_sample_CMF_output.py'
CMF_hab_fo = r'path/to/sampled_CMF_output/QCMF_at_HardingeBridge.txt'
CMF_bah_fo = r'path/to/sampled_CMF_output/QCMF_at_Bahadurabad.txt'
CMF_hab = np.loadtxt(CMF_hab_fo)
CMF_bah = np.loadtxt(CMF_bah_fo)

# PCR dischage
# using dump-files here obtained from ncview
PCR_hab_fo = r'path/to/sampled_PCR_output/QPCR_at_HardingeBridge.dump'
PCR_bah_fo = r'path/to/sampled_PCR_output/QPCR_at_Bahadurabad.dump'
PCR_hab = np.loadtxt(PCR_hab_fo, skiprows=4)[:,1]
PCR_bah = np.loadtxt(PCR_bah_fo, skiprows=4)[:,1]

# FIND COMMON TIME PERIOD
# in case time lengths differ
minLength = min(len(PCR_bah), len(CMF_bah))
lag1 = len(CMF_bah) - minLength

# LIMIT ARRAYS TO COMMON TIME PERIOD
obs_hab_mL = obs_hab[:minLength]
obs_bah_mL = obs_bah[:minLength]
CMF_hab_mL = CMF_hab[lag1:minLength+lag1]
CMF_bah_mL = CMF_bah[lag1:minLength+lag1]
PCR_hab_mL = PCR_hab[:minLength]
PCR_bah_mL = PCR_bah[:minLength]
date = date[:minLength]

# DEFINE TIME PERIOD FOR PLOT AND TIMESERIES OUTPUT
d1 = datetime.date(2004,01,01)
d2 = datetime.date(2009,01,01)
delta_d1 = (d1 - start).days
delta_d2 = (d2 - start).days

# PLOT OBSERVED AND SIMULATED DISCHAGE
# at Harbinge Bridge
fig = plt.figure(1, figsize=(40,20))
ax1 = fig.add_subplot(211)
l_obs, = ax1.plot_date(x=date, y=obs_hab_mL/1000, fmt=':', c='k', linewidth=8)
l_pcr, = ax1.plot_date(x=date, y=PCR_hab_mL/1000, fmt='-', c='r', linewidth=8)
l_cmf, = ax1.plot_date(x=date, y=CMF_hab_mL/1000, fmt='-', c='b', linewidth=8)
ax1.legend((l_obs, l_pcr, l_cmf, ), ('OBS', 'PCR-DynRout', 'PCR->CMF', ), loc=1, ncol=3, frameon=True)
ax1.set_yticks(np.arange(20,60.1,20))
ax1.set_ylim(0,60)
ax1.set_xlim(d1, d2)
ax1.xaxis.labelpad = 1
plt.setp(ax1.get_yticklabels())
plt.setp(ax1.get_xticklabels(), visible=False)

# at Bahadurabad
ax2 = fig.add_subplot(212)
ax2.plot_date(x=date, y=obs_bah_mL/1000, fmt=':', c='k', linewidth=8)
ax2.plot_date(x=date, y=PCR_bah_mL/1000, fmt='-', c='r', linewidth=8)
ax2.plot_date(x=date, y=CMF_bah_mL/1000, fmt='-', c='b', linewidth=8)
ax2.set_yticks(np.arange(30,90.1,30))
ax2.set_ylim(0,90)
ax2.set_xlim(d1, d2)
ax2.xaxis.labelpad = 1
plt.setp(ax2.get_yticklabels())
plt.setp(ax2.get_xticklabels(), visible=True)
plt.text(0.05,0.9,'b)',fontweight='bold')

pos1 = ax1.get_position()
posX0leg = pos1.x0 - 0.125
fig.text(posX0leg, 0.5, 'discharge Q [$10^3 m^3 s^{-1}$]', va='center', rotation='vertical', fontweight=22)

plt.tight_layout()

# SAVE ALIGNED TXT FILES FOR ANALYSIS IN R-STUDIO
# observed discharge
np.savetxt(r'path/to/validationFolder/hab_iwm_mL.txt', obs_hab_mL[delta_d1:delta_d2])
np.savetxt(r'path/to/validationFolder/bah_iwm_mL.txt', obs_bah_mL[delta_d1:delta_d2])
# CMF discharge
np.savetxt(r'path/to/validationFolder/CMF_HardingeBridge_aligned.txt', CMF_hab_mL[delta_d1:delta_d2])
np.savetxt(r'path/to/validationFolder/CMF_Bahadurabad_aligned.txt', CMF_bah_mL[delta_d1:delta_d2])
# PCR discharge
np.savetxt(r'path/to/validationFolder/PCR_HardingeBridge_aligned.txt', PCR_hab_mL[delta_d1:delta_d2])
np.savetxt(r'path/to/validationFolder/PCR_Bahadurabad_aligned.txt', PCR_bah_mL[delta_d1:delta_d2])


