# Script for time series analysis of simulated discharge at Obidos
# input time series have prevoiusly been prepared in Python
# to make this script work, first set working directory to file directory
# Session -> Set Working Directory -> To Source File Location

# Install stable version 
#install.packages("hydroGOF")

setwd('path/to/workingDirectory')
# for instance 'validationFolder' based on py-scripts

## HARDINGE BRIDGE
PCR_hab = read.csv('PCR_HardingeBridge_aligned.txt')
CMF_hab = read.csv('CMF_HardingeBridge_aligned.txt')
OBS_hab = read.csv('hab_iwm_mL.txt')

KGE_PCR_hab=hydroGOF::KGE(PCR_hab, OBS_hab, method="2009", out.type="full")
KGE_CMF_hab=hydroGOF::KGE(CMF_hab, OBS_hab, method="2009", out.type="full")

NSE_PCR_hab=hydroGOF::NSE(PCR_hab,OBS_hab)
NSE_CMF_hab=hydroGOF::NSE(CMF_hab,OBS_hab)

## BAHADURABAD
PCR_bah = read.csv('PCR_Bahadurabad_aligned.txt')
CMF_bah = read.csv('CMF_Bahadurabad_aligned.txt')
OBS_bah = read.csv('bah_iwm_mL.txt')

KGE_PCR_bah=hydroGOF::KGE(PCR_bah, OBS_bah, method="2009", out.type="full")
KGE_CMF_bah=hydroGOF::KGE(CMF_bah, OBS_bah, method="2009", out.type="full")

NSE_PCR_bah=hydroGOF::NSE(PCR_bah,OBS_bah)
NSE_CMF_bah=hydroGOF::NSE(CMF_bah,OBS_bah)

