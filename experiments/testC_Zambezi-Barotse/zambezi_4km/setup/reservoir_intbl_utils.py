# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 10:37:29 2018

@author: verseve

Adapted by Arjen Haag, March 2019, to work as part of wflow setup script.

Most notable adaptations include:
- removed all unused output data (and their calculations)
- added an additional option for highest accuracy (in which the area is allowed to be missing, as this is not relevant for wflow)
- checks on output min_cap, to make sure it never exceeds norm_cap (was already in place for norm_cap vs. max_cap)
- an option to use two different methods
    - original method (with the only change being the checks mentioned above)
    - other method, where the potential increase in maximum level (when GranD dam height exceeds HydroLAKES depth average) is only applied on this maximum level
                    (and not forced upon other output variables through the 'factor_shape' and 'lin_coeff' variables).
The reasoning behind the other method is twofold:
    1. When checking various Wflow catchments is was found that the resulting increases very often led to normal and/or minimum values exceeding maximum values,
       which caused them to be capped by checks on the output data (implying that their calculation was incorrect). Without these increases this was no longer the case.
    2. The increases were thought to be applied inconsistently. For example, 'factor_shape' was only applied on the normal level and only when both normal area and volume were not missing.
       The increase was a result of a discrepancy between the datasets of GranD and HydroLAKES but could also be applied to data from JRC, which might not have the same discrepancy, e.g.:
       norm_level = dam_height[GranD]/max_level[HydroLAKES]*(norm_cap[GRanD]/norm_area[JRC])
"""

import setup_logging


def reservoir_parameters(min_cap,norm_cap,max_cap,min_area,norm_area,max_areaJRC,max_areaHL,dam_height,max_level,mv,method=1,logger=None):
    """
    :function is based on ComputationDatabase.py script of Alessia Matano
    
    :param min_cap: minimum reservoir capacity (CAP_MIN from GRanD) Mcm
    :param norm_cap: normal reseroir capacity (CAP_REP from GRanD) Mcm
    :param max_cap: maximum reseroir capacity (Vol_total from HydroLAKES) Mcm
    :param min_area: minumum reservoir area (from JRC) ha
    :param norm_area: normal reservoir area (from JRC) ha
    :param max_areaJRC: maximum reservoir area (from JRC) ha 
    :param max_areaHL: maximum reservoir area (Lake_area from HydroLAKES) ha 
    :param dam_height: dam height (Dam_hgt_m from GRanD) m
    :param max_level: Depth_avg from HydroLAKES m
    :param mv: missing value (e.g. -99) - [should be <= 0 !]
    :param method: 0 = original (with factor), 1 = without factor
    :param logger: instance of Python logging module
    
    :return: 
    max_area, min_cap, norm_cap, max_cap, factor, accuracy_min, accuracy_norm
    """
    
    # Replace all the missing values (-99) with 0
    min_cap     = max(min_cap, 0)
    norm_cap    = max(norm_cap, 0)
    max_cap     = max(max_cap, 0)
    min_area    = max(min_area, 0)
    norm_area   = max(norm_area, 0)
    max_areaJRC = max(max_areaJRC, 0)
    max_areaHL  = max(max_areaHL, 0)
    dam_height  = max(dam_height, 0)
    max_level   = max(max_level, 0)
    
    mv = 0
    
    # Maximum area
    # the source considered more reliable is JRC, so the HydroLAKES value is used only if JRC is missing
    if max_areaJRC != mv:
        max_area = max_areaJRC
    else:
        max_area = max_areaHL
    
    # Maximum level
    # a validation has shown that GRanD dam height is a more reliable value than HydroLAKES depth average (when it is a larger value)
    if dam_height > max_level:
        max_level_f  = dam_height
        factor_shape = dam_height/max_level
        setup_logging.showInfo(logger, 'GRanD dam height used as max. level instead of HydroLAKES depth average. Difference factor: ' + str(round(factor_shape,2)))
    else:
        max_level_f  = max_level
        factor_shape = 1.0
    
    # coefficient for linear relationship
    lin_coeff = max_area/max_level_f #[ha/m]
    
    # adjust factor based on chosen method
    if method == 0:
        factor_used = factor_shape
    elif method == 1:
        factor_used = 1.0
    
    # Operational (norm) level
    if (norm_cap != mv and norm_area != mv):
        norm_level    = factor_used*(norm_cap/norm_area)*100 #[m]
        norm_area_f   = norm_area
        norm_cap_f    = norm_cap
        accuracy_norm = 1
    elif norm_cap != mv:
        norm_level    = ((norm_cap/lin_coeff)**(1/2))*10 #[m]
        norm_area_f   = (norm_cap/norm_level)*100 #[ha]
        norm_cap_f    = norm_cap
        accuracy_norm = 1
    elif norm_area != mv:
        norm_level    = norm_area/lin_coeff #[m]
        norm_area_f   = norm_area
        norm_cap_f    = (norm_area*norm_level)/100 #[Mcm]
        accuracy_norm = 2
    else:
        #the calculation is based on the max area (and not max level or max capacity) as it is the most reliable value
        norm_area_f   = max_area*0.666 #[ha]
        norm_level    = norm_area_f/lin_coeff #[m]
        norm_cap_f    = (norm_area_f*norm_level)/100 #[Mcm]
        accuracy_norm = 3
    
    # CHECK norm level (1)
    if (accuracy_norm == 1 and norm_level > max_level_f):
        #it is assumed that the norm capacity value is not reliable, so this value is delated and the linear relationship assumption is introduced
        norm_level    = norm_area_f/lin_coeff #[m]
        norm_cap_f    = (norm_area_f*norm_level)/100 #[Mcm]
        accuracy_norm = 21
    elif(accuracy_norm == 2 and norm_level > max_level_f):
        norm_area_f   = max_area*0.666 #[ha]
        norm_level    = norm_area_f/lin_coeff #[m]
        norm_cap_f    = (norm_area_f*norm_level)/100 #[Mcm]
        accuracy_norm = 31
    
    # Minimum level
    if (min_area != mv and min_cap != mv):
        min_level    = (min_cap/min_area)*100 #[m]
        min_area_f   = min_area
        min_cap_f    = min_cap
        accuracy_min = 1
    elif min_cap != mv:
        min_level    = ((min_cap/lin_coeff)**(1/2))*10 #[m]
        min_area_f   = (min_cap/min_level)*100 #[ha]
        min_cap_f    = min_cap #[Mcm]
        accuracy_min = 1
    elif min_area != mv:
        min_level    = (min_area/lin_coeff) #[m]
        min_area_f   = min_area #[ha]
        min_cap_f    = (min_area*min_level)/100 #[Mcm]
        accuracy_min = 2
    else:
        #the calculation is based on the max area (and not max level or max capacity) as it is the most reliable value 
        min_area_f   = max_area*0.333 #[ha]  
        min_level    = min_area_f/lin_coeff #[m]
        min_cap_f    = (min_area_f*min_level)/100 #[Mcm]
        accuracy_min = 3
    
    # CHECK minumum level (1)
    if (accuracy_min == 1 and min_level > norm_level):
        accuracy_min = 21
        min_level    = min_area_f/lin_coeff #[ha]
        min_cap_f    = (min_area_f*min_level)/100 #[Mcm]
    elif (accuracy_min == 2 and min_level > norm_level):
        accuracy_min = 31
        min_area_f   = max_area*0.333 #[ha]  
        min_level    = min_area_f/lin_coeff #[m]
        min_cap_f    = (min_area_f*min_level)/100 #[Mcm]
    elif (min_level > norm_level):
        min_area_f   = norm_area_f*0.45 #[ha]  
        min_level    = min_area_f/lin_coeff #[m]
        min_cap_f    = (min_area_f*min_level)/100 #[Mcm]
        accuracy_min = 4
    
    # CHECK norm level (2)
    if norm_cap_f > max_cap:
        setup_logging.showWarning(logger, 'norm_cap > max_cap! setting norm_cap equal to max_cap.' )
        norm_cap_f    = max_cap
        accuracy_norm = 5
    
    # CHECK minumum level (2)
    if min_cap_f > norm_cap_f:
        setup_logging.showWarning(logger, 'min_cap > norm_cap! setting min_cap equal to norm_cap.')
        min_cap_f    = norm_cap_f
        accuracy_min = 5
    
    return max_area, min_cap_f, norm_cap_f, max_cap, factor_shape, accuracy_min, accuracy_norm