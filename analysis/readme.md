# The analysis scripts

Here the python and R scripts applied for the analysis of the model results are provided.
You can find four files whose content and objective are briefly outlined hereafter. For more detailed information, please consult the scripts themselves.

1. **0_preprocess_and_sample_CMF_output.py**: CaMa-Flood (CMF) can be compiled with or without netCDF-support. In case it is compiled without netCDF, output files will be written as bindary-files (.bin) per model year. To be able to retrieve output information at a specific location from those files, it is necessary to reconstruct the geometry, to locate the corresponding grid indices, and to sample the model results at those indices for all model steps. In case CMF is compiled with netCDF, this file may not be needed.

1. **1_plot_and_align_model_output_PCR_CMF.py**: Model output from PCR-GLOBWB (PCR) and CMF may be present for different simulation periods. To clip the output timeseries as well as the time series of observed values to the common period, this script can be used. Also, it plots the results from both models for a visual inspection of the results. Time series of the common period are saved to file to then be used for a analysis of the Kling-Gupta Efficiency (KGE) in the file '3_hydroGOF_KGEanalysis.R'.

1. **2_plot_benchmark_PCR_CMF_LFP.py**: For a visual benchmarking of the differences between output obtained from PCR, CMF, and Lisflood-FP (LFP), this script plots output from all three models for a consistent time period.

1. **3_hydroGOF_KGEanalysis.R**: This R-script requires time series of observed and simulated discharge, such as obtained with script 1. It then can calculate the KGE including its individual components (linear correlation r, bias ratio β, and variability α) to provide a qualitative assessment of model output accuracy.

The developers tried to document the scripts as clear and extensive as possible. In case problems are encountered or aspects of the scripts are unclear, please do not hesitate to post an issue in the corresponding GitHub-repository.

Please also note that the scripts were specifically designed to analyse the output of the test cases presented in the accompanying NHESS-article. In case you want to apply them for other set-ups than that, please be aware that they may not work without making necessary adjustments. 