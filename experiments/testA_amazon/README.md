All model and forcing files for a coupled PCR - CMF run are provided in this folder. 

BEFORE RUNNING
--------------
0. Create a glofrim_exp environment in conda using the environment.yml file. Make sure that PCR and CMF have been installed. See glofrim [manual](https://glofrim.readthedocs.io/en/latest/).
1. Untar the model and forcing data archives (PCR.tar.gz, CMR.tar.gz and forcing.tar.gz)
2. Change /path/to/libcama.so in glofrim_PCR2CMF.ini
3. Set an absolute input (inputDir) and output directory (outputDir) in PCR/setup_PCR_30min_Amazon.ini 
4. Compile CMF/generate_inpmat.F

HOW TO RUN
---------- 
There are two options to run the coupled PCR>CMF model:
1. Using glofrim_runner.py in glofrim/script folder. See example in run_glofrim.sh shell script. 
2. Using provided python notebooks

The output will be saved in the OUT folder in this directory.