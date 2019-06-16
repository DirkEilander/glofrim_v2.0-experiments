All model and forcing files for a coupled PCR - CMF run are provided in this folder. 

BEFORE RUNNING
--------------
0. Make sure glofrim, PCR and CMF have been installed, see ../../README.md
1. Untar PCD.tar.gz and CMR.tar.gz archive
2. Change /path/to/libcama.so in glofrim_PCR2CMF.ini
3. Set an absolute input (inputDir) and output directory (outputDir) in PCR/setup_PCR_30min_Amazon.ini 
4. Compile CMF/generate_inpmat.F

HOW TO RUN
---------- 
There are two options to run the coupled PCR>CMF model:
1. Using glofrim_runner.py in glofrim/script folder. See example in run_glofrim.sh shell script. 
2. Using provided python notebooks