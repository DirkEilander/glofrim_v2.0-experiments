All model and forcing files for a coupled PCR - CMF run are provided in this folder. 

BEFORE RUNNING
--------------
0. Make sure glofrim, PCR and CMF have been installed, see ../../README.md
1. Change /path/to/libcama.so in glofrim_PCR2CMF.ini
2. Set an absolute input (inputDir) and output directory (outputDir) in PCR/setup_PCR_30min_Amazon.ini 
3. Compile CMF/generate_inpmat.F

HOW TO RUN
----------
4. option a) Using glofrim_runner.py in glofrim/script folder. See example in run_glofrim.sh shell script; option b) using provided python notebooks