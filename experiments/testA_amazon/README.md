All model and forcing files for a coupled PCR - CMF run are provided in this folder. 

BEFORE RUNNING
--------------
0. Setup environment: Create a glofrim_exp environment in conda using the environment.yml file. Make sure that PCR and CMF have been installed. See glofrim [manual](https://glofrim.readthedocs.io/en/latest/).
1. Prepare GLOFRIM ini file: Change /path/to/libcama.so in glofrim_PCR2CMF.ini
2. Prepare CMF model files: Untar the CMF/hires.tar.gz archive and compile CMF/generate_inpmat.F using the CMF/Makefile
3. Prepare PCR ini file: Set an absolute input (inputDir) and output directory (outputDir) in PCR/setup_PCR_30min_Amazon.ini 

HOW TO RUN
---------- 
1. Using glofrim_runner.py in glofrim/script folder. See example in run_glofrim_testA.sh shell script. You will need to change /path/to/glofrim/scripts 

OUTPUTS
-------
The output will be saved in the OUT folder in this directory.

NOTE
----
As the PCRGLOB-WB DynRout module is not available in the version with BMI adaptor, the run including this module was done seperately using a different ini file (PCR/DynRout_30min_Amazon_noWaterbodies_2005_2010_para.ini) compared the coupled GLOFRIM run.
