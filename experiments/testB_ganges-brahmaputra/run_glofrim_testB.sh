export OMP_NUM_THREADS=10
# activate python conda environment based on the environment.yml file
source activate glofrim_exp

# SET path to glofrim/scripts folder
# GLOFRIM_PATH=/path/to/glofrim/scripts
GLOFRIM_PATH=/home/dirk/repos/glofrim/scripts

# create output directories
mkdir -p ./OUT
mkdir -p ./OUT/PCR
mkdir -p ./OUT/CMF
mkdir -p ./OUT/LFP

# run glofrim
python "$GLOFRIM_PATH"/glofrim_runner.py run glofrim_PCR2CMF2LFP.ini -s 2001-01-01 -e 2009-12-31