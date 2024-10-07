#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -p shared
#SBATCH -A idpp
#SBATCH -J idpp

# python
module load python/miniconda3.9
source /share/apps/python/miniconda3.9/etc/profile.d/conda.sh
conda deactivate
conda deactivate
conda deactivate
module load java/1.8.0_31
module load cuda/11.1
module load gcc/9.4.0
conda activate idpp

# create base env
# conda env remove -n ms-pred
cd /qfs/projects/idpp/repos/ms-pred
# mamba env create -f environment.yml
mamba activate ms-gen
pip install -r requirements.txt
python3 setup.py develop

# Scarf quickstart
. quickstart/scarf/download_model.sh
. quickstart/scarf/run_model.sh

# Iceberg quickstart
. quickstart/iceberg/download_model.sh
. quickstart/iceberg/run_model.sh

# Download formatted gnps
. data_scripts/download_gnps.sh

# ICEBERG processing
# . data_scripts/dag/run_magma.sh

# Add run params to output log
sacct --state RUNNING --format=jobid,start,end,CPUTImeRaw,NodeList,ReqCPUS,ReqMem,Elapsed,Timelimit
