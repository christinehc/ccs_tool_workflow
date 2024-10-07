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
module load cuda/11.3
module load gcc/9.4.0
module load singularity/4.1.2
conda activate snekmer

# create base env
mamba env create -f environment.yml
conda activate idpp
pip install -r requirements.txt

# install DeepCCS
cd /qfs/projects/idpp/repos/DeepCCS/core
python setup.py install

# create GrAFF-MS conda environment
cd /qfs/projects/idpp/repos/idpp/workflow/envs
conda env create -f envs/graff-ms.yml

sacct --state RUNNING --format=jobid,start,end,CPUTImeRaw,NodeList,ReqCPUS,ReqMem,Elapsed,Timelimit