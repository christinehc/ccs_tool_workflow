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
conda activate idpp

# create base env
conda env remove -n graff-ms
cd /qfs/projects/idpp/repos/idpp/workflow/envs
mamba env create -f graff-ms.yml

sacct --state RUNNING --format=jobid,start,end,CPUTImeRaw,NodeList,ReqCPUS,ReqMem,Elapsed,Timelimit