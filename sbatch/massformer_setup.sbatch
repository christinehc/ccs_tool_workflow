#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -p shared
#SBATCH -A idpp
#SBATCH -J massformer

# python
module load python/miniconda3.9
source /share/apps/python/miniconda3.9/etc/profile.d/conda.sh
conda deactivate
conda deactivate
conda deactivate

module load cuda/11.3
module load gcc/9.4.0
conda activate massformer-gpu

# create base env
cd /qfs/projects/idpp/repos/massformer
pip install -I -e
bash download_scripts/download_public.sh

sacct --state RUNNING --format=jobid,start,end,CPUTImeRaw,NodeList,ReqCPUS,ReqMem,Elapsed,Timelimit