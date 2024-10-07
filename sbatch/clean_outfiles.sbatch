#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH -N 1
#SBATCH -p shared
#SBATCH -A idpp
#SBATCH -J outfiles

# python
module load python/miniconda3.9
source /share/apps/python/miniconda3.9/etc/profile.d/conda.sh
conda deactivate
conda deactivate
conda deactivate
conda activate idpp

# move outfiles to separate directory
for dir in /qfs/projects/idpp/batches_for_prediction/[^_]*/; do mv ${dir%*/}/*.out /qfs/projects/idpp/batches_for_prediction/_outfiles/.; done

sacct --state RUNNING --format=jobid,start,end,CPUTImeRaw,NodeList,ReqCPUS,ReqMem,Elapsed,Timelimit