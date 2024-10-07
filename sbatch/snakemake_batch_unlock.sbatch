#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH -N 1
#SBATCH -p shared
#SBATCH -A idpp
#SBATCH -J smkunlock

# python
module load python/miniconda3.9
source /share/apps/python/miniconda3.9/etc/profile.d/conda.sh
conda deactivate
conda deactivate
conda deactivate
conda activate idpp

# unlock every directory in batches
for dir in /qfs/projects/idpp/batches_for_prediction/[^_]*/; do snakemake -j 20 --scheduler greedy --latency-wait 100 --use-conda --use-envmodules --directory ${dir} --configfile config/config.yaml --cluster-config config/cluster.yaml --unlock; done

sacct --state RUNNING --format=jobid,start,end,CPUTImeRaw,NodeList,ReqCPUS,ReqMem,Elapsed,Timelimit