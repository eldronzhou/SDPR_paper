#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 1
#SBATCH --mem 5g
#SBATCH --time 1-00:00:00
#SBATCH --job-name ldpred2

i=${SLURM_ARRAY_TASK_ID}
j=Scene1A
h2=0.5
k=10K

module load R

mkdir -p sim_${i}/

export OPENBLAS_NUM_THREADS=1

#Rscript ldpred2_pipeline.R ${i} 0.5 ${j} ${k}

rm -f sim_${i}/.sbk

mkdir -p ../predict/ldpred2

Rscript cv.R ${i} ${h2} ${j}
