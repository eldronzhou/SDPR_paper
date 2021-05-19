#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 3
#SBATCH --mem 8g
#SBATCH --time 1-00:00:00
#SBATCH --job-name lassosum

i=${SLURM_ARRAY_TASK_ID}
j=Scene4
k=50K
h2=0.5

module load R

mkdir -p sim_${i}/

Rscript lasso_pipeline.R ${i} 0.5 ${j} ${k}

mkdir -p ../predict/lassosum

Rscript cv.R ${i} ${h2} ${j}

rm -rf sim_${i}
