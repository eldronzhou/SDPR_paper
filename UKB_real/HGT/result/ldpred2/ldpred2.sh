#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 10
#SBATCH --mem 35g
#SBATCH --time 1-00:00:00
#SBATCH --job-name ldpred2

module load R

export OPENBLAS_NUM_THREADS=1

Rscript ldpred2_pipeline.R 

