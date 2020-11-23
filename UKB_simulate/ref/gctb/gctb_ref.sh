#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 1
#SBATCH --mem 30g
#SBATCH --time 05:00:00
#SBATCH --job-name gctb_ref

i=${SLURM_ARRAY_TASK_ID}

gctb --bfile ../../genotype/validate/validate_10K --chr ${i} --make-full-ldm --out ./validate_10K_chr${i}

