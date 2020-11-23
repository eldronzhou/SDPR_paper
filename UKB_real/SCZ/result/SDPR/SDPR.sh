#!/usr/bin/bash
#SBATCH --partition pi_zhao,general,scavenge
#SBATCH -c 3
#SBATCH --mem 4g
#SBATCH --time 5:00:00
#SBATCH --job-name SDPR

i=${SLURM_ARRAY_TASK_ID}

~/SDPR/SDPR -ref_dir ../../../ref/SDOR -valid ../../../genotype/Ukb_imp_v2_hm3.bim -ss ../../summ_stats/PRS_cs.txt  -N 65955 -chr ${i} -out ./res_${i}.txt -a 0.1 -thin 1 -n_threads 3 -mcmc

