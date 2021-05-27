#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 3
#SBATCH --mem 4g
#SBATCH --time 5:00:00
#SBATCH --job-name SDPR

i=${SLURM_ARRAY_TASK_ID}

~/SDPR/SDPR -ref_dir ~/SDPR/ref/ -valid /ysm-gpfs/pi/zhao/gz222/UKB_real/genotype/Ukb_imp_v2_hm3.bim -ss ../../summ_stats/SDPR.txt -opt_llk 2 -N 94288 -chr ${i} -out ./res_${i}.txt -a 0.1 -thin 1 -n_threads 3 -mcmc
