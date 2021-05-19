#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 3
#SBATCH --mem 4g
#SBATCH --time 5:00:00
#SBATCH --job-name SDPR

i=${SLURM_ARRAY_TASK_ID}
j=Scene4_chr1
k=10K


~/SDPR/SDPR -ref_dir ~/scratch60/UKB_simulate_new/ref/SDPR/ -valid ../../../../genotype/test/test_10K.bim -ss ../../../../summary_stat/${j}/h2_0.3/${k}/sim_${i}_prscs.txt -N 10000 -chr 1 -a 0 -out ./sim_${i}.txt -n_threads 3 -mcmc

module load PLINK/1.90-beta5.3
plink --bfile ../../../../genotype/test/test_10K --chr 1 --score sim_${i}.txt 1 2 3 header --out ../predict/SDPR/sim_${i}

