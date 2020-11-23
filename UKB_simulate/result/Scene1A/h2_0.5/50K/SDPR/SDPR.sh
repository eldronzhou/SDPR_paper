#!/usr/bin/bash
#SBATCH --partition pi_zhao,general,scavenge
#SBATCH -c 3
#SBATCH --mem 15g
#SBATCH --time 1-00:00:00
#SBATCH --job-name SDPR

i=${SLURM_ARRAY_TASK_ID}
h2=0.5
sz=50K
k=Scene1A

source activate python2

python ../../../../../software/SDPR.py --ss ../../../../../summary_stat/${k}/h2_${h2}/${sz}/sim_${i}.PHENO1.glm.linear --load_ld ../../../../../ref/SDPR/validate_10K_r300_full_r.1.gz --N 50000 --M 1000 --mcmc_samples 1000 --burn 200 --VS True --out sim_${i}

module load PLINK/1.90-beta5.3

mkdir -p ../predict/SDPR/
plink --bfile ../../../../../genotype/test/test_10K --score sim_${i}.txt 1 2 4 header --out ../predict/SDPR/sim_${i}
