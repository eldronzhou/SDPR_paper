#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 1
#SBATCH --mem 100g
#SBATCH --time 3:00:00
#SBATCH --job-name gctb

i=${SLURM_ARRAY_TASK_ID}
j=Scene1A
h2=0.5
sz=100K

gctb --sbayes R \
     --mldm ./validate_10K.mldmlist \
     --pi 0.95,0.02,0.02,0.01 \
     --gamma 0.0,0.01,0.1,1 \
     --gwas-summary ../../../../../summary_stat/${j}/h2_${h2}/${sz}/sim_${i}_cojo.txt \
     --chain-length 4000 \
     --burn-in 2000 \
     --exclude-mhc \
     --out-freq 200 \
     --out ./sim_${i}_scaled

rm -f sim_${i}*.parRes sim_${i}*.covRes sim_${i}*.Par sim_${i}*.CovEffects sim_${i}*.SnpEffects

mkdir -p ../predict/gctb/

# scale the effect size by allele frequencies
module load R
Rscript convert_basic.R sim_${i}_scaled.snpRes sim_${i}_scaled

module load PLINK/1.90-beta5.3
plink --bfile ../../../../../genotype/test/test_10K --score sim_${i}_scaled.txt 1 2 3 header --out ../predict/gctb/sim_${i};




