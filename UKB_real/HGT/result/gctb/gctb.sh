#!/usr/bin/bash
#SBATCH --partition general,scavenge,bigmem
#SBATCH -c 1
#SBATCH --mem 50g
#SBATCH --time 1-00:00:00
#SBATCH --job-name gctb

gctb --sbayes R \
     --mldm ./Ukb_hm3_5k_ldm_spase.mldmlist \
     --pi 0.95,0.02,0.02,0.01 \
     --gamma 0.0,1e-4,1e-3,1 \
     --gwas-summary ../../summ_stats/gctb.ma \
     --chain-length 4000 \
     --burn-in 2000 \
     --out-freq 100 \
     --exclude-mhc \
     --rsq 0.9 \
     --p-value 0.4 \
     --out ./gctb_res


#rm -f *.parRes sim_${i}*.covRes sim_${i}*.Par sim_${i}*.CovEffects sim_${i}*.SnpEffects

mkdir -p ../predict/gctb
module load PLINK/1.90-beta5.3
plink --bfile ../../../genotype/Ukb_imp_v2_hm3 --score gctb_res.snpRes 2 5 8 header --out ../predict/gctb 

