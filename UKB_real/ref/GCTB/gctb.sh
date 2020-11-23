#!/usr/bin/bash
#SBATCH --partition pi_zhao,general,scavenge,bigmem
#SBATCH -c 1
#SBATCH --mem 70g
#SBATCH --time 1-00:00:00
#SBATCH --job-name gctb

i=${SLURM_ARRAY_TASK_ID}

gctb --bfile ../../genotype/ref5000/Ukb_imp_v2_hm3_5k \
    --gen-map ../../genotype/ref5000/Ukb_imp_v2_hm3_5k.map \
    --make-shrunk-ldm \
    --chr ${i} \
    --out ./Ukb_hm3_5k_chr${i} 

gctb --ldm Ukb_hm3_5k_chr${i}.ldm.shrunk --make-sparse-ldm --chisq 0 --out Ukb_hm3_5k_chr${i}


