#!/usr/bin/bash
#SBATCH --partition pi_zhao,general,scavenge
#SBATCH -c 3
#SBATCH --mem 4g
#SBATCH --time 5:00:00
#SBATCH --job-name SDPR

mkdir -p ../predict/SDPR
module load PLINK/1.90-beta5.3

cat *.txt > res.txt
plink --bfile ../../../genotype/Ukb_imp_v2_hm3 --score res.txt 1 2 3 --out ../predict/SDPR/SDPR

