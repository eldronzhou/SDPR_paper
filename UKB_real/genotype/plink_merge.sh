#!/usr/bin/bash
#SBATCH --partition pi_zhao,general 
#SBATCH -c 1
#SBATCH --mem 40g
#SBATCH --time 2-00:00:00
#SBATCH --job-name Plink


plink --bfile path_to_your_UKB --keep sample.list --extract Ukb_v2_hm3.txt --make-bed --out Ukb_imp_v2_hm3