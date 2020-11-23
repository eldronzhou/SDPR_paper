#!/usr/bin/bash
#SBATCH --partition general 
#SBATCH -c 1
#SBATCH --mem 40g
#SBATCH --time 07:00:00
#SBATCH --job-name Plink


module load PLINK/1.90-beta5.3

# replace file with path to your  UKB genotype data
plink --bfile /ysm-gpfs/pi/zhao/yy496/ukb_imp/ukb_imp_qc1_set --keep test_10K.list --extract ../../snp_list/Ukb_hm3_1kG_noMHC.snplist --keep-allele-order --make-bed --out validate_10K

