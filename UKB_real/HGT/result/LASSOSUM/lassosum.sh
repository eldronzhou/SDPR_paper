#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 3
#SBATCH --mem 8g
#SBATCH --time 1-00:00:00
#SBATCH --job-name lassosum


module load R

Rscript lasso_pipeline.R

mkdir -p ../predict/lassosum

module load PLINK

plink2 --bfile /ysm-gpfs/pi/zhao/gz222/UKB_real/genotype/Ukb_imp_v2_hm3 --score lassosum.txt 1 2 --score-col-nums 3-82 --threads 3 --memory 5120 --out ../predict/lassosum/lassosum

gzip ../predict/lassosum/lassosum



