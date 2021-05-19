#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 3
#SBATCH --mem 10g
#SBATCH --time 1-00:00:00
#SBATCH --job-name ldpred2

module load R

export OPENBLAS_NUM_THREADS=1

mkdir -p tmp-data/

#Rscript ldpred2_pipeline.R 

rm -rf tmp-data/*.sbk

mkdir -p ../predict/ldpred2

module load PLINK

plink2 --bfile /ysm-gpfs/pi/zhao/gz222/UKB_real/genotype/Ukb_imp_v2_hm3 --score ldpred2.txt 1 2 --score-col-nums 3-67 --threads 3 --memory 5120 --out ../predict/ldpred2/ldpred2

gzip ../predict/ldpred2/*.sscore
