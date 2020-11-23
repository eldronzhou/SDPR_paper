#!/usr/bin/bash
#SBATCH --partition pi_zhao,general,scavenge
#SBATCH -c 3
#SBATCH --mem 3g
#SBATCH --time 1-00:00:00
#SBATCH --job-name PRS_cS

i=${SLURM_ARRAY_TASK_ID}

source activate python2

python2 ~/software/PRScs/PRScs.py --ref_dir=../../../ref/PRS_CS/ldblk_1kg_eur --bim_prefix=../../../genotype/Ukb_imp_v2_hm3 --sst_file=../../summ_stats/PRS_cs.txt --n_gwas=233766 --chrom=${i} --out_dir=./

for phi in 1e-6 1e-4 1e-2 1
do
mkdir -p ${phi}/
python2 ~/software/PRScs/PRScs.py --ref_dir=../../../ref/PRS_CS/ldblk_1kg_eur --bim_prefix=../../../genotype/Ukb_imp_v2_hm3 --sst_file=../../summ_stats/PRS_cs.txt --n_gwas=233766 --chrom=${i} --out_dir=${phi}/ --phi=${phi}

done
