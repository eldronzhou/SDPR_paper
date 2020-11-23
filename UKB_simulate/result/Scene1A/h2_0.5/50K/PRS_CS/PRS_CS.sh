#!/usr/bin/bash
#SBATCH --partition pi_zhao,general,scavenge
#SBATCH -c 3
#SBATCH --mem 3g
#SBATCH --time 2-00:00:00
#SBATCH --job-name PRS_cS

i=${SLURM_ARRAY_TASK_ID}
j=Scene1A
k=50K

source activate python2

# change your path of your PRSCS
for p in {1..10}; do
python2 ~/software/PRScs/PRScs.py --ref_dir=../../../../../ref/PRS_CS/ldblk_1kg_eur --bim_prefix=../../../../../genotype/test/test_10K --sst_file=../../../../../summary_stat/${j}/h2_0.5/${k}/sim_${p}_cojo.txt --n_gwas=50000 --beta_std=True --chrom ${i} --out_dir=./sim_${p}

for phi in 1e-6 1e-4 1e-2 1;
do
mkdir -p ${phi}
python2 ~/software/PRScs/PRScs.py --ref_dir=../../../../../ref/PRS_CS/ldblk_1kg_eur --bim_prefix=../../../../../genotype/test/test_10K --sst_file=../../../../../summary_stat/${j}/h2_0.5/${k}/sim_${p}_cojo.txt --n_gwas=50000 --beta_std=True --phi=${phi} --chrom ${i} --out_dir=${phi}/sim_${p}
done
done

