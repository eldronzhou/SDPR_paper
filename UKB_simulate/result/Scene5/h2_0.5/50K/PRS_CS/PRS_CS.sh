#!/usr/bin/bash
#SBATCH --partition pi_zhao,general,scavenge
#SBATCH -c 3
#SBATCH --mem 2g
#SBATCH --time 1-00:00:00
#SBATCH --job-name PRS_cS

i=${SLURM_ARRAY_TASK_ID}
j=Scene5
k=50K

source activate python2

for p in {1..10}; do
#python2 ~/software/PRScs/PRScs.py --ref_dir=/ysm-gpfs/pi/zhao/gz222/UKB_simulate/ref/PRS_CS/ldblk_1kg_eur --bim_prefix=../../../../../genotype/test/test_10K --sst_file=../../../../../summary_stat/${j}/h2_0.5/${k}/sim_${p}_prscs.txt --n_gwas=50000 --beta_std=True --chrom ${i} --out_dir=./sim_${p}

for phi in 1e-6 1e-4 1e-2 1;
do
mkdir -p ${phi}
#python2 ~/software/PRScs/PRScs.py --ref_dir=/ysm-gpfs/pi/zhao/gz222/UKB_simulate/ref/PRS_CS/ldblk_1kg_eur --bim_prefix=../../../../../genotype/test/test_10K --sst_file=../../../../../summary_stat/${j}/h2_0.5/${k}/sim_${p}_prscs.txt --n_gwas=50000 --beta_std=True --phi=${phi} --chrom ${i} --out_dir=${phi}/sim_${p}
done
done

cat sim_${i}_pst_eff_a1_b0.5_*_chr*.txt > sim_${i}.txt
cat 1e-6/sim_${i}_pst_eff_a1_b0.5_*_chr*.txt > 1e-6/sim_${i}.txt
cat 1e-4/sim_${i}_pst_eff_a1_b0.5_*_chr*.txt > 1e-4/sim_${i}.txt
cat 1e-2/sim_${i}_pst_eff_a1_b0.5_*_chr*.txt > 1e-2/sim_${i}.txt
cat 1/sim_${i}_pst_eff_a1_b0.5_*_chr*.txt > 1/sim_${i}.txt
module load R
mkdir -p ../predict/PRS_CS
h2=0.5
Rscript cv.R ${i} ${h2} ${j}
