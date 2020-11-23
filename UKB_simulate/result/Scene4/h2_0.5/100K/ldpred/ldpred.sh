#!/usr/bin/bash
#SBATCH --partition pi_zhao
#SBATCH -c 1
#SBATCH --mem 12g
#SBATCH --time 15:00:00
#SBATCH --job-name ldpred

i=${SLURM_ARRAY_TASK_ID}
j=Scene4
h2=0.5
sz=100K

mkdir -p sim_${i}

# you need to change path to your LDpred
python ~/software/ldpred/LDpred.py coord --ssf ../../../../../summary_stat/${j}/h2_${h2}/$sz/sim_${i}.PHENO1.glm.linear --N 100000 --beta --gf ../../../../../genotype/validate/validate_10K --ssf-format CUSTOM --chr \#CHROM --pos POS --A1 A1 --A2 REF --rs ID --eff BETA --pval P --out sim_${i}/sim_${i}_coord

for p in {1e-5,3e-5,1e-4,3e-4,0.001,0.003,0.01,0.03,0.1,0.3,1}; do 
python ~/software/ldpred/LDpred.py gibbs --cf sim_${i}/sim_${i}_coord --ldr 227 --ldf sim_${i}/ld --out sim_${i}/${p} --f ${p} --N 100000 --n-iter 100 
done

rm -f sim_${i}/sim_${i}_coord

module load R
mkdir -p ../predict/ldpred
Rscript cv.R ${i} ${h2} ${j}

rm -rf sim_${i}
