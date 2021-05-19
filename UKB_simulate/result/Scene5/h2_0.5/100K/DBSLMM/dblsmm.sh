#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 3
#SBATCH --mem 2g
#SBATCH --time 1-00:00:00
#SBATCH --job-name dblsmm

i=${SLURM_ARRAY_TASK_ID}
j=Scene5
k=100K
h2=0.5

module load R
module load PLINK/1.90-beta5.3

dbslmm=~/software/DBSLMM/dbslmm
DBSLMM=~/software/DBSLMM/software/DBSLMM.R
m=681828
blockf=~/software/DBSLMM/block/EUR/fourier_ls-all

mkdir -p sim_${i}/

for pv in 1e-05 1e-06 1e-07 1e-08
do
    for r2 in 0.05 0.1 0.15 0.2 0.25
    do
	Rscript ${DBSLMM} --summary ../../../../../summary_stat/${j}/h2_0.5/${k}/sim_${i}_gemma.assoc.txt --outPath sim_${i}/ --plink plink --dbslmm ${dbslmm} --ref ../../../../../genotype/validate/validate_10K --pv ${pv} --r2 ${r2} --mafMax 1 --n 100000 --nsnp ${m} --type d --block ${blockf}.bed --h2 0.5 --thread 3
    	mv sim_${i}/sim_${i}_gemma.dbslmm.txt sim_${i}/sim_${i}_pv_${pv}_r2_${r2}.txt
    done
done

mkdir -p ../predict/dblsmm
Rscript cv.R ${i} ${h2} ${j}

rm -rf sim_${i}
