#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 3
#SBATCH --mem 5g
#SBATCH --time 1-00:00:00
#SBATCH --job-name dblsmm

module load R
module load PLINK/1.90-beta5.3

dbslmm=~/software/DBSLMM/dbslmm
DBSLMM=~/software/DBSLMM/software/DBSLMM.R
m=`cat ../../summ_stats/gemma.assoc.txt | wc -l` 
blockf=~/software/DBSLMM/block/EUR/fourier_ls-all


for pv in 1e-05 1e-06 1e-07 1e-08
do
    for r2 in 0.05 0.1 0.15 0.2 0.25
    do
	Rscript ${DBSLMM} --summary ./gemma.assoc.txt --outPath ./ --plink plink --dbslmm ${dbslmm} --ref /ysm-gpfs/pi/zhao/gz222/1000g_phase3/genotype_1KG_eur_SNPmaf5/hm3/eur_SNPmaf5_nomhc --pv ${pv} --r2 ${r2} --mafMax 1 --n 227688 --nsnp ${m} --type d --block ${blockf}.bed --h2 0.13 --thread 3    	
    	mv gemma.dbslmm.txt pv_${pv}_r2_${r2}.txt
    done
done

mkdir -p ../predict/dbslmm

module load PLINK/1.90-beta5.3
for i in `ls *.txt`; do 
        plink --bfile /ysm-gpfs/pi/zhao/gz222/UKB_real/genotype/Ukb_imp_v2_hm3 --score ${i} 1 2 3 header --out ../predict/dbslmm/${i} 
done

gzip ../predict/dbslmm/*.profile

