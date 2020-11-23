#!/usr/bin/bash
#SBATCH --partition pi_zhao
#SBATCH -c 1
#SBATCH --mem 10g
#SBATCH --time 2-00:00:00
#SBATCH --job-name clumping

module load PLINK/1.90-beta5.3

mkdir -p ../predict/P+T

for i in {0.2,0.4,0.6,0.8}; do
plink --bfile ../../../ref/1000G/eur_SNPmaf5_nomhc --clump ../../summ_stats/ldpred.txt --clump-kb 250 --clump-p1 1 --clump-p2 1 --clump-r2 ${i} --clump-snp-field SNP --clump-field P --out plink_r2_${i}

awk 'NR > 1 {print $3}' plink_r2_${i}.clumped > clumped_r2_${i}.snplist
rm -f plink_r2_${i}.clumped

plink --bfile ../../../genotype/Ukb_imp_v2_hm3 --extract clumped_r2_${i}.snplist --out ../predict/P+T/r2_${i} --q-score-range range.txt ../../summ_stats/ldpred.txt 3 7 header --score ../../summ_stats/ldpred.txt 3 4 6 header
rm -f *.nopred *.nosex
done




