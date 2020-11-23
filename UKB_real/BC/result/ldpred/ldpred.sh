#!/usr/bin/bash
#SBATCH --partition pi_zhao
#SBATCH -c 1
#SBATCH --mem 20g
#SBATCH --time 2-00:00:00
#SBATCH --job-name ldpred

module load PLINK/1.90-beta5.3

python ~/software/ldpred/LDpred.py coord --ssf ../../summ_stats/ldpred.txt --N 227688 --beta --gf ../../../ref/1000G/eur_SNPmaf5_nomhc --ssf-format CUSTOM --chr CHR --A1 A1 --A2 A2 --rs SNP --eff b --pval P --pos POS --out coord

for p in {1e-5,3e-5,1e-4,3e-4,0.001,0.003,0.01,0.03,0.1,0.3,1}; do 
python ~/software/ldpred/LDpred.py gibbs --cf coord --ldr 300 --ldf ld --out ${p} --f ${p} --N 227688
done

rm -f coord

mkdir -p ../predict/ldpred
for i in `ls *.txt`; do 
    plink --bfile ../../../genotype/Ukb_imp_v2_hm3 --score ${i} 3 4 7 header --out ../predict/ldpred/${i}; 
done 

