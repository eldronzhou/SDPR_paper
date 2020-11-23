#!/usr/bin/bash
#SBATCH --partition pi_zhao,general,scavenge
#SBATCH -c 3
#SBATCH --mem 3g
#SBATCH --time 1-00:00:00
#SBATCH --job-name PRS_cS

cat *.txt > PRS_cs.txt
cat 1e-6/*.txt > 1e-6/PRS_cs.txt
cat 1e-4/*.txt > 1e-4/PRS_cs.txt
cat 1e-2/*.txt > 1e-2/PRS_cs.txt
cat 1/*.txt > 1/PRS_cs.txt

module load PLINK/1.90-beta5.3
mkdir -p ../predict/PRS_CS/
plink --bfile ../../../genotype/Ukb_imp_v2_hm3 --score PRS_cs.txt 2 4 6 --out ../predict/PRS_CS/PRS_CS

for phi in 1e-6 1e-4 1e-2 1
do
plink --bfile ../../../genotype/Ukb_imp_v2_hm3 --score ${phi}/PRS_cs.txt 2 4 6 --out ../predict/PRS_CS/${phi}
done


