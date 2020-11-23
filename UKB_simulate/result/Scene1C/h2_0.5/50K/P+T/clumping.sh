#!/usr/bin/bash
#SBATCH --partition general
#SBATCH -c 1
#SBATCH --mem 1g
#SBATCH --time 01:00:00
#SBATCH --job-name clumping

j=${SLURM_ARRAY_TASK_ID}
k=Scene1C
h2=0.5
sz=50K

module load PLINK/1.90-beta5.3
module load R

mkdir -p sim_${j}

for i in {0.01,0.05,0.1,0.2,0.4,0.6,0.8}; do
 
plink --bfile ../../../../../genotype/validate/validate_10K --clump ../../../../../summary_stat/${k}/h2_${h2}/${sz}/sim_${j}.PHENO1.glm.linear --clump-snp-field ID --clump-p1 1 --clump-p2 1 --clump-kb 250 --clump-r2 ${i} --out sim_${j}/plink_r2_${i};
awk 'NR > 1 {print $3}' sim_${j}/plink_r2_${i}.clumped > sim_${j}/clumped_r2_${i}.snplist
rm -f sim_${j}/plink_r2_${i}.clumped

plink --bfile ../../../../../genotype/validate/validate_10K --extract sim_${j}/clumped_r2_${i}.snplist --score ../../../../../summary_stat/${k}/h2_${h2}/${sz}/sim_${j}.PHENO1.glm.linear 3 6 9 header --q-score-range range.txt  ../../../../../summary_stat/${k}/h2_${h2}/${sz}/sim_${j}.PHENO1.glm.linear 3 12 header --out sim_${j}/r2_${i}

rm -f sim_${j}/*.nopred sim_${j}/*.nosex
done

mkdir -p ../predict/P+T
Rscript cv.R ${j} ${h2} ${sz} ${k}
rm -rf sim_${j}/



