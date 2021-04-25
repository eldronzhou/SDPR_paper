#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 1
#SBATCH --mem 20g
#SBATCH --time 1-00:00:00
#SBATCH --job-name GCTA-simu

i=${SLURM_ARRAY_TASK_ID}

# iterate over j=Scene1A, Scene1B, Scene1C, Scene4
# here is an example for Scene1A
j=Scene1A

for h2 in 0.5; do

mkdir -p ./phenotype/${j}/h2_${h2}/validate
gcta64 --bfile ./genotype/validate/validate_10K --simu-qt --simu-causal-loci ./effect_size/${j}/h2_${h2}/sim_${i}.txt --simu-hsq ${h2} --simu-rep 1 --out ./phenotype/${j}/h2_${h2}/validate/sim_${i}

mkdir -p ./phenotype/${j}/h2_${h2}/test
gcta64 --bfile ./genotype/test/test_10K --simu-qt --simu-causal-loci ./effect_size/${j}/h2_${h2}/sim_${i}.txt --simu-hsq ${h2} --simu-rep 1 --out  ./phenotype/${j}/h2_${h2}/test/sim_${i}

for sz in 10K 50K 100K; do

mkdir -p ./phenotype/${j}/h2_${h2}/discover/${sz}
gcta64 --bfile ./genotype/discover/${sz}/discover_${sz} --simu-qt --simu-causal-loci ./effect_size/${j}/h2_${h2}/sim_${i}.txt --simu-hsq ${h2} --simu-rep 1 --out  ./phenotype/${j}/h2_${h2}/discover/${sz}/sim_${i}

module load PLINK # require PLINK2

mkdir -p summary_stat/${j}/h2_${h2}/${sz}/
plink2 --bfile ./genotype/discover/${sz}/discover_${sz} --pheno ./phenotype/${j}/h2_${h2}/discover/${sz}/sim_${i}.phen --glm a0-ref cols=chrom,pos,alt,ref,a1freq,orbeta,se,tz,p,nobs --out summary_stat/${j}/h2_${h2}/${sz}/sim_${i}

module load R # require R

# convert to COJO format for SBayesR
awk '{print $3,$6,$4,$7,$9,$10,$12,$8}' summary_stat/${j}/h2_${h2}/${sz}/sim_${i}.PHENO1.glm.linear > summary_stat/${j}/h2_${h2}/${sz}/sim_${i}_cojo.txt

# convert to PRS-CS format for PRS-CS and SDPR
awk 'BEGIN{ print "SNP","A1","A2","BETA","P"} NR>1 {print $3,$6,$4,$9,$12}' summary_stat/${j}/h2_${h2}/${sz}/sim_${i}.PHENO1.glm.linear > summary_stat/${j}/h2_${h2}/${sz}/sim_${i}_prscs.txt

# convert to GEMMA format for DBSLMM
awk 'BEGIN {OFS="\t"}; NR>1 {print $1,$3,$2,$8,$8,$6,$4,$7,$9,$10,$11}' summary_stat/${j}/h2_${h2}/${sz}/sim_${i}.PHENO1.glm.linear > summary_stat/${j}/h2_${h2}/${sz}/sim_${i}_gemma.assoc.txt

done
done
