#!/usr/bin/bash
#SBATCH --partition pi_zhao,general,scavenge
#SBATCH -c 3
#SBATCH --mem 3g
#SBATCH --time 2-00:00:00
#SBATCH --job-name PRS_cS

i=${SLURM_ARRAY_TASK_ID}
j=Scene4
k=50K

cat sim_${i}_pst_eff_a1_b0.5_*_chr*.txt > sim_${i}.txt
cat 1e-6/sim_${i}_pst_eff_a1_b0.5_*_chr*.txt > 1e-6/sim_${i}.txt
cat 1e-4/sim_${i}_pst_eff_a1_b0.5_*_chr*.txt > 1e-4/sim_${i}.txt
cat 1e-2/sim_${i}_pst_eff_a1_b0.5_*_chr*.txt > 1e-2/sim_${i}.txt
cat 1/sim_${i}_pst_eff_a1_b0.5_*_chr*.txt > 1/sim_${i}.txt
module load R
mkdir -p ../predict/PRS_CS
h2=0.5
Rscript cv.R ${i} ${h2} ${j}

