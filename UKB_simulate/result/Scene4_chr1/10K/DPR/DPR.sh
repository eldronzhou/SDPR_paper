#!/usr/bin/bash
#SBATCH --partition scavenge
#SBATCH -c 1
#SBATCH --mem 15g
#SBATCH --time 1-00:00:00
#SBATCH --job-name DPR

i=${SLURM_ARRAY_TASK_ID}

DPR -bfile ~/scratch60/tmp/sim_${i} -dpr 2 -nk 4 -w 2000 -s 4000 -o sim_${i}

DPR -bfile test/test_10K -epm output/sim_${i}.param.txt -emu output/sim_${i}.log.txt -predict -outdir ../predict/DPR -o sim_${i}
