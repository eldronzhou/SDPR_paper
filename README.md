# SDPR_paper
We are not able to release data used in the analysis because of the restricted access to the UK Biobank genotype and phenotype data. However, if you have access to UK Biobank, you can follow the instructions to reproduce the results in our paper. The path to the software and datasets in these scripts may not be correct. We use slurm to schedule jobs to HPC, and you need to change the header if you use another system. If you have issues about these scripts, please report to the [issue](https://github.com/eldronzhou/SDPR_paper/issues) page.

# Simulations
## Requirements
* PRS-CS (python2.7, numpy, scipy, recommend to install anaconda)
* GCTB
* GCTA
* PLINK-1.90 
* PLINK-2.0
* R
* about 340 Gb Desk space 

**1. Obtaining Genotype data**

We provide the list of SNPs and individuals used in our simulations. You can obtain the exact genotype data we used if you have access the to UK Biobank.

```
git clone https://github.com/eldronzhou/SDPR_paper.git
cd UKB_simulate/genotype
cd discover/10K 
sbatch plink_merge.sh # you need to change the path to the original UK Biobank genotype

# repeat the above procedure for discover/50K, discover/100K, validate/ and test/
```

**2. Simulating phenotype and generating summary statistics**

Next run the simulation script to generate summary statistics. You need to change to the directory of `UKB_simulate/` to submit the script.
```
# for Scene1A
sbatch GCTA_sim.sh

# change the j in GCTA_sim.sh and resubmit
# for Scene1B, Scene1C and Scene4
```

**3. Constructing the reference LD matrix**

```
# download reference of PRS_CS
cd ref/PRS_CS; sh get_ref.sh 

# download reference of SDPR
cd ../SDPR; sh get_ref.sh 

# estimate ref of gctb for each chromosome
cd ../gctb; sbatch --array=1-22 gctb.sh 
```

**4. Running the analysis**

We will use Scene1A as the example for demonstration. You can repeat the same procedure for Scene1B, Scene1C and Scene4.

```
cd result/Scene1A/h2_0.5/10K/

# SBayesR
cd gctb/; sbatch gctb.sh

# PRS-CS
cd ../PRS_CS/; sbatch --array=1-22 PRS_CS.sh
# after all jobs finish
sbatch PRS_CS_res.sh

# LDpred
cd ../ldpred/; sbatch ldpred.sh

# P+T
cd ../P+T/; sbatch clumping.sh

# SDPR
cd ../SDPR/; sbatch SDPR.sh

# repeat the above procedures for 50K and 100K

# after all jobs finished, make the plot
# under the directory of UKB_simulate/result/Scene1A/h2_0.5/
Rscript get_Res.R
```

# Real data applications





