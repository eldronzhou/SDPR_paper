# SDPR_paper
We are not able to release data used in the analysis because of the restricted access to the UK Biobank genotype and phenotype data. However, if you have access to UK Biobank, you can follow the instructions to reproduce the results in our paper. The path to the software and datasets in these scripts may not be correct. We use slurm to schedule jobs to HPC, and you need to change the header if you use another system. If you have issues about these scripts, please report to the [issue](https://github.com/eldronzhou/SDPR_paper/issues) page.

# Table of Contents
- [Simulations](#sim)
  - [Requirements](#sim-req)
  - [Obtaining genotype data](#sim-geno)
  - [Simulating phenotype and generating summary statistics](#sim-ss)
  - [Constructing the reference LD matrix](#sim-ref)
  - [Running the analysis](#sim-analysis)
- [Real data applications](#real)
  - [Requirements](#real-req)
  - [Obtaining genotype data](#real-geno)
  - [Constructing the reference LD matrix](#real-ref)
  - [Obtaining and cleaning the summary statistics](#real-ss)
  - [Running the analysis](#real-analysis)


# <a name="sim"></a>Simulations
## <a name="sim-req"></a>Requirements
* PRS-CS (python2.7, numpy, scipy, recommend to install anaconda)
* GCTB
* GCTA
* LDpred 
* PLINK-1.90 
* PLINK-2.0
* R
* about 340 Gb desk space 

## Workflow
**<a name="sim-geno"></a>1. Obtaining Genotype data**

We provide the list of SNPs used in our simulations. Sample list used in the simulation can be obtained upon request so that you can obtain the exact genotype data we used if you have the access to UK Biobank.

```
git clone https://github.com/eldronzhou/SDPR_paper.git
cd UKB_simulate/genotype
cd discover/10K 
sbatch plink_merge.sh # you need to change the path to the original UK Biobank genotype

# repeat the above procedure for discover/50K, discover/100K, validate/ and test/
```

**<a name="sim-ss"></a>2. Simulating phenotype and generating summary statistics**

Next run the simulation script to generate summary statistics. You need to change to the directory of `UKB_simulate/` to submit the script.
```
# for Scene1A
sbatch GCTA_sim.sh

# change the j in GCTA_sim.sh and resubmit
# for Scene1B, Scene1C and Scene4
```

**<a name="sim-ref"></a>3. Constructing the reference LD matrix**

```
# download reference of PRS_CS
cd ref/PRS_CS; sh get_ref.sh 

# download reference of SDPR
cd ../SDPR; sh get_ref.sh 

# estimate ref of gctb for each chromosome
cd ../gctb; sbatch --array=1-22 gctb.sh 
```

**<a name="sim-analysis"></a>4. Running the analysis**

We will use Scene1A as the example for demonstration. You can repeat the same procedure for Scene1B, Scene1C and Scene4.

```
cd result/Scene1A/h2_0.5/10K/

# SBayesR
cd gctb/; sbatch gctb.sh

# PRS-CS
cd ../PRS_CS/; sbatch --array=1-22 PRS_CS.sh
# after all jobs finish
sbatch --array=1-10 PRS_CS_res.sh

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

# <a name="real"></a>Real data applications

## <a name="real-req"></a>Requirements

* same requirements as simulations except for PLINK-2.0 and GCTA
* LDSC 
* Installed SDPR
* about 150 Gb desk space

## Workflow

Sample list of UKB_real and 5000 ref for SBayesR can be obtained upon request.

**<a name="real-geno"></a>1. Obtaining Genotype data**

```
cd UKB_real/genotype; sbatch plink_merge.sh

# 5000 UKB individuals for ref matrix of SBayesR
cd ref_5000/; sbatch plink_merge.sh
```

**<a name="real-ref"></a>2. Constructing the reference LD matrix**

```
cd UKB_real/ref

# get genotype of 1000G EUR
cd 1000G/; sh get_ref.sh

# get reference of PRS_CS
cd ../PRS_CS/; sh get_ref.sh

# estimate reference for SBayesR
cd ../GCTB/; sbatch --array=1-22 gctb.sh

# estimate reference for SDPR
cd ../SDPR/; sbatch --array=1-22 SDPR_ref.sh
```

**<a name="real-ss"></a>3. Obtaining and cleaning the summary statistics**

Many GWAS consortium publishes summary statistics. However, due to the data access agreement, we are not able to directly provide the original copy. You can find the study of the summary statistics in the Table 1 of the manuscript and download the summary statistics on your own. If you need assistance, feel free to submit to the issue. Here we provide an example on downloading and processing the height summary statistics from the GIANT consortium. The procedure is similar for other traits. 

```
cd UKB_real/HGT/summ_stats/

# download
wget https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz

# make sure you check the name of summary statistics is right
# you also need to change the path of LDSC, munge_summstats.py, and reference of LDSC
Rscript clean1.R
```

**<a name="real-analysis"></a>4. Running the analysis**

Here is the example on running analysis for height. The procedure is similar for other traits. 

```
cd UKB_real/HGT/result/

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
cd ../SDPR/; sbatch --array=1-22 SDPR.sh
# after all jobs finish
sbatch SDPR_res.sh
```

We provide the script to make the figure under the directory `UKB_real/HGT/result/predict`, although the input file is not available because the restricted access to the UK Biobank phenotype information. 

