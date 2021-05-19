
library(stringr)

args <- commandArgs(trailingOnly = T)
j = args[1] # repeat
h2 = args[2] # h2
k = args[3] # scene

pheno = read.table(paste0("/ysm-gpfs/pi/zhao/gz222/UKB_simulate/phenotype/",k,"/h2_",h2,"/validate/sim_",j,".phen"), header=F)

#file = list.files(path=paste0("sim_",j), pattern="*.txt$", full.names=T)
file = paste0(c("sim_","1e-6/sim_","1e-4/sim_","1e-2/sim_","1/sim_"), j, ".txt")
r2 = rep(0, length(file))

for (i in seq_along(file)) {
	command = paste0("module load PLINK/1.90-beta5.3; plink --bfile /ysm-gpfs/pi/zhao/gz222/UKB_simulate/genotype/validate/validate_10K --score ",file[i]," 2 4 6 --out ",file[i] )
	system(command, intern=F, wait=T)
	dat = read.table(paste0(file[i],".profile"), header=T)
	r2[i] = cor(pheno[,3], dat$SCORE)^2
}

idx = which.max(r2)
print(paste0(file[idx], ": ", r2[idx]))


#name = str_extract(file[idx], "r2_0.[0-9]*")

command = paste0("module load PLINK/1.90-beta5.3; plink --bfile /ysm-gpfs/pi/zhao/gz222/UKB_simulate/genotype/test/test_10K --score ", file[idx], " 2 4 6  --out ../predict/PRS_CS/sim_",j)

system(command, wait=T, intern=F)


