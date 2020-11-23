
library(stringr)

args <- commandArgs(trailingOnly = T)
j = args[1] # repeat
h2 = args[2] # h2
k = args[3] # scene

pheno = read.table(paste0("../../../../../phenotype/",k,"/h2_",h2,"/validate/sim_",j,".phen"), header=F)

file = list.files(path=paste0("sim_",j), pattern="*.txt$", full.names=T)
r2 = rep(0, length(file))

# you may want to change your path to plink
for (i in seq_along(file)) {
	command = paste0("plink --bfile ../../../../../genotype/validate/validate_10K --score ",file[i]," 3 4 7 header --out ",file[i] )
	system(command, intern=F, wait=T)
	dat = read.table(paste0(file[i],".profile"), header=T)
	r2[i] = cor(pheno[,3], dat$SCORE)^2
}

idx = which.max(r2)
print(paste0(file[idx], ": ", r2[idx]))


#name = str_extract(file[idx], "r2_0.[0-9]*")

command = paste0("plink --bfile ../../../../../genotype/test/test_10K --score ", file[idx], " 3 4 7 header --out ../predict/ldpred/sim_",j)

system(command, wait=T, intern=F)


