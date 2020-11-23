
library(stringr)

args <- commandArgs(trailingOnly = T)
j = args[1] # repeat
h2 = args[2] # h2
sz = args[3] # sample size
k = args[4] # sceneario

pheno = read.table(paste0("../../../../../phenotype/",k,"/h2_",h2,"/validate/sim_",j,".phen"), header=F)

file = list.files(path=paste0("sim_",j), pattern="*.profile", full.names=T)
r2 = rep(0, length(file))

for (i in seq_along(file)) {
	dat = read.table(file[i], header=T)
	r2[i] = cor(pheno[,3], dat$SCORE)^2
}

idx = which.max(r2)
print(paste0(file[idx], ": ", r2[idx]))

name = str_extract(file[idx], "r2_0.[0-9]*")

command = paste0("plink --bfile ../../../../../genotype/test/test_10K --extract sim_",j,"/clumped_",name,".snplist --score ../../../../../summary_stat/",k,"/h2_",h2,"/",sz,"/sim_",j,".PHENO1.glm.linear 3 6 9 header --q-score-range range.txt ../../../../../summary_stat/",k,"/h2_",h2,"/",sz,"/sim_",j,".PHENO1.glm.linear 3 12 header --out sim_",j,"/",name)

system(command, wait=T, intern=F)

command = paste0("cp ", file[idx],"  ../predict/P+T/sim_",j,".profile")
system(command, wait=T, intern=F)

