
library(stringr)
library(lassosum)
library(data.table)
library(parallel)

args <- commandArgs(trailingOnly = T)
j = args[1] # repeat
h2 = args[2] # h2
k = args[3] # scene
sz = args[4] # sample size

ss = fread(paste0("../../../../../summary_stat/", k, "/h2_", h2, 
		  "/", sz, "/sim_", j, ".PHENO1.glm.linear"))

ref.bfile <- "../../../../../genotype/validate/validate_10K"

LDblocks <- "EUR.hg19"

cor <- p2cor(p=ss$P, n=10e3, sign=ss$BETA)

idx = which(is.na(cor))
if (length(idx) != 0 ) {
    cor = cor[-idx]
    ss = ss[-idx.]
}

s = c(0.2, 0.5, 0.9, 1)

lambda = exp(seq(log(0.001), log(0.1), length.out=20))

cl <- makeCluster(3)

out <- lassosum.pipeline(cor=cor, chr=ss$"#CHROM", pos=ss$POS,
	A1=ss$A1, ref.bfile=ref.bfile, LDblocks=LDblocks,
	s=s, lambda=lambda, cluster=cl)

for (i in seq_along(s)) {
    for (p in seq_along(lambda)) {
	res = data.frame(SNP=ss$ID, A1=ss$A1,
			 BETA=out$beta[[i]][,p])
	if (!all(res$BETA == 0)) {
		write.table(res, file=paste0("sim_",j,"/s_",s[i],"_lambda_",lambda[p],".txt"),
		row.names=F, col.names=T, quote=F, append=F)
	}
    }
}


