
library(stringr)
library(lassosum)
library(data.table)
library(parallel)

ss = fread("../../summ_stats/ldpred.txt")
ss = ss[order(ss$CHR, ss$POS),]

ref.bfile <- "/ysm-gpfs/pi/zhao/gz222/1000g_phase3/genotype_1KG_eur_SNPmaf5/hm3/eur_SNPmaf5_nomhc"

LDblocks <- "EUR.hg19"

bim = read.table(paste0(ref.bfile, ".bim"), header=F, stringsAsFactors=F)

cor <- p2cor(p=ss$P, n=233766, sign=ss$b)

idx = which(is.na(cor))

if (length(idx) != 0) {
    cor = cor[-idx]
    ss = ss[-idx,]
}

s = c(0.2, 0.5, 0.9, 1)

lambda = exp(seq(log(0.001), log(0.1), length.out=20))

cl <- makeCluster(3)

out <- lassosum.pipeline(cor=cor, chr=ss$CHR, pos=ss$POS,
	A1=ss$A1, A2=ss$A2, ref.bfile=ref.bfile, LDblocks=LDblocks,
	s=s, lambda=lambda, cluster=cl)

for (i in seq_along(s)) {
    if (i == 1) {
	beta = out$beta[[1]]
    }
    else {
	beta = cbind(beta, out$beta[[i]])
    }
}

res = data.frame(SNP=ss[out$sumstats$order,]$SNP, 
		 A1=out$sumstats$A1,
		 BETA=beta)

write.table(res, file="lassosum.txt",
	row.names=F, col.names=T, quote=F, append=F)




