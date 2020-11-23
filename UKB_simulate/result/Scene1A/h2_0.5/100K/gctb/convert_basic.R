#!/usr/bin/R

input = commandArgs(trailingOnly=T)

freq = read.table("../../../../../genotype/discover/100K/discover_100K.frq", header=T, stringsAsFactors=F)
res = read.table(input[1], header=T, stringsAsFactors=F)
freq = freq[freq$SNP %in% res[,2],]
stopifnot(all.equal(freq$SNP, res$Name))
res$beta = res$Effect/sqrt(2*res$GeneFrq*(1-res$GeneFrq))
res$A1 = freq$A1
res = res[,c("Name","A1","beta")]

write.table(res, file=paste0(input[2], ".txt"), sep="\t", row.names=F, quote=F)
