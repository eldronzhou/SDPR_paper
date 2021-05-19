#!/usr/bin/R

input = commandArgs(trailingOnly=T)

freq = read.table("/ysm-gpfs/pi/zhao/gz222/UKB_simulate/genotype/discover/10K/discover_10K.frq", header=T, stringsAsFactors=F)
#freq = freq[1:5166,]
res = read.table(input[1], header=T, stringsAsFactors=F)
freq = freq[freq$SNP %in% res[,2],]
stopifnot(all.equal(freq$SNP, res$Name))
res$beta = res$A1Effect/sqrt(2*res$A1Frq*(1-res$A1Frq))
res$A1 = freq$A1
res = res[,c("Name","A1","beta")]

write.table(res, file=paste0(input[2], ".txt"), sep="\t", row.names=F, quote=F)
