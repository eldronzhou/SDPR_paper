
library(data.table)

a = fread("../../summ_stats/gemma.assoc.txt")

b = fread("/ysm-gpfs/pi/zhao/gz222/1000g_phase3/genotype_1KG_eur_SNPmaf5/hm3/eur_SNPmaf5.bim")

a = a[order(a$V1, a$V3),]

idx = match(a$V2, b$V2)
 
flip = which(a$V6 == b[idx,]$V6 & a$V7 == b[idx,]$V5)
a[flip,9] = -a[flip,9]
tmp = a[flip,6]
a[flip,6] = a[flip,7]
a[flip,7] = tmp

write.table(a, file="gemma.assoc.txt", row.names=F, col.names=F, quote=F, append=F, sep="\t")
