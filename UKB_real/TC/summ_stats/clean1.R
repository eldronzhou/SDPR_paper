a = read.table("jointGwasMc_TC.txt", header=T, stringsAsFactors=F)
a = a[!duplicated(a$SNP_hg19),]
rm_idx = which((a$beta/a$se)^2 > 80)
a = a[-rm_idx,]
a$P = 2*pnorm(-abs(a$beta/a$se))

write.table(a[,c("rsid","A1","A2","beta","P","N")], file="clean1.txt", row.names=F, col.names=T, quote=F, append=F)

system("source activate ldsc")

system("python ~/software/ldsc/munge_sumstats.py --sumstats clean1.txt --out clean")

system("~/software/ldsc/ldsc.py --h2 clean.sumstats.gz --ref-ld-chr /ysm-gpfs/pi/zhao/gz222/prs_comparison/LDSC/eur_w_ld_chr/ --out ./clean_h2 --w-ld-chr /ysm-gpfs/pi/zhao/gz222/prs_comparison/LDSC/eur_w_ld_chr/")

# prepare for PRS-CS
b = read.table("clean.sumstats.gz", header=T, stringsAsFactors=F)
snp = read.table("../../snp_list/1kgma5_prscs_inter.txt", stringsAsFactors=F)
b = b[b$SNP %in% snp$V1,]
b$P = 2*pnorm(-abs(b$Z))
colnames(b) = c("SNP","A1","A2","BETA","N","P")
write.table(b[,c(1,2,3,4,6)], file="PRS_cs.txt", append=F, sep="\t", quote=F, row.names=F, col.names=T)

# prepare for LDpred
bim = read.table("../../ref/1000G/eur_SNPmaf5_nomhc.bim", header=F, stringsAsFactors=F)
b = dplyr::left_join(b, bim, by=c("SNP"="V2"))
tmp = b[,c(7,9,1:4,6)]
colnames(tmp)[1:2] = c("CHR", "POS")
a = a[,c("rsid","beta")]
tmp = dplyr::left_join(tmp, a, by=c("SNP"="rsid"))
table(is.na(tmp[,"beta"]))

write.table(tmp[,c("CHR","POS","SNP","A1","A2","beta","P")], file="ldpred.txt", append=F, sep="\t", quote=F, row.names=F, col.names=T)

system("rm -rf clean1.txt")

print(median(b$N))



