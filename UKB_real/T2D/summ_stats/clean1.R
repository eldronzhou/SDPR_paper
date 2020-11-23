a = read.table("Diagram_1000G_Diabetes_2017.txt", header=T, stringsAsFactors=F)
bim = read.table("/ysm-gpfs/pi/zhao/gz222/1000g_phase3/genotype_1KG_eur_SNPmaf5/hm3/eur_SNPmaf5_nomhc.bim", header=F, stringsAsFactors=F)
bim$Chr.Position = paste0(bim$V1, ":", bim$V4)
a = dplyr::left_join(a, bim)
a = a[!is.na(a$V1),]
a = a[,c("V2","Allele1","Allele2","Effect","Pvalue","TotalSampleSize")]
colnames(a) = c("SNP","A1","A2","log_odds","P","N")
write.table(a, file="clean1.txt", row.names=F, col.names=T, quote=F, append=F)

system("source activate ldsc")

system("python ~/software/ldsc/munge_sumstats.py --sumstats clean1.txt --out clean")

system("~/software/ldsc/ldsc.py --h2 clean.sumstats.gz --ref-ld-chr /ysm-gpfs/pi/zhao/gz222/prs_comparison/LDSC/eur_w_ld_chr/ --out ./clean_h2 --w-ld-chr /ysm-gpfs/pi/zhao/gz222/prs_comparison/LDSC/eur_w_ld_chr/")

# prepare for PRS-CS
b = read.table("clean.sumstats.gz", header=T, stringsAsFactors=F)
snp = read.table("~/scratch60/UKB_real/snp_list/1kgma5_prscs_inter.txt", stringsAsFactors=F)
b = b[b$SNP %in% snp$V1,]
b$P = 2*pnorm(-abs(b$Z))
colnames(b) = c("SNP","A1","A2","BETA","N","P")
write.table(b[,c(1,2,3,4,6)], file="PRS_cs.txt", append=F, sep="\t", quote=F, row.names=F, col.names=T)

# prepare for LDpred
bim = read.table("../../ref/1000G/eur_SNPmaf5_nomhc.bim", header=F, stringsAsFactors=F)
b = dplyr::left_join(b, bim, by=c("SNP"="V2"))
tmp = b[,c(7,9,1:4,6)]
colnames(tmp)[1:2] = c("CHR", "POS")
a = a[,c("SNP","log_odds")]
tmp = dplyr::left_join(tmp, a, by="SNP")
table(is.na(tmp[,"log_odds"]))

write.table(tmp[,c("CHR","POS","SNP","A1","A2","log_odds","P")], file="ldpred.txt", append=F, sep="\t", quote=F, row.names=F, col.names=T)

system("rm -rf clean1.txt")

print(median(b$N))



