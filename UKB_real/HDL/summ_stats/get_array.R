a = read.table("SDPR.txt", header=T, stringsAsFactors=F)

b = read.table("Mc_HDL.txt", header=T, stringsAsFactors=F)

# array 2 specific
idx2 = which(a$SNP %in% b$rsid & a$N<150e3)

# array 1 specific
idx1 = which(!a$SNP %in% b$rsid & a$N<95e3)

# typed on both array
a$ARRAY = 0

# array 1 specific
a[idx1,]$ARRAY = 1

# array 2 specific
a[idx2,]$ARRAY = 2

write.table(a, file="SDPR.txt", row.names=F, col.names=T, quote=F, appen=F, sep="\t")
