#!/usr/bin/R

snp = read.table("../../snp_list/Ukb_hm3_1kG_noMHC.snplist")
h2s = c(0.5); M = length(snp[,1])

set.seed(1234)
n_rep = 10
for (h2 in h2s) {
sigma2 = h2/M 

for (i in 1:n_rep) {
repeat{
beta = rnorm(M ,0, sqrt(sigma2))
if (abs(M*var(beta) - h2) < 1e-3) break
}
res = data.frame(snp, beta)
res = res[res$beta != 0, ]
write.table(res, file=paste0("h2_",h2,"/sim_",i,".txt"), append=F, quote=F, sep="\t", row.names=F,
            col.names=F)
}
}
