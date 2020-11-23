#!/usr/bin/R

snp = read.table("../../snp_list/Ukb_hm3_1kG_noMHC.snplist")
h2s = c(0.5); M = length(snp[,1])
group = c(100, 100, 1000)
pi = c(group, M-sum(group)) / M

set.seed(1234)
n_rep = 10
for (h2 in h2s) {
const = (h2/M) / sum( pi[1:3]*c(1,.1,.01))
sigma2 = const*c(1,.1, .01)

for (i in 1:n_rep) {
repeat{
id = sample(1:4, size=M, replace=T, prob=pi)
beta = rep(0, M)
beta[id == 1] = rnorm(sum(id == 1), 0, sqrt(sigma2[1]))
beta[id == 2] = rnorm(sum(id == 2), 0, sqrt(sigma2[2]))
beta[id == 3] = rnorm(sum(id == 3), 0, sqrt(sigma2[3]))
if (abs(M*var(beta) - h2) < 1e-3) break
}
res = data.frame(snp, beta)
res = res[res$beta != 0, ]
write.table(res, file=paste0("h2_",h2,"/sim_",i,".txt"), append=F, quote=F, sep="\t", row.names=F,
            col.names=F)
}
}
