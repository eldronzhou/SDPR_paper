
library(stringr)
library(data.table)
library(bigsnpr)
library(methods)

args <- commandArgs(trailingOnly = T)
j = args[1] # repeat
h2 = args[2] # h2
k = args[3] # scene
sz = args[4] # sample size

ss = fread(paste0("../../../../../summary_stat/", k, "/h2_", h2, 
		  "/", sz, "/sim_", j, ".PHENO1.glm.linear"))

ukb = snp_attach("../../../../../ref/ldpred2/validate_10K.rds")
G = ukb$genotypes

for (chr in 1:22) {
    # read corr
    load(paste0("../../../../../ref/ldpred2/UKB_ldpred2_chr", chr, ".rds"))

    ind.chr = which(ss$"#CHROM" == chr)

    if (chr == 1) {
	df_beta = ss[ind.chr, c("BETA","SE","OBS_CT")]
	ld = Matrix::colSums(corr^2)
	corr1 = as_SFBM(corr, paste0("sim_",j,"/"))
    }
    else {
	df_beta = rbind(df_beta, ss[ind.chr, c("BETA","SE","OBS_CT")])
	ld = c(ld, Matrix::colSums(corr^2))
	corr1$add_columns(corr, nrow(corr1))
    }
}

colnames(df_beta) = c("beta","beta_se","n_eff")

(ldsc = with(df_beta, snp_ldsc(ld, length(ld), chi2=(beta/beta_se)^2, 
			       sample_size = n_eff, blocks = NULL)))

h2_est <- ldsc[["h2"]]

# LDpred2-inf
beta_inf = snp_ldpred2_inf(corr1, df_beta, h2_est)

res = data.frame(SNP=ss$ID, A1=ss$A1, BETA=beta_inf)

write.table(res, file=paste0("sim_",j,"/ldpred_inf.txt"), row.names=F, col.names=F, quote=F, append=F)

# LDpred2-grid
#(h2_seq = round(h2_est * c(0.7, 1, 1.4), 4))
#(p_seq = signif(seq_log(1e-5, 1, length.out = 21), 2))
#(params = expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE)))

#beta_grid = snp_ldpred2_grid(corr1, df_beta, params, ncores=10)

#for (i in 1:ncol(beta_grid)) {
#    res = data.frame(SNP=ss$ID, A1=ss$A1, BETA=beta_grid[,i])
#    write.table(res, file=paste0("sim_",j,"/ldpred_",i,".txt"), row.names=F, col.names=F, quote=F, append=F)
#}

#rm(beta_grid)
#gc()

# LDpred2-auto
#multi_auto <- snp_ldpred2_auto(corr1, df_beta, h2_init=h2_est,
#	vec_p_init = seq_log(1e-4, 0.9, 10), ncores=10)

#beta_auto = sapply(multi_auto, function(auto) auto$beta_est)
#G2 = snp_fastImputeSimple(G)
#pred_auto = big_prodMat(G2, beta_auto)
#sc = apply(pred_auto, 2, sd)
#keep = abs(sc - median(sc)) < 3 * mad(sc)
#final_beta_auto = rowMeans(beta_auto[,keep])

#res = data.frame(SNP=ss$ID, A1=ss$A1, BETA=final_beta_auto)
#write.table(res, file=paste0("sim_",j,"/ldpred_auto.txt"), row.names=F, col.names=F, quote=F, append=F)




