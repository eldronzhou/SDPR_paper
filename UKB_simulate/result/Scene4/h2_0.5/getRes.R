library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)


sz = c("10K", "50K", "100K")
method = c("SDPR","PRS-CS","SBayesR","LDpred","P+T")

res = data.frame(0, nrow=10*3*length(method), ncol=3)
res2 = data.frame(0, nrow=10*3*length(method), ncol=3)
for (j in 1:3) {
        for (i in 1:10) {
                pheno = read.table(paste0("../../../phenotype/Scene1C/h2_0.5/test/sim_",i,".phen"), header=F, stringsAsFactors=F)
                SDPR = read.table(paste0(sz[j],"/predict/SDPR/sim_",i,".profile"), header=T, stringsAsFactors=F)
                stopifnot(all.equal(SDPR$FID, pheno$V1))
                res[50*(j-1)+5*(i-1)+1,1] = sz[j]
                res[50*(j-1)+5*(i-1)+1,2] = method[1]
                res[50*(j-1)+5*(i-1)+1,3] = cor(pheno$V3, SDPR$SCORE)^2
                
                prs_cs = read.table(paste0(sz[j],"/predict/PRS_CS/sim_",i,".profile"), header=T, stringsAsFactors=F)
                stopifnot(all.equal(prs_cs$FID, pheno$V1))
                res[50*(j-1)+5*(i-1)+2,1] = sz[j]
                res[50*(j-1)+5*(i-1)+2,2] = method[2]
                res[50*(j-1)+5*(i-1)+2,3] = cor(pheno$V3, prs_cs$SCORE)^2
                
                gctb = read.table(paste0(sz[j],"/predict/gctb/sim_",i,".profile"), header=T, stringsAsFactors=F)
                stopifnot(all.equal(gctb$FID, pheno$V1))
                res[50*(j-1)+5*(i-1)+3,1] = sz[j]
                res[50*(j-1)+5*(i-1)+3,2] = method[3]
                res[50*(j-1)+5*(i-1)+3,3] = cor(pheno$V3, gctb$SCORE)^2
                
                ldpred = read.table(paste0(sz[j],"/predict/ldpred/sim_",i,".profile"), header=T, stringsAsFactors=F)
                stopifnot(all.equal(gctb$FID, pheno$V1))
                res[50*(j-1)+5*(i-1)+4,1] = sz[j]
                res[50*(j-1)+5*(i-1)+4,2] = method[4]
                res[50*(j-1)+5*(i-1)+4,3] = cor(pheno$V3, ldpred$SCORE)^2
                
                clumping = read.table(paste0(sz[j],"/predict/P+T/sim_",i,".profile"), header=T, stringsAsFactors=F)
                stopifnot(all.equal(gctb$FID, pheno$V1))
                res[50*(j-1)+5*(i-1)+5,1] = sz[j]
                res[50*(j-1)+5*(i-1)+5,2] = method[5]
                res[50*(j-1)+5*(i-1)+5,3] = cor(pheno$V3, clumping$SCORE)^2
        }
}
colnames(res) = c("sample_size","method","R2")
res$sample_size = factor(res$sample_size, levels=c("10K","50K","100K"))
res$method = factor(res$method, levels=method)

tiff("Scenario_4.tiff", units="in", width=6, height=4, res=300)
res = res %>% group_by(sample_size, method) %>% mutate(med=median(R2))
ggplot(data=res, aes(x=sample_size, y=R2,fill=method)) + 
        geom_boxplot(outlier.shape=NA, position=position_dodge(width=.75), width=.4) +
        geom_point(position=position_jitterdodge(dodge.width=.75), color="black", size=.1) +
        ggtitle("Scenario 4") + theme(plot.title=element_text(hjust = 0.5),
                                       text=element_text(size=10))
dev.off()
