library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)

bd = read.table("../reportdate.tab", sep="\t", header=T, stringsAsFactors=F)

lvl.1001 <- c(-3,-1,1,2,3,4,8,6,1001,1002,1003,2001,2002,2003,2004,3001,3002,3003,3004,4001,4002,4003)
lbl.1001 <- c("Prefer not to answer","Do not know","White","Mixed","Asian or Asian British","Black or Black British","Chinese","Other ethnic group","British","Irish","Any other white background","White and Black Caribbean","White and Black African","White and Asian","Any other mixed background","Indian","Pakistani","Bangladeshi","Any other Asian background","Caribbean","African","Any other Black background")
bd$f.21000.0.0 <- ordered(bd$f.21000.0.0, levels=lvl.1001, labels=lbl.1001)
bd$f.21000.1.0 <- ordered(bd$f.21000.1.0, levels=lvl.1001, labels=lbl.1001)
bd$f.21000.2.0 <- ordered(bd$f.21000.2.0, levels=lvl.1001, labels=lbl.1001)
lvl.0009 <- c(0,1)
lbl.0009 <- c("Female","Male")
bd$f.22001.0.0 <- ordered(bd$f.22001.0.0, levels=lvl.0009, labels=lbl.0009)

pcs = read.table("../pcs.txt", header=T)

SDPR = read.table("SDPR/res_UKB.profile", header=T, stringsAsFactors=F)

gctb = read.table("gctb/gctb.profile", header=T, stringsAsFactors=F)

read_profile = function(dir) {
        ldpred_files = list.files(path=dir, pattern="*.profile.gz$", full.names=T)
        for (i in 1:length(ldpred_files)) {
                if (i == 1) {
                        ldpred = read.table(ldpred_files[i], header=T)
                        ldpred = ldpred[,c(1,6)]
                        colnames(ldpred) = c("FID", ldpred_files[i])
                }
                else {
                        tmp = read.table(ldpred_files[i], header=T)
                        tmp = tmp[,c(1,6)]
                        colnames(tmp) = c("FID", ldpred_files[i])
                        ldpred = inner_join(ldpred, tmp, by="FID") 
                }
        }
        ldpred
}

PRS_cs = read_profile("./PRS_CS/")

ldpred = read_profile("./ldpred/")

ldpred2 = read.table("ldpred2/ldpred2.sscore.gz", 
                     header=F, stringsAsFactors=F)[,-c(2:8)]
colnames(ldpred2)[1] = "FID"

lassosum = read.table("lassosum/lassosum.sscore.gz",
                      header=F, stringsAsFactors=F)[,-c(2:8)]
colnames(lassosum)[1] = "FID"

dbslmm = read_profile("dbslmm/")

clumping = read_profile("./P+T/")

height = inner_join(bd, SDPR, by=c("f.eid"="FID")) %>% 
        inner_join(pcs, by=c("f.eid"="FID"))

rm_idx = is.na(height$f.50.0.0) | is.na(height$f.22001.0.0) | is.na(height$f.21003.0.0)

height = height[!(rm_idx), c("f.eid","f.50.0.0","f.22001.0.0","f.21003.0.0",paste0("PC",1:10))]
colnames(height) = c("ID","height","sex","age", paste0("PC",1:10))
lm_fit_base = lm(height~sex+age, data=height)
height$adjust_height = residuals(lm_fit_base)
fid = (height$ID)

get_cor = function(dat, i=-1) {
        if (i == -1)
                a = summary(lm(height~sex+age+SCORE+PC1+PC2+PC3+PC4+PC8+PC6+PC7+PC8+PC9+PC10, data=dat))$residuals
        else
                a = summary(lm(height~sex+age+dat[,i]+PC1+PC2+PC3+PC4+PC8+PC6+PC7+PC8+PC9+PC10, data=dat))$residuals
        b = summary(lm(height~sex+age+PC1+PC2+PC3+PC4+PC8+PC6+PC7+PC8+PC9+PC10, data=dat))$residuals
        1 - sum(a^2)/sum(b^2)
}

set.seed(1)
res = replicate(10, {
        idx = sort(sample.int(length(fid), size=round(length(fid)/2), replace=F))
        validate_idx = fid[idx]
        test_idx = fid[-idx]
        
        # ldpred
        validate_dat = height[height$ID %in% validate_idx, ] %>%
                left_join(ldpred, by=c("ID"="FID"))
        
        validate_cor = rep(0, ncol(validate_dat)-16)
        for (i in 16:ncol(validate_dat)) {
                validate_cor[i] = get_cor(validate_dat, i)
        }
        max_idx = which.max(validate_cor)
        print(paste0(colnames(validate_dat)[max_idx], ": ",
                     validate_cor[max_idx]))
        test_dat = height[height$ID %in% test_idx, ] %>%
                left_join(ldpred, by=c("ID"="FID"))
        
        cor_ldpred = get_cor(test_dat, max_idx)
        
        # ldpred2
        validate_dat = height[height$ID %in% validate_idx, ] %>%
                left_join(ldpred2, by=c("ID"="FID"))

        validate_cor = rep(0, ncol(validate_dat)-16)
        for (i in 16:ncol(validate_dat)) {
                validate_cor[i] = get_cor(validate_dat, i)
        }
        max_idx = which.max(validate_cor)
        print(paste0(colnames(validate_dat)[max_idx], ": ",
                     validate_cor[max_idx]))
        test_dat = height[height$ID %in% test_idx, ] %>%
                left_join(ldpred2, by=c("ID"="FID"))

        cor_ldpred2 = get_cor(test_dat, max_idx)
        
        # dbslmm
        validate_dat = height[height$ID %in% validate_idx, ] %>%
                left_join(dbslmm, by=c("ID"="FID"))
        
        validate_cor = rep(0, ncol(validate_dat)-16)
        for (i in 16:ncol(validate_dat)) {
                validate_cor[i] = get_cor(validate_dat, i)
        }
        max_idx = which.max(validate_cor)
        print(paste0(colnames(validate_dat)[max_idx], ": ",
                     validate_cor[max_idx]))
        test_dat = height[height$ID %in% test_idx, ] %>%
                left_join(dbslmm, by=c("ID"="FID"))
        
        cor_dbslmm = get_cor(test_dat, max_idx)
        
        # lassosum
        validate_dat = height[height$ID %in% validate_idx, ] %>%
                left_join(lassosum, by=c("ID"="FID"))
        
        validate_cor = rep(0, ncol(validate_dat)-16)
        for (i in 16:ncol(validate_dat)) {
                validate_cor[i] = get_cor(validate_dat, i)
        }
        max_idx = which.max(validate_cor)
        print(paste0(colnames(validate_dat)[max_idx], ": ",
                     validate_cor[max_idx]))
        test_dat = height[height$ID %in% test_idx, ] %>%
                left_join(lassosum, by=c("ID"="FID"))
        
        cor_lassosum = get_cor(test_dat, max_idx)
        
        validate_dat = height[height$ID %in% validate_idx, ] %>%
                left_join(clumping, by=c("ID"="FID"))
        validate_cor = rep(0, ncol(validate_dat)-16)
        for (i in 16:ncol(validate_dat)) {
                validate_cor[i] = get_cor(validate_dat, i)
        }
        max_idx = which.max(validate_cor)
        print(paste0(colnames(validate_dat)[max_idx], ": ",
                     validate_cor[max_idx]))
        test_dat= height[height$ID %in% test_idx, ] %>%
                left_join(clumping, by=c("ID"="FID"))
        
        cor_clumping = get_cor(test_dat, max_idx)
        
        test_dat = height[height$ID %in% test_idx, ] %>%
                left_join(SDPR, by=c("ID"="FID"))
        
        cor_SDPR = get_cor(test_dat)
        
        test_dat = height[height$ID %in% test_idx, ] %>%
                left_join(gctb, by=c("ID"="FID"))
        
        cor_gctb = get_cor(test_dat)
        
        test_dat = height[height$ID %in% test_idx, ] %>%
                left_join(PRS_cs, by=c("ID"="FID"))
        
        validate_dat = height[height$ID %in% validate_idx, ] %>%
                left_join(PRS_cs, by=c("ID"="FID"))
        
        validate_cor = rep(0, ncol(validate_dat)-16)
        for (i in 16:ncol(validate_dat)) {
                validate_cor[i] = get_cor(validate_dat, i)
        }
        max_idx = which.max(validate_cor)
        print(paste0(colnames(validate_dat)[max_idx], ": ",
                     validate_cor[max_idx]))
        test_dat = height[height$ID %in% test_idx, ] %>%
                left_join(PRS_cs, by=c("ID"="FID"))
        
        cor_PRSCS = get_cor(test_dat, max_idx)
        c(cor_SDPR, cor_PRSCS, cor_gctb, cor_ldpred, cor_clumping,
          cor_ldpred2, cor_lassosum, cor_dbslmm)
        print(c(cor_SDPR, cor_PRSCS, cor_gctb, cor_ldpred, cor_clumping,
                cor_ldpred2, cor_lassosum, cor_dbslmm))
})

res = t(res)
res2 = data.frame(nrow=10*8, ncol=2)
for (i in 1:10) {
        res2[8*(i-1)+1,1] = res[i,1]; res2[8*(i-1)+1,2] = "SDPR"
        res2[8*(i-1)+2,1] = res[i,2]; res2[8*(i-1)+2,2] = "PRS-CS"
        res2[8*(i-1)+3,1] = res[i,3]; res2[8*(i-1)+3,2] = "SBayesR"
        res2[8*(i-1)+4,1] = res[i,4]; res2[8*(i-1)+4,2] = "LDpred"
        res2[8*(i-1)+5,1] = res[i,5]; res2[8*(i-1)+5,2] = "P+T"
        res2[8*(i-1)+6,1] = res[i,6]; res2[8*(i-1)+6,2] = "LDpred2"
        res2[8*(i-1)+7,1] = res[i,7]; res2[8*(i-1)+7,2] = "Lassosum"
        res2[8*(i-1)+8,1] = res[i,8]; res2[8*(i-1)+8,2] = "DBSLMM"
}

colnames(res2) = c("R2","method")
res2$method = factor(res2$method, levels=c("SDPR","PRS-CS","SBayesR","LDpred","P+T",
                                           "LDpred2","Lassosum","DBSLMM"))
res2 = res2 %>% group_by(method) %>% mutate(med=mean(R2))

res3 = data.frame(mean=colMeans(res), sd=apply(res, 2, sd))
method = c("SDPR","PRS-CS","SBayesR","LDpred","P+T",
           "LDpred2","Lassosum","DBSLMM")
res3$method = factor(method, levels=method)

tiff("height.tiff", units="in", width=6, height=4, res=300)
ggplot(res3, aes(x=method, y=mean, fill=method)) + 
        geom_bar(stat="identity") + 
        geom_errorbar(aes(x=method, ymin=mean-sd, ymax=mean+sd), width=.5) + 
        ggtitle("Height") + ylim(c(0,.3)) + ylab("R2") + theme_bw(10) +
        theme(plot.title=element_text(hjust=0.5), text=element_text(size=10))  
dev.off()

set.seed(1)
res_auto = replicate(10, {
        idx = sort(sample.int(length(fid), size=round(length(fid)/2), replace=F))
        validate_idx = fid[idx]
        test_idx = fid[-idx]
        
        test_dat = height[height$ID %in% test_idx, ] %>%
                left_join(SDPR, by=c("ID"="FID"))
        
        cor_SDPR = get_cor(test_dat)
        
        test_dat = height[height$ID %in% test_idx, ] %>%
                left_join(PRS_cs[,c(1,ncol(PRS_cs))], by=c("ID"="FID"))
        colnames(test_dat)[ncol(test_dat)] = "SCORE"
        cor_PRSCS = get_cor(test_dat)
        
        test_dat = height[height$ID %in% test_idx, ] %>%
                left_join(gctb, by=c("ID"="FID"))
        cor_gctb = get_cor(test_dat)
        
        test_dat = height[height$ID %in% test_idx, ] %>%
                left_join(ldpred2[,c(1,ncol(ldpred2))], by=c("ID"="FID"))
        colnames(test_dat)[ncol(test_dat)] = "SCORE"
        cor_ldpred2 = get_cor(test_dat)
        
        c(cor_SDPR, cor_PRSCS, cor_gctb, cor_ldpred2)
        print(c(cor_SDPR, cor_PRSCS, cor_gctb, cor_ldpred2))
})
colMeans(t(res_auto))


