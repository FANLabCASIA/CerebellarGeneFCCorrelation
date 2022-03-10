# -------------------------------
# -------------------------------
# Code to run Behavioe-FC-Gene analyses in:
# 
# Uncovering the Genetic Profiles Underlying the Intrinsic Organization of the Human Cerebellum
# Yaping Wang, Lin Chai, Congying Chu, Deying Li, Chaohong Gao, Xia Wu, Zhengyi Yang, Yu Zhang,
# Junhai Xu, Jens Randel Nyengaard, Bing Liu, Kristoffer Hougaard Madsen, Tianzi Jiang, Lingzhong Fan
#
#
# Written by: Yaping Wang
# Contact:    wangyaping19@mails.ucas.ac.cn
# Noted: Part of this code is adapted from Jie Lisa Ji, "Mapping brain-behavior space relationships along the psychosis spectrum"
# -------------------------------
# -------------------------------


## R package needed
# ----------------
# rm(list = ls())
library(ggplot2)
library(scales)
library(matrixStats)
library(plot3D)
library(rgl)
library(ggridges)
# library(reshape2)
# library(fmsb)
# library(gplots)
# library(Hmisc)
# library(grid)
# ----------------


## Function to fo the permutation
# ----------------
sign.pc<-function(x,R=10000,s=59,...){
  # run PCA
  pc.out<-prcomp(x, scale=FALSE,...)
  # the proportion of variance of each PC
  pve=(pc.out$sdev^2/sum(pc.out$sdev^2))[1:s]
  # a matrix with R rows and s columns that contains the proportion of variance explained by each pc for each randomization replicate.
  pve.perm<-matrix(NA,ncol=s,nrow=R)
  for(i in 1:R){
    # permutation each column
    x.perm<-apply(x,2,sample)
    # run PCA
    pc.perm.out<-prcomp(x.perm,scale=FALSE,...)
    # the proportion of variance of each PC.perm
    pve.perm[i,]=(pc.perm.out$sdev^2/sum(pc.perm.out$sdev^2))[1:s]
  }
  # calcalute the p-values
  pval<-apply(t(pve.perm)>pve,1,sum)/R
  return(list(pve=pve,pval=pval, pve.perm=pve.perm))
}
# ----------------


# Read in data
# ----------------
setwd(" ") # Enter your path
base.dir <- ' ' # Enter your path
bh_name <- read.csv(paste(base.dir, '59behavior.csv', sep = ''))
bh_dat <- read.csv(paste(base.dir, 'HCP_behavior.csv', sep = ''))
bh_dat <- bh_dat[bh_dat$X3T_RS.fMRI_Count==4,]
rownames(bh_dat) <- bh_dat[,1]

# Just chose the 218 unrelated subjects to preserve the exchangebility during permutation test
Be_218 <- read.csv(paste(base.dir, '218_unrelated.csv', sep = ''))
Dat <- bh_dat[rownames(bh_dat) %in% Be_218$Subject, which(colnames(bh_dat) %in% bh_name$Formal.Name)]
Dat$Dexterity_AgeAdj <- Dat$Dexterity_AgeAdj -1
Dat_n <- scale(Dat,center=TRUE,scale=TRUE)
Dat_n <- Dat_n[complete.cases(Dat_n),]
write.csv(x = Dat_n, file = 'ToPCA_211.csv')
# ----------------


## Compute PCA and test significance using permutation
# ----------------
pcaA<- prcomp(as.matrix(Dat_n), scale. = FALSE)
write.csv(pcaA$rotation, file = paste0("PCArotations.csv"))
prop_var<- data.frame(cbind((1:59), pcaA$sdev^2/sum(pcaA$sdev^2)))
names(prop_var) <- c("PC", "Variance")

# Permutation 10000 times to obtaine the significance of each PC
signif <- sign.pc(as.matrix(Dat_n),R = 10000,s = 59)
pve.perm <- signif$pve.perm
pve <- signif$pve
pval <- signif$pval

# Write out the PCA score for subjects to palm for fc
n <- length(pval[pval< 0.05])
write.csv(pcaA$rotation[, 1:n], 
          file = " ")

for (j in 1:n) {
  x_cur <- pcaA$x
  write.table( as.numeric(x_cur[,j]), file = paste0("PCAscore_PC", j, ".csv"), 
               row.names = FALSE, col.names = FALSE)
}
# ----------------

