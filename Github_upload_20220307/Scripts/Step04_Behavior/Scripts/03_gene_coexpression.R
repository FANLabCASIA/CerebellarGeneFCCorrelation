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
# -------------------------------
# -------------------------------


## R package need
# ----------------
library(data.table) 
library(WGCNA) 
library(psych) 
library(doBy) 
library(limma) 
library(plyr) 
library(gplots) 
library(ggplot2) 
library(biomartr) 
library(pbapply) 
library(biomaRt) 
library(edgeR) 
library(readxl) 
library(GEOquery) 
library(curl)
library(httr)
library(Hmisc)
library(ltm)
library(pheatmap)
library(corrplot)
library(RColorBrewer)
library(Glimma)
library(ggpubr)
library(ggthemes)
# ----------------


## Modify these filepaths for your local directory structure
# ----------------
setwd("/n02dat01/users/ypwang/AHBA/Github_20220301/Scripts/Step04_Behavior/Results/")
base.dir <- '/n02dat01/users/ypwang/AHBA/Github_20220301' 
load(paste(base.dir, '/Output/RData/all_data_cere_net17.Rdata', sep = '')) # READ AHBA DATA
function.lib <- paste(base.dir, '/Scripts/function_library.R', sep = '')
source(function.lib) # SET BASE DIRECTORY, source function library
donor.nums <- c('9861', '10021', '12876', '14380', '15496', '15697') # SUBJECT LIST
# ----------------


## 0. Prepare before 
# ----------------
# Get info/sample data about all cortical and cerebellar samples
# CORTEX
cort_in     <- read.csv(paste(base.dir, '/Reference_files/cort_regions.csv', sep = ''), header = FALSE) # Ontology IDs corresponding to cortical samples
cortex      <- as.numeric(cort_in$V1)
all_data    <- get_region_info_cort(all_data=all_data, filenames=donor.nums, reg_IDs=cortex, name='cort')

# Cerebellum (the code get the sample about the cerebellum should be include in the loop)
all_data    <- get_region_info_cere(all_data=all_data, filenames=donor.nums, name='cere')

# Read previously defined info about the overlap of each sample to the cortical atlases
all_data <- readAtlasOverlap(all_data, filenames=donor.nums, atlas.dir=paste0(base.dir, '/Output/Atlas_overlap'))

# Read info about split label names for cortex, IDs, and color labels
# this splits the Yeo atlas into 57/114 spatially contiguous regions, depending on whether it's the 7- or 17-network atlas
atlas.key.7  <- read.table(paste0(base.dir, '/Data/Yeo_JNeurophysiol11_SplitLabels/Yeo_JNeurophysiol11_SplitLabels/MNI152/7Networks_ColorLUT_freeview.txt'), col.names=c('ID','Name','R','G','B','A'))
atlas.key.17 <- read.table(paste0(base.dir, '/Data/Yeo_JNeurophysiol11_SplitLabels/Yeo_JNeurophysiol11_SplitLabels/MNI152/17Networks_ColorLUT_freeview.txt'), col.names=c('ID','Name','R','G','B','A'))

# Naming information about the Buckner cerebellar regions (e.g. Default Mode=7)
buckner.names <- NULL
buckner.names[['sev']]     <- as.character(read.csv(paste0(base.dir, '/Reference_files/buckner7_names.CSV'), header=F)$V1)
buckner.names[['sevteen']] <- as.character(read.csv(paste0(base.dir, '/Reference_files/buckner17_names.csv'), header=F)$V1)

# Seperately for each donor, average expression of samples in the same cerebral and cerebellar functional parcel
use.regions.7 <- getCortRegions(all_data, donor.nums, atlas.num='7', thresh=2) 
all_data <-  avg_parcel_expression(all_data, 
                                   cerebellar.atlas='BucknerMNI152_7', # We opt for the 7 network parcellation in cerebellum to minimize functional regions with sparse sampling. 
                                   cerebellar.num='sev',
                                   cortical.atlas='splitLabel_7',
                                   cortical.num=51,
                                   cort.use.regions=use.regions.7)

# Array to match each of the 51 spatially contiguous parcels to a Yeo cortical network name (e.g. Parcel 23=Default)
reg.names.7 <- as.character(buckner.names[['sev']][2:8]) # Left hemispheres
reg.names.17 <- as.character(buckner.names[['sevteen']][2:18]) # both hemispheres
# ----------------


## 1. Link pc-fc to gene
# ----------------
# Construct the gene co-expression matrix for GCI+ and GCI-
plot.order.17 <- read.csv(paste0(base.dir, '/Reference_files/17network_names.csv'), header = FALSE)
plot.order.17 <- plot.order.17[-1,]
twohemi.net17.regions <- getCereRegions(all_data, donor.nums[1:2], '17', 2) # 10 parcel ID 

cere.use.regions.17 <- twohemi.net17.regions 
for ( donor in donor.nums[1:2]){
  print(paste('Averaging cere data for: ', donor, sep = ''))
  
  cerebellar.atlas <- 'BucknerMNI152_17'
  cerebellar.num   <- 17
  cere.use.regions <- twohemi.net17.regions
  all_data[[donor]]$cere_expr_17         <- averageCereExpr(cere.use.regions, all_data[[donor]], buckner.names[['sevteen']], 
                                                            all_data[[donor]][[cerebellar.atlas]], 'cere_meanNorm')
  all_data[[donor]]$cere_expr_nonorm_17  <- averageCereExpr(cere.use.regions, all_data[[donor]], buckner.names[['sevteen']], 
                                                            all_data[[donor]][[cerebellar.atlas]], 'all_cere_micros')
} 

# Read in gene name
GCI_De <- read.csv(file = paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet23_GCI-_List_n197.csv', sep=''))
GCI_De <- GCI_De$X
GCI_In <- read.csv(file = paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet22_GCI+_List_n246.csv', sep=''))
GCI_In <- GCI_In$X

use_net <- buckner.names[['sevteen']][as.integer(cere.use.regions)+1]
gene <- NULL
gene$De <- as.character(GCI_De)
gene$In <- as.character(GCI_In)

gene.idxs <- NULL
gene.idxs$De <- which(rownames(all_data[['10021']]$cere_expr) %in% GCI_De)
gene.idxs$In <- which(rownames(all_data[['10021']]$cere_expr) %in% GCI_In)

# Calculate each gene co-experssion matrix seperately for each subject
out.mats.17  <- NULL
avg.mat <- NULL
for ( i in c("De", "In")) {
  for ( donor in donor.nums[1:2] ){
    cur.mat           <- all_data[[donor]]$cere_expr_17[gene.idxs[[i]], ]  
    untransformed     <- cor(cur.mat, method = 'spearman')
    scaled  <- untransformed 
    
    z.corrs <- fisherz(untransformed[upper.tri(untransformed)])  
    scaled[upper.tri(untransformed)] <- z.corrs
    
    z.corrs <- fisherz(untransformed[lower.tri(untransformed)])   
    scaled[lower.tri(untransformed)] <- z.corrs
    
    out.mats.17[[donor]] <- scaled  
    rownames(out.mats.17[[donor]]) <- plot.order.17[twohemi.net17.regions]
    colnames(out.mats.17[[donor]]) <- plot.order.17[twohemi.net17.regions]
  }
  # Average the correlation matrices
  avg.mat[[i]]   <- (out.mats.17[[1]] + out.mats.17[[2]])/2                          
  write.csv(avg.mat[[i]], file = paste0(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet25_Behavior-FC-Gene_GCI_', i, ".csv"))
}

# Chose which PC is significant
Summ <- matrix(NA, 8, 7)
rownames(Summ) <- paste0("PC", 1:8)
colnames(Summ) <- c("uncp",  "GCI-", "r_-", "p_-", "GCI+", "r_+", "p_+")

fc_pc_p <- NULL
for ( j in   paste0("PC", 1:8)) {
  fc_pc_p[[j]] <- read.csv(paste0(base.dir,'/Scripts/Step04_Behavior/Results/FC_n211_cere_', j, '_dat_ztstat_uncp.csv'),
                           header = FALSE)
  # print(sum(fc_pc_p[[j]]<0.05))
  Summ[j,1]<- sum(fc_pc_p[[j]]<0.05)
}

fc_pc_pm <- NULL
for ( j in   paste0("PC", 1:8)) {
  fc_pc_pm[[j]] <- read.csv(paste0('FC_n211_cere_', j, '_matrix_z_uncp.csv'),
                           header = FALSE)
  colnames(fc_pc_pm[[j]] ) <-  c("VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB", "SalVentAttnA", "SalVentAttnB", "Limbic_Pole", "Limbic_OFC", "ContC", "ContA", "ContB", "TempPar", "DefaultC", "DefaultA", "DefaultB")
  rownames(fc_pc_pm[[j]] ) <-  c("VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB", "SalVentAttnA", "SalVentAttnB", "Limbic_Pole", "Limbic_OFC", "ContC", "ContA", "ContB", "TempPar", "DefaultC", "DefaultA", "DefaultB")
}

# Read in  fc-pc zscore
fc_pc<- NULL
for ( i in paste0("PC", 1:8)) {
  fc_pc[[i]] <- read.csv(paste0('FC_n211_cere_',i,'_matrix_z.csv'),
                         header = FALSE)
  colnames(fc_pc[[i]] ) <-  c("VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB", "SalVentAttnA", "SalVentAttnB", "Limbic_Pole", "Limbic_OFC", "ContC", "ContA", "ContB", "TempPar", "DefaultC", "DefaultA", "DefaultB")
  rownames(fc_pc[[i]] ) <-  c("VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB", "SalVentAttnA", "SalVentAttnB", "Limbic_Pole", "Limbic_OFC", "ContC", "ContA", "ContB", "TempPar", "DefaultC", "DefaultA", "DefaultB")
}

# Calculate the correlation of gene-fc-behaviorPC
for ( i in 1:8) {
  print(i)
  use_idxs        <- colnames(fc_pc[[i]]) %in% rownames(avg.mat$De)  # avg.mat is the mRNA coorelation martix of two donors which both have left and right hemisphere 49*49
  plot.func       <- fc_pc[[i]][use_idxs, use_idxs]
  func.1d.arr     <- plot.func[upper.tri(plot.func)]
  mrna.1d.arr     <- avg.mat$De[upper.tri(avg.mat$De)]
  cor_cur <- cor.test(func.1d.arr, mrna.1d.arr)
  Summ[i, 2] <- cor_cur$p.value < 0.05
  Summ[i, 3] <- cor_cur$estimate
  Summ[i, 4] <- cor_cur$p.value 
}

for ( i in 1:8) {
  print(i)
  use_idxs        <- colnames(fc_pc[[i]]) %in% rownames(avg.mat$In)  # avg.mat is the mRNA coorelation martix of two donors which both have left and right hemisphere 49*49
  plot.func       <- fc_pc[[i]][use_idxs, use_idxs]
  func.1d.arr     <- plot.func[upper.tri(plot.func)]
  mrna.1d.arr     <- avg.mat$In[upper.tri(avg.mat$In)]
  cor_cur <- cor.test(func.1d.arr, mrna.1d.arr)
  Summ[i, 5] <- cor_cur$p.value < 0.05
  Summ[i, 6] <- cor_cur$estimate
  Summ[i, 7] <- cor_cur$p.value 
}
# ----------------


## 2. Permutation for PC3 and GCI-
# ----------------
P =10000
Cor.matrix <- matrix(NA, 3, 1)
rownames(Cor.matrix) <- c("r", "p_unc", "p_perm")
use_idxs        <- colnames(fc_pc[["PC3"]]) %in% rownames(avg.mat[["De"]])  # avg.mat is the mRNA coorelation martix of two donors which both have left and right hemisphere 49*49
plot.func       <- fc_pc[["PC3"]][use_idxs, use_idxs]
func.1d.arr     <- plot.func[upper.tri(plot.func)]
mrna.1d.arr     <- avg.mat[["De"]][upper.tri(avg.mat[["De"]])]
cor_cur <-cor.test(func.1d.arr, mrna.1d.arr)
Cor.matrix[1, ] <- cor_cur$estimate
Cor.matrix[2, ] <- cor_cur$p.value
# permutation test
cor_perm <- NULL
for(i in 1:P){ # i= 4
  # print("permutation=")
  print(i)
  # permutation each column
  gene_perm <- sample(mrna.1d.arr)
  # run cor.test
  cor_perm <-cbind(cor_perm, cor.test(func.1d.arr, gene_perm)$estimate)
}
Cor.matrix[3, ] <- sum(cor_perm>cor_cur$estimate)/P
write.csv(Cor.matrix, file = paste0('Gene_fc_PC3_GCI+.csv'))

# Output the supplementary data
for ( i in c("PC1", "PC3", "PC4")) { 
  p.mat <- fc_pc_pm[[i]]
  r.mat <- fc_pc[[i]]
  write.csv(p.mat, file = paste0(base.dir, "/Output/SupplementaryData/SUPPDATA_sheet25_Behavior-FC-Gene_p_", i, ".csv"))
  write.csv(p.mat, file = paste0(base.dir, "/Output/SupplementaryData/SUPPDATA_sheet25_Behavior-FC-Gene_r_", i, ".csv"))
}
# ----------------


## 3. Plot
# ----------------
col1 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                           "#D1E5F0", "#FFFFFF"))
col2 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                           "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                           "#D6604D", "#B2182B", "#67001F"))
pdf(file = paste0(base.dir, '/Output/Figure/Fig4e_corrplot_GCI-.pdf'),
    width = 10, height = 10)
colnames(avg.mat$De) <- rownames(avg.mat$De)
corrplot(corr = as.matrix(avg.mat$De), type="lower", method = "color",  insig = "label_sig",  diag=FALSE,
         col = col2(100), tl.pos = "ld", tl.cex = 1.2, tl.col = "black", cl.pos = "r")
corrplot(corr = as.matrix(avg.mat$De), add=TRUE, type="lower", method = "color", diag=FALSE, cl.pos = "n",
         tl.pos = "n", tl.cex = 0.9, tl.col = "black", col = col2(100),
         addCoef.col = "black", number.cex = 1)
dev.off()


pdf(file = paste0(base.dir, '/Output/Figure/Fig4b_coefficient_pc3.pdf'),
    width = 10, height = 10)
diag(plot.func) = 1
corrplot(corr = as.matrix(plot.func), is.corr = FALSE, type="lower", method = "color",  insig = "label_sig",  diag=FALSE,
         col = col1(100), tl.pos = "ld", tl.cex = 1.2, tl.col = "black", cl.pos = "r")
corrplot(corr = as.matrix(plot.func), add=TRUE, is.corr = FALSE, type="lower", method = "color", diag=FALSE, cl.pos = "n",
         tl.pos = "n", tl.cex = 0.9, tl.col = "black", col = col1(100),
         addCoef.col = "black", number.cex = 1)
dev.off()
# ----------------