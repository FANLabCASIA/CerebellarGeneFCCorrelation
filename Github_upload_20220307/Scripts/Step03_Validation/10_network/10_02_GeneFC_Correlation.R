# -------------------------------
# -------------------------------
# Code to Validate Gene-FC correlation analyses in MDTB atlas:
# 
# Uncovering the Genetic Profiles Underlying the Intrinsic Organization of the Human Cerebellum
# Yaping Wang, Lin Chai, Congying Chu, Deying Li, Chaohong Gao, Xia Wu, Zhengyi Yang, Yu Zhang,
# Junhai Xu, Jens Randel Nyengaard, Bing Liu, Kristoffer Hougaard Madsen, Tianzi Jiang, Lingzhong Fan
#
#
# Written by: Yaping Wang
# Contact:    wangyaping19@mails.ucas.ac.cn
# Noted: this code is adapted from Kevin M. Anderson, "Gene expression links functional networks across cortex and striatum"
# -------------------------------
# -------------------------------


# R packages need
# -------------------------------
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
library(stringr)
# ----------------


# Modify these filepaths for your local directory structure
# ----------------
base.dir <- '/n02dat01/users/ypwang/AHBA/Github_20220301'
load(paste(base.dir, '/Output/RData/all_data_cere_net10.Rdata', sep = ''))
function.lib <- paste(base.dir, '/Scripts/function_library.R', sep = '')
source(function.lib)
donor.nums <- c('9861', '10021', '12876', '14380', '15496', '15697')
# ----------------


# 0. Prepare before
# ----------------
# Cerebellum the code get the sample about the cerebellum should be include in the loop
all_data    <- get_region_info_cere(all_data=all_data, filenames=donor.nums, name='cere')

# Read in the cerebellar rs-fcMRI correlation matrix
cere2cere.fcmri    <- read.csv(paste0(base.dir,'/Data/FC/gsr_cere10_all.csv'), header = FALSE)
mdtb_names  <- read.csv(paste(base.dir, '/Reference_files/MDTB_10network_names.CSV', sep = ''), header=FALSE)
mdtb_names <- as.character(mdtb_names$V1) # contains None
region.names <- mdtb_names[-1]            # dosen't contain None
colnames(cere2cere.fcmri)   <- region.names
rownames(cere2cere.fcmri)   <- region.names


# Read previously defined info about the overlap of each sample to the cortical atlases
# add the structure: "MDTB_10" to every donor in the all_data_cere.R
atlas.dir=paste0(base.dir, '/Output/Atlas_overlap')
for ( donor in donor.nums ){ # donor <-'9861'
  atlas.name <- paste0('/', 'MDTB_10cere_', donor, '_10net.csv')
  atlas.in   <- read.csv(file = paste(atlas.dir, atlas.name, sep=''))
  atlas.out  <- atlas.in$x
  cur.n      <- strsplit(atlas.name, '_')[[1]][1]
  use_name   <- paste(gsub('/', '', cur.n), '10', sep = '_')
  print(use_name)
  all_data[[donor]][[use_name]] <- atlas.out
}


# Read info about split label names for cortex, IDs, and color labels
atlas.key  <- read.table(paste0(base.dir, '/Reference_files/MDTB_10Regions_Color.txt'), col.names=c('ID','R','G','B'))

# Seperately for each donor, average expression of samples in the same cerebral and cerebellar functional parcel
all_use <- NULL
for ( donor in donor.nums){ # donor <-'9861'
  # split.name        <- paste('BucknerMNI152_', atlas.num, sep = '') # Each individual cortical parcel (values = 1-114), eg., 466 samples for donor 9861
  split.regions     <- unique(all_data[[donor]][["MDTB_10"]])       # split_region index (corresponding to the functional atlas) for each sample
  use.split.regions <- NULL
  all_use <- c(all_use, unique(split.regions[split.regions != 0]))
}
# regions that are present in at least X number of subjects (set with thresh variable)
reg_counts   <- sort(table(all_use))
regions_min  <- reg_counts[reg_counts >= 2]
regions.out  <- as.numeric(rownames(as.matrix(regions_min)))
use.regions.10  <- sort(regions.out) # All 10 regions

# Average the expression of samples in the same cerebral and cerebellar functional parcel
for ( donor in names(all_data)){
  print(paste('Averaging cerebellun data for: ', donor, sep = ''))
  cere.regions.tmp <- unique(all_data[[donor]][["MDTB_10"]]) # the region ID 
  cere.regions     <- sort(cere.regions.tmp[cere.regions.tmp != 0]) # Don't count areas that weren't assigned (i.e not 0)
  
  # Mean Normalized cerebellar expression values, error in the difference between ceres_meanNorm and cere_meanNorm
  all_data[[donor]]$cere_expr <- averageCereExpr(cere.region=cere.regions, 
                                                 data.struct=all_data[[donor]], 
                                                 buckner.names=mdtb_names, 
                                                 samp.labels=all_data[[donor]][["MDTB_10"]], 
                                                 data.type='ceres_meanNorm')
  # non-mean Normalized cerebellar expression values
  all_data[[donor]]$cere_expr_nonorm <- averageCereExpr(cere.region=cere.regions, 
                                                        data.struct=all_data[[donor]], 
                                                        buckner.names=mdtb_names, 
                                                        samp.labels=all_data[[donor]][["MDTB_10"]], 
                                                        data.type='all_cere_micros')
}
# ----------------


## 1. Segregation: Network-specific genes
# ----------------
# Calculate cerebellar differential expression of genes in all 6 donors 
# First get the use regions which satisfy the level
all_use <- NULL
for ( donor in donor.nums ){ # donor <-'9861'
  # split.name        <- paste('BucknerMNI152_', atlas.num, sep = '') # Each individual cortical parcel (values = 1-114), eg., 466 samples for donor 9861
  split.regions     <- unique(all_data[[donor]][["MDTB_10"]])       # split_region index (corresponding to the functional atlas) for each sample
  use.split.regions <- NULL
  all_use <- c(all_use, unique(split.regions[split.regions != 0]))
}
# regions that are present in at least X number of subjects (set with thresh variable)
reg_counts   <- sort(table(all_use))
regions_min  <- reg_counts[reg_counts >= 2]
regions.out  <- as.numeric(rownames(as.matrix(regions_min)))
cere.use.regions_num  <- sort(regions.out) # All 10 regions
cere.use.regions  <- mdtb_names[-1][cere.use.regions_num]


out <- averageWithinCereNetworks(all_data=all_data, 
                                 donor.nums=donor.nums, 
                                 use.cere.networks=cere.use.regions, 
                                 type='all_cere_micros',
                                 net_names=mdtb_names,
                                 atlas_field='MDTB_10')


expr    <- out[[1]]  #  all.expr = num [1:20738, 1:28], rows means the 20738 genes, colums means the 4 donor * 7 networks name
regions <- out[[2]]  #  region_arr = donor_arr = chr[1: 28], 28 = 7 network * 4 donor, 28 is the region ID, 4 donor * 7 networks name
donors  <- out[[3]]  #  So here we get the averaged expression value for each gene(20738) in each networks(7) in each donors(4)
colnames(expr) <- paste(donors, regions, sep = '_')

# Use limma to calculate differential expression for each network, relative to all others
fac              <- as.factor(regions)                # the network type of each column in 'expr', get 7 levels which corresponde to 7 networks name
design           <- model.matrix(~0 + fac)            # design matrix, 0 means no intercept
colnames(design) <- gsub('fac', '', colnames(design)) # colnames(design) = fac+networks name, eg., facCont, after this step, we get the colname euqal to the networks name
corfit           <- duplicateCorrelation(expr, design, block=donors) 

# Calculate Differential expression for each region
region.names <- mdtb_names[-1]
sig.genes       <- NULL
cere.foldchange <- NULL
cere.foldchange[['q05_genes']] <- list()
for ( net in cere.use.regions ) {  
  print(paste('Getting preferential expression for: ', net, sep = ''))
  
  # negative weight for the contrast matrix, depends on number of comparison networks
  mult.term    <- round(1/(length(cere.use.regions)-1),6) # round(a,b),四舍五入 , a 是四舍五入的对象， b是保留的小数位
  o.nets       <- cere.use.regions[cere.use.regions != net]   # name of the other networks, 
  cur.contrast <- paste('1*', net, '-', mult.term, '*', paste(o.nets, collapse = paste('-', mult.term, '*', sep = '')), sep = '')
  
  # Make the contrast matrix
  cmtx <- NULL
  cmtx <- makeContrasts(contrasts=cur.contrast, levels=colnames(design)) # cmtx is our contrast matrix, eg., for vis, is1*Vis-0.166667*Default-0.166667*Cont-0.166667*Limbic-0.166667*VentAttn-0.166667*DorsAttn-0.166667*SomMot
  tmplm   <- lmFit(expr, design, block=donors, correlation=corfit$consensus.correlation ) # Fit the linear model to the data
  fit     <- eBayes(contrasts.fit( tmplm, cmtx ) ) # contrast.fit: fit the linear model to estimate a set of contrast
  cere.foldchange[['fit_df']][[net]] <- fit        # cere.foldchange has three list: fit_df = fit; stats= ordered topTable(fit, number=Inf); q01_genes = row names of selected stats
  tmp     <- topTable(fit, number=Inf)             # A number of summary statistics are presented by topTable() for the top genes and the selected contrast.
  cere.foldchange[['stats']][[net]] <- tmp[order(rownames(tmp)),] 
  pos.idxs       <- which(cere.foldchange$stats[[net]]$logFC > 0)         # The logFC column gives the value of the contrast. Usually this represents a log2-fold change  between two or more experimental 
  adjusted.ps    <- which(cere.foldchange$stats[[net]]$adj.P.Val <= .05)  # adj.P.Value is the p-value adjusted for multiple testing
  inter     <- intersect(adjusted.ps, pos.idxs)
  genes     <- rownames(cere.foldchange$stats[[net]])[inter ]
  print(length(genes)) 
  sig.genes <- c(sig.genes, genes)
  logFC    <- cere.foldchange[['stats']][[net]]$logFC[inter]
  P.Value   <- cere.foldchange[['stats']][[net]]$P.Value[inter ] 
  adj.P.Val <-  cere.foldchange$stats[[net]]$adj.P.Val[inter ]
  IDmatrix <- all_data[["9861"]][["probes_collapse"]]
  genes <- IDmatrix[(IDmatrix[,"gene_symbol"]%in%genes), ]
  cere.foldchange[['q05_genes']][[net]] <- cbind(genes,  logFC, P.Value, adj.P.Val)
}
cere.foldchange.n6 <- cere.foldchange
cere.n6.sig.genes  <- sig.genes
cere.n6.sig.genes <- unique(cere.n6.sig.genes) # 481

# Overlap with main strategy
sig.All <- read.csv(paste(base.dir,"/Output/SupplementaryData/SUPPDATA_sheet3_Gene_cere_GeneList.csv", sep = ""), header = FALSE)
sig.overlap <- intersect(cere.n6.sig.genes, as.matrix(sig.All)) 
write.csv(file = paste(base.dir, 
                       '/Output/SupplementaryData/SUPPDATA_sheet4_Gene_cere_Specificity_Overlap_net10.csv', 
                       sep=''), x = sig.overlap )

name.out     <- paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet4_Gene_cere_Specificity_geneList_net10.csv', sep = '') # same as generated in the cere script
write.table(x=cere.n6.sig.genes, file=name.out, quote=TRUE, col.names=FALSE, row.names=FALSE)

# SUPPLEMENTAL DATA TABLE
# format the differential expression information into one big table
out.matrix  <- NULL
col.headers <- NULL
for ( reg in cere.use.regions ) {
  out.dat     <- cbind(cere.foldchange.n6$stats[[reg]]$logFC, cere.foldchange.n6$stats[[reg]]$P.Value, cere.foldchange.n6$stats[[reg]]$adj.P.Val)
  out.matrix  <- cbind(out.matrix, out.dat)
  col.headers <- c(col.headers, c(paste(reg, 'logFC', sep='_'), paste(reg, 'pval', sep='_'), paste(reg, 'adj.pval', sep='_')))
}
colnames(out.matrix) <- col.headers
rownames(out.matrix) <- rownames(cere.foldchange.n6$stats[[reg]])
write.csv(x=out.matrix, file=paste0(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet4_Gene_cere_Specificity_diffExpr_net10.csv')) # same as generated in the cere script

# Genes that are significantly differentially expressed in one of the cerebellar networks, n=6
write.me    <- NULL
write.genes <- vector()
write.nums  <- NULL
for (buckner in cere.use.regions){
  print(buckner)
  genes.tmp   <- as.character(cere.foldchange.n6[['q05_genes']][[buckner]]$gene_symbol)
  genes.in    <- paste(genes.tmp, collapse = ',')
  write.genes <- rbind.fill.matrix(write.genes, t(genes.tmp)) # t(): change the rows and columns; rbind.fill.matrix: Bind matrices by row, and fill missing columns with NA.
  write.nums  <- c(write.nums, length(genes.tmp))
  write.me    <- rbind(write.me, c(buckner, length(genes.tmp), paste(genes.tmp, collapse = ',')))
}
# Write the significant genes in readable format
write.genes[is.na(write.genes)] <- ''
t.write.genes <- t(write.genes)
write.me      <- rbind(cere.use.regions, rbind(write.nums, t.write.genes))
write.csv(x=write.me, file=paste0(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet4_Gene_cere_Specificity_Network_net10.csv'),row.names=FALSE)
# ----------------


## 2. Integration: Gene-FC correlation in MDTB parcellation
# ----------------
# Make the mRNA cerebellum-cerebellum correlation matrix
plot.order <- region.names

# caclulate the region will be used
all_use <- NULL
for ( donor in donor.nums[1:2] ){ # donor <-'9861'
  split.regions     <- unique(all_data[[donor]][["MDTB_10"]])       # split_region index (corresponding to the functional atlas) for each sample
  use.split.regions <- NULL
  all_use <- c(all_use, unique(split.regions[split.regions != 0]))
}
# regions that are present in at least X number of subjects (set with thresh variable)
reg_counts   <- sort(table(all_use))
regions_min  <- reg_counts[reg_counts >= 2]
regions.out  <- as.numeric(rownames(as.matrix(regions_min)))
twohemi.net10.regions  <- sort(regions.out) #
# twohemi.net10.regions  <- mdtb_names[-1][twohemi.net10.regions]

# Average expression within each parcel that contains data from each bi-hemispheric subject. 
cere.use.regions.10 <- twohemi.net10.regions # 10 parcel ID from both bi-hemispheric AHBA donors, same as the twohemi.net17.regions
for ( donor in donor.nums[1:2]){
  print(paste('Averaging cere data for: ', donor, sep = ''))
  
  # Set the modifiable atlas parameters, for instance you can run cerebellum 17/cortical 17
  cerebellar.atlas <- 'MDTB_10'
  # cerebellar.num   <- 17
  cere.use.regions <- twohemi.net10.regions
  # Mean Normalized cerebellar expression values
  all_data[[donor]]$cere_expr_10         <- averageCereExpr(cere.use.regions, all_data[[donor]], mdtb_names, 
                                                            all_data[[donor]][["MDTB_10"]], 'ceres_meanNorm')
  
  # non-Mean Normalized cerebellar expression values
  all_data[[donor]]$cere_expr_nonorm_10  <- averageCereExpr(cere.use.regions, all_data[[donor]], mdtb_names, 
                                                            all_data[[donor]][["MDTB_10"]], 'all_cere_micros')
} # here for the cere_expr_17, we just get 10 column, which correspond to the 10 region ID derived from the both bi-hemisphere donor
# Network-associated genes used for the cortico-cortical correlations

gene.idxs <- which(rownames(all_data[['10021']]$cere_expr_10) %in% unique(cere.n6.sig.genes))
out.mats.10  <- NULL
res <- NULL
p_correct <- NULL
p <- NULL
# Calculate each gene co-experssion matrix seperately for each subject
for ( donor in donor.nums[1:2] ){
  # donor <- donor.nums[[1]]
  cur.mat           <- all_data[[donor]]$cere_expr_10[gene.idxs, ] # the difference between cort_expr_17 and cere_expr_17
  untransformed     <- cor(cur.mat, method = 'spearman')
  res[[donor]]      <- rcorr(as.matrix(cur.mat), type = "spearman" )
  
  scaled  <- untransformed # z-tranfsorm spearman correlations
  
  z.corrs <- fisherz(untransformed[upper.tri(untransformed)]) # Upper portion of corr grid
  scaled[upper.tri(untransformed)] <- z.corrs
  
  z.corrs <- fisherz(untransformed[lower.tri(untransformed)])  # Lower portion of corr grid
  scaled[lower.tri(untransformed)] <- z.corrs
  
  # Label the rows/columns
  out.mats.10[[donor]] <- scaled # scaled sig genes expression
  rownames(out.mats.10[[donor]]) <- plot.order[twohemi.net10.regions]
  colnames(out.mats.10[[donor]]) <- plot.order[twohemi.net10.regions]
}
# Average the correlation matrices
avg.mat   <- (out.mats.10[[1]] + out.mats.10[[2]])/2                            # average the z-transformed correlations of each donors
rownames(avg.mat) <- plot.order[twohemi.net10.regions]
colnames(avg.mat) <- plot.order[twohemi.net10.regions]
p.adj <- 0.05/length(avg.mat[upper.tri(avg.mat)])
test1 <- res$`9861`$P> p.adj
test2 <- res$`10021`$P> p.adj
p.mat <- test1 + test2 

# plot
col2 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                           "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                           "#D6604D", "#B2182B", "#67001F"))

diag(avg.mat) <- 0
pdf(file = paste(base.dir, 'Output/SupplementaryData/SupplementaryFigure/FigS2a_corrplot_481genes_net10.pdf', sep = ''),
    width = 10, height = 10)
rownames(avg.mat)
# corrplot(corr = as.matrix(avg.mat), type="full", method = "color",  # diag=FALSE,
#          col = col2(100), tl.pos = "lt", tl.cex = 1.2, tl.col = "black",
#          col.lim = c(-0.75,0.75))
corrplot(corr = as.matrix(avg.mat), type="upper", method = "color",  # diag=FALSE,
         col = col2(100), tl.pos = "lt", tl.cex = 1.2, tl.col = "black",
         col.lim = c(-0.75, 0.75))
corrplot(corr = as.matrix(avg.mat), add=TRUE, type="upper", method = "color", diag=FALSE, cl.pos = "n",
         tl.pos = "n", tl.cex = 0.9, tl.col = "black", col = col2(100),
         addCoef.col = "black",  number.digits = 2, number.cex = 1.0)
corrplot(corr = as.matrix(avg.mat), add=TRUE, type="lower", method="color",diag=FALSE,tl.pos="n", cl.pos="n", 
         col = col2(100), 
         p.mat = as.matrix(p.mat), insig = "label_sig", sig.level = .05, pch.cex = 1.0, pch.col = "black")
dev.off()
diag(avg.mat) <- 1

use_idxs        <- colnames(cere2cere.fcmri) %in% rownames(avg.mat)  # avg.mat is the mRNA coorelation martix of two donors which both have left and right hemisphere 49*49
plot.func       <- cere2cere.fcmri[use_idxs, use_idxs]
diag(plot.func) <- 0 

pdf(file = paste(base.dir, "/Output/SupplementaryData/SupplementaryFigure/FigS2b_corrplot_fc_net10.pdf", sep=''),
    width = 10, height = 10)
# corrplot(corr = as.matrix(plot.func), type="full", method = "color",  # diag=FALSE,
#          col = col2(100), tl.pos = "lt", tl.cex = 1.2, tl.col = "black",
#          col.lim = c(0,0.6))
corrplot(corr = as.matrix(plot.func), type="full", method = "color",  # diag=FALSE,
         col = col2(100), tl.pos = "lt", tl.cex = 1.2, tl.col = "black",
         col.lim = c(0,0.7))
corrplot(corr = as.matrix(plot.func), add=TRUE, type="upper", method = "color", diag=FALSE, cl.pos = "n",
         tl.pos = "n", tl.cex = 0.9, tl.col = "black", col = col2(100), col.lim = c(0,0.7),
         addCoef.col = "black",  number.digits = 2, number.cex = 1.0)
dev.off()
diag(plot.func) <- 1

func.1d.arr   <- plot.func[upper.tri(plot.func)]
mrna.1d.arr <- avg.mat[upper.tri(avg.mat)]
cor.test(func.1d.arr, mrna.1d.arr)    
write.csv(file = paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet15_Gene_cere_Corr_net10.csv', sep=''), x = avg.mat)
write.csv(file = paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet16_FC_cere_Corr_net10.csv', sep=''), x = plot.func)
# ----------------