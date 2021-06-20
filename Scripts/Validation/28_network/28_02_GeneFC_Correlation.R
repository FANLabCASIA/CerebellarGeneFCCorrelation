# -------------------------------
# -------------------------------
# Code to Validate Gene-FC correlation analyses in Lobular atlas:
# 
# Uncovering the Genetic Profiles Underlying the Intrinsic Organization of the Human Cerebellum
# Yaping Wang, Lin Chai, Deying Li, Chaohong Gao, Congying Chu, Zhengyi Yang, 
# Yu Zhang, Junhai Xu, Jens Randel Nyengaard1, Bing Liu, Kristoffer Hougaard Madsen, Tianzi Jiang, Lingzhong Fan
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
base.dir <- '/n01dat01/ypwang/AHBA/CerebellarGeneFCCorrelation' 
load(paste(base.dir, '/Data/AHBA/AHBA_original_data/all_data_cere_net28.Rdata', sep = ''))
function.lib <- paste(base.dir, '/Scripts/function_library.R', sep = '')
source(function.lib)
donor.nums <- c('9861', '10021', '12876', '14380', '15496', '15697')
# ----------------

# prepare before
# ----------------
# Cerebellum the code get the saple about the cerebellum should be include in the loop
all_data    <- get_region_info_cere(all_data=all_data, filenames=donor.nums, name='cere')

# Read previously defined info about the overlap of each sample to the cortical atlases
atlas.dir=paste0(base.dir, '/Output/Atlas_overlap')
for ( donor in donor.nums ){ # donor <-'9861'
  atlas.name <- paste('/28net_', donor, '_cere.csv', sep='')
  atlas.in   <- read.csv(file = paste(atlas.dir, atlas.name, sep=''))
  atlas.out  <- atlas.in$x
  cur.n      <- strsplit(atlas.name, '_')[[1]][1]
  use_name   <- paste(gsub('/', '', cur.n), 'cere', sep = '_')
  print(use_name)
  all_data[[donor]][[use_name]] <- atlas.out
}

# Seperately for each donor, average expression of samples in the same cerebellar functional parcel
cere28_names_in  <- read.csv(paste(base.dir, '/Reference_files/28parcellation_names.CSV', sep = ''), header=FALSE)
network_names <- as.character(cere28_names_in$V2)
# ----------------

# Calculate cerebellar differential expression of genes in 6 donors
# ------------------------
all_use <- NULL
for ( donor in donor.nums){ # donor <-'9861'
  split.regions     <- unique(all_data[[donor]][["28net_cere"]])       # split_region index (corresponding to the functional atlas) for each sample
  # use.split.regions <- NULL
  all_use <- c(all_use, unique(split.regions[split.regions != 0]))
}

reg_counts   <- sort(table(all_use))
regions_min  <- reg_counts[reg_counts >= 2]
regions.out  <- as.numeric(rownames(as.matrix(regions_min)))
use.regions.28_mum  <- sort(regions.out) 
use.regions.28      <- network_names[-1][use.regions.28_mum]

out <- averageWithinCereNetworks(all_data=all_data, 
                                 donor.nums=donor.nums, 
                                 use.cere.networks=use.regions.28, # here the  use regions need to be the name, not number 
                                 type='all_cere_micros',
                                 net_names=network_names,
                                 atlas_field='28net_cere')


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
sig.genes       <- NULL
cere.foldchange <- NULL
cere.foldchange[['q05_genes']] <- list()
for ( net in use.regions.28 ) {  
  print(paste('Getting preferential expression for: ', net, sep = ''))
  
  # negative weight for the contrast matrix, depends on number of comparison networks
  mult.term    <- round(1/(length(use.regions.28)-1),6) # round(a,b),四舍五入 , a 是四舍五入的对象， b是保留的小数位
  o.nets       <- use.regions.28[use.regions.28 != net]   # name of the other networks, 
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
cere28.foldchange.n6 <- cere.foldchange
cere28.n6.sig.genes  <- sig.genes
cere28.n6.sig.genes  <- unique(cere28.n6.sig.genes)
name.out     <- paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet25_Gene_geneList_net28.csv', sep = '') # same as generated in the cere script
write.table(x=cere28.n6.sig.genes, file=name.out, quote=TRUE, col.names=FALSE, row.names=FALSE)
# format the differential expression information into one big table
out.matrix  <- NULL
col.headers <- NULL
for ( reg in use.regions.28 ) {
  out.dat     <- cbind(cere28.foldchange.n6$stats[[reg]]$logFC, cere28.foldchange.n6$stats[[reg]]$P.Value, cere28.foldchange.n6$stats[[reg]]$adj.P.Val)
  out.matrix  <- cbind(out.matrix, out.dat)
  col.headers <- c(col.headers, c(paste(reg, 'logFC', sep='_'), paste(reg, 'pval', sep='_'), paste(reg, 'adj.pval', sep='_')))
}
colnames(out.matrix) <- col.headers
rownames(out.matrix) <- rownames(cere28.foldchange.n6$stats[[reg]])
write.csv(x=out.matrix, file=paste0 (base.dir, '/Output/SupplementaryData/SUPPDATA_sheet24_Gene_diffExpr_net10.csv')) 

# Genes that are significantly differentially expressed in one of the cerebellar networks, n=6
write.me    <- NULL
write.genes <- vector()
write.nums  <- NULL
for (buckner in use.regions.28){
  print(buckner)
  genes.tmp   <- as.character(cere28.foldchange.n6[['q05_genes']][[buckner]]$gene_symbol)
  genes.in    <- paste(genes.tmp, collapse = ',')
  write.genes <- rbind.fill.matrix(write.genes, t(genes.tmp)) # t(): change the rows and columns; rbind.fill.matrix: Bind matrices by row, and fill missing columns with NA.
  write.nums  <- c(write.nums, length(genes.tmp))
  write.me    <- rbind(write.me, c(buckner, length(genes.tmp), paste(genes.tmp, collapse = ',')))
}
# Write the significant genes in readable format
write.genes[is.na(write.genes)] <- ''
t.write.genes <- t(write.genes)
write.me      <- rbind(use.regions.28, rbind(write.nums, t.write.genes))
write.csv(x=write.me, file=paste0(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet24_Gene_Network_net10.csv'),row.names=FALSE)
# --------------------------

# calculate the gene-fc scross 28 network
# --------------------------
all_use <- NULL
for ( donor in donor.nums[1:2]){ # donor <-'9861'
  split.regions     <- unique(all_data[[donor]][["28net_cere"]])       # split_region index (corresponding to the functional atlas) for each sample
  # use.split.regions <- NULL
  all_use <- c(all_use, unique(split.regions[split.regions != 0]))
}

reg_counts   <- sort(table(all_use))
regions_min  <- reg_counts[reg_counts >= 2]
regions.out  <- as.numeric(rownames(as.matrix(regions_min)))
use.regions.n2_mum  <- sort(regions.out) 
use.regions.n2      <- network_names[-1][use.regions.n2_mum]

for ( donor in donor.nums[1:2]){
  print(paste('Averaging 28 net cerebellun data for: ', donor, sep = ''))
  # cere.regions.tmp <- unique(all_data[[donor]][["28net_cere"]]) # the region ID 
  # cere.regions     <- sort(cere.regions.tmp[cere.regions.tmp != 0]) # Don't count areas that weren't assigned (i.e not 0)
  
  # Mean Normalized cerebellar expression values, error in the difference between ceres_meanNorm and cere_meanNorm
  all_data[[donor]]$cere28_expr <- averageCereExpr(cere.region=use.regions.n2_mum, 
                                                   data.struct=all_data[[donor]], 
                                                   buckner.names=network_names, 
                                                   samp.labels=all_data[[donor]][["28net_cere"]], 
                                                   data.type='ceres_meanNorm')
  # non-mean Normalized cerebellar expression values
  all_data[[donor]]$cere28_expr_nonorm <- averageCereExpr(cere.region=use.regions.n2_mum, 
                                                          data.struct=all_data[[donor]], 
                                                          buckner.names=network_names, 
                                                          samp.labels=all_data[[donor]][["28net_cere"]], 
                                                          data.type='all_cere_micros')
}


gene.idxs <-  which(rownames(all_data[['10021']]$cere28_expr) %in% unique(cere28.n6.sig.genes)) 
out.mats.28  <- NULL
p_correct <- NULL
p <- NULL
res <- NULL
# Calculate each gene co-experssion matrix seperately for each subject
for ( donor in donor.nums[1:2] ){
  # donor <- donor.nums[[1]]
  cur.mat           <- all_data[[donor]]$cere28_expr[gene.idxs, ] # the difference between cort_expr_17 and cere_expr_17
  untransformed     <- cor(cur.mat, method = 'spearman')
  # untransformed     <- cor(cur.mat, method = 'pearson')
  res[[donor]]      <- rcorr(as.matrix(cur.mat), type = "spearman" )
  # pvalue_adjust <- p.adjust(res2$P, method = "BH")
  # p_correct[[donor]] <- cbind(pvalue_adjust[1:10],pvalue_adjust[11:20],pvalue_adjust[21:30], pvalue_adjust[31:40],
  #                          pvalue_adjust[41:50], pvalue_adjust[51:60], pvalue_adjust[61:70], pvalue_adjust[71:80],
  #                          pvalue_adjust[81:90], pvalue_adjust[91:100])
  # p[[donor]] <- res2$P
  # test <- corr.test(as.matrix(cur.mat), use = "complete", method = 'spearman', adjust = 'BH')
  # rcor.test(as.matrix(cur.mat),p.adjust = FALSE, p.adjust.methods = "BH")
  scaled  <- untransformed # z-tranfsorm spearman correlations
  
  z.corrs <- fisherz(untransformed[upper.tri(untransformed)]) # Upper portion of corr grid
  scaled[upper.tri(untransformed)] <- z.corrs
  
  z.corrs <- fisherz(untransformed[lower.tri(untransformed)])  # Lower portion of corr grid
  scaled[lower.tri(untransformed)] <- z.corrs
  
  # Label the rows/columns
  out.mats.28[[donor]] <- scaled # scaled sig genes expression

}
# Average the correlation matrices
avg.mat   <- (out.mats.28[[1]] + out.mats.28[[2]])/2                            # average the z-transformed correlations of each donors
p.adj <- 0.05/120
test1 <- res$`9861`$P> p.adj
test2 <- res$`10021`$P> p.adj
p.mat <- test1 + test2 


# plot
col2 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                           "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                           "#D6604D", "#B2182B", "#67001F"))
diag(avg.mat) <- 0

pdf(file = paste(base.dir, '/Output/SupplementaryData/Fig/Validation/corrplot_1663genes_net28.pdf', sep = ''),
    width = 10, height = 10)
corrplot(corr = as.matrix(avg.mat), type="full", method = "color",  # diag=FALSE,
         col = col2(100), tl.pos = "lt", tl.cex = 1.2, tl.col = "black",
         cl.lim = c(-0.75,0.75))
dev.off()
diag(avg.mat) <- 0
# Select the rows that correspond to the mRNA data, oj
cere2cere.fcmri.28  <- read.csv(paste0(base.dir,'/Data/FC/gsr_cere28_all.csv'), header = FALSE)
colnames(cere2cere.fcmri.28) <-  network_names[-1]
rownames(cere2cere.fcmri.28) <-  network_names[-1]

use_idxs        <- colnames(cere2cere.fcmri.28) %in% rownames(avg.mat)  # avg.mat is the mRNA coorelation martix of two donors which both have left and right hemisphere 49*49
plot.func       <- cere2cere.fcmri.28[use_idxs, use_idxs]
# plot by corrplot
diag(plot.func) <- 0 
pdf(file = paste(base.dir, "/Output/SupplementaryData/Fig/Validation/corrplot_fc_net28.pdf", sep=''),
    width = 10, height = 10)
corrplot(corr = as.matrix(plot.func), type="full", method = "color",  # diag=FALSE,
         col = col2(100), tl.pos = "lt", tl.cex = 1.2, tl.col = "black",
         cl.lim = c(0,0.6))
dev.off()
diag(plot.func) <- 1
# Global correspondence of mrna/fcmri, here just use the bi-hemisphere donors: 9861, 10021
func.1d.arr <- plot.func[upper.tri(plot.func)]
mrna.1d.arr <- avg.mat[upper.tri(avg.mat)]
test <- cor.test(func.1d.arr, mrna.1d.arr)
write.csv(file = paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet26_Gene_correlation_net28.csv', sep=''), x = avg.mat)
write.csv(file = paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet27__FC_WithinCere_net28.csv', sep=''), x = plot.func)
# --------------------------
