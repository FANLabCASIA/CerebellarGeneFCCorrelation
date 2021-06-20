# -------------------------------
# -------------------------------
# Code to run Permutation test in:
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

# R package need
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
# -------------------------------

# Modify these filepaths for your local directory structure
# -------------------------------
base.dir <- '/n01dat01/ypwang/AHBA/CerebellarGeneFCCorrelation' 

# READ AHBA DATA
load(paste(base.dir, '/Data/AHBA/AHBA_original_data/all_data_cere_net17.Rdata', sep = ''))

# SET BASE DIRECTORY, source function library
function.lib <- paste(base.dir, '/Scripts/function_library.R', sep = '')
source(function.lib)

# SUBJECT LIST
donor.nums <- c('9861', '10021', '12876', '14380', '15496', '15697')
# -------------------------------

## Prepare before
# ----------------
# Naming information about the Buckner cerebellar regions (e.g. Default Mode=7)
buckner.names <- NULL
buckner.names[['sev']]     <- as.character(read.csv(paste0(base.dir, '/Reference_files/buckner7_names.CSV'), header=F)$V1)
buckner.names[['sevteen']] <- as.character(read.csv(paste0(base.dir, '/Reference_files/buckner17_names.csv'), header=F)$V1)

reg.names.7 <- as.character(buckner.names[['sev']][2:8]) # Left hemispheres
reg.names.17 <- as.character(buckner.names[['sevteen']][2:18]) # both hemispheres

# genes/regions in cortical analysis
gene.list    <- rownames(all_data[[1]]$all_cere_micros)
region.names <- c('Default', 'Cont', 'Limbic', 'VentAttn', 'DorsAttn', 'SomMot', 'Vis')

# Cerebellum the code get the saple about the cerebellum should be include in the loop
all_data    <- get_region_info_cere(all_data=all_data, filenames=donor.nums, name='cere')
atlas_dir <- paste(base.dir, '/Output/Atlas_overlap/', sep = '')
for ( donor in donor.nums ){ # donor <- "9861"
  types  <- 'BucknerMNI152_cere_'
  atlas.name <- paste0('/', types, donor, '_', '17', 'net.csv')
  atlas_in   <- read.csv(file = paste(atlas_dir, atlas.name, sep=''))
  cur_n    <- strsplit(atlas.name, '_')[[1]][1]
  use_name <- paste(gsub('/', '', cur_n), '17', sep = '_')
  all_data[[donor]][[use_name]] <- atlas_in$x 
}


plot.order.17 <- read.csv(paste0(base.dir, '/Reference_files/17network_names.csv'), header = FALSE)
plot.order.17 <- plot.order.17[-1,]
twohemi.net17.regions <- getCereRegions(all_data, donor.nums[1:2], '17', 2) # 10 parcel ID 

# Average expression within each parcel that contains data from each bi-hemispheric subject. 
cere.use.regions.17 <- twohemi.net17.regions 
# -------------------------------

# Permutation test
# -------------------------------
Sig <- NULL
R_random <- NULL
for (i in 1:10000){ # i=1
  for ( donor in donor.nums ){ # donor <- "9861"
    types  <- 'BucknerMNI152_cere_'
    atlas.name <- paste0('/', types, donor, '_', '7', 'net.csv')
    atlas_in   <- read.csv(file = paste(atlas_dir, atlas.name, sep=''))
    cur_n    <- strsplit(atlas.name, '_')[[1]][1]
    use_name <- paste(gsub('/', '', cur_n), '7', sep = '_')
    all_data[[donor]][[use_name]] <- sample(atlas_in$x ,length(atlas_in$x))
    # print(use_name)
  }
  use.regions.7 <- getCortRegions(all_data, donor.nums, atlas.num='7', thresh=2) 
  all_data <- avg_parcel_expression(all_data, 
                                     cerebellar.atlas='BucknerMNI152_7', # We opt for the 7 network parcellation in cerebellum to minimize functional regions with sparse sampling. 
                                     cerebellar.num='sev',
                                     cortical.atlas='splitLabel_7',
                                     cortical.num=51,
                                     cort.use.regions=use.regions.7)
  
  ## Calculate the network-specific genes across 7 network
  cere.lh.use.regions.7 <- getCereRegions(all_data, donor.nums[3:6], atlas.num='7', thresh=2)  # 7 parcel ID (each parcel have some samples) of left hemisphere donor used to define the differential gene 
  cere.lh.use.regions.7 <- c('Default', 'Cont', 'Limbic', 'VentAttn', 'DorsAttn', 'SomMot', 'Vis') # because for all cere networks, they have more than 2 donors samples, which get from before line
  region.names         <- c('Vis','VentAttn','DorsAttn','SomMot','Limbic','Cont','Default')
  use.n6.net7.regions  <- getCereRegions(all_data, donor.nums, '7', 2) # both have , so same as done in 4 donors
  out <- averageWithinCereNetworks(all_data=all_data, 
                                   donor.nums=donor.nums, 
                                   use.cere.networks=cere.lh.use.regions.7, 
                                   type='all_cere_micros',
                                   net_names=buckner.names$sev,
                                   atlas_field='BucknerMNI152_7')
  
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
  cere.foldchange[['q01_genes']] <- list()
  for ( net in cere.lh.use.regions.7 ) { # rest.networks = region.names, is the array of the 7 networks name
    # print(paste('Getting preferential expression for: ', net, sep = ''))
    
    # negative weight for the contrast matrix, depends on number of comparison networks
    mult.term    <- round(1/(length(region.names)-1),6) # round(a,b),四舍五入 , a 是四舍五入的对象， b是保留的小数位
    o.nets       <- region.names[region.names != net]   # name of the other networks, 
    cur.contrast <- paste('1*', net, '-', mult.term, '*', paste(o.nets, collapse = paste('-', mult.term, '*', sep = '')), sep = '')
    
    # Make the contrast matrix
    cmtx <- NULL
    cmtx <- makeContrasts(contrasts=cur.contrast, levels=colnames(design)) # cmtx is our contrast matrix, eg., for vis, is1*Vis-0.166667*Default-0.166667*Cont-0.166667*Limbic-0.166667*VentAttn-0.166667*DorsAttn-0.166667*SomMot
    tmplm   <- lmFit(expr, design, block=donors, correlation=corfit$consensus.correlation ) # Fit the linear model to the data
    fit     <- eBayes(contrasts.fit( tmplm, cmtx ) ) # contrast.fit: fit the linear model to estimate a set of contrast
    cere.foldchange[['fit_df']][[net]] <- fit        # cere.foldchange has three list: fit_df = fit; stats= ordered topTable(fit, number=Inf); q01_genes = row names of selected stats
    tmp     <- topTable(fit, number=Inf)             # A number of summary statistics are presented by topTable() for the top genes and the selected contrast.
    cere.foldchange[['stats']][[net]] <- tmp[order(rownames(tmp)),] 
    
    # Positive fold change, FDR corrected p<0.01
    pos.idxs       <- which(cere.foldchange$stats[[net]]$logFC > 0)         # The logFC column gives the value of the contrast. Usually this represents a log2-fold change  between two or more experimental 
    adjusted.ps    <- which(cere.foldchange$stats[[net]]$adj.P.Val <= .05)  # adj.P.Value is the p-value adjusted for multiple testing
    genes.tmp      <- rownames(cere.foldchange$stats[[net]])[intersect(adjusted.ps, pos.idxs)]
    cere.foldchange[['q01_genes']][[net]] <- genes.tmp
    sig.genes <- c(sig.genes, genes.tmp)
  }
  
  cere.foldchange.n6 <- cere.foldchange
  cere.n6.sig.genes  <- sig.genes
  cere.n6.sig.genes <- unique(cere.n6.sig.genes)
  
  # save the network-specific genes number in one table
  Sig_cur <- matrix(c(length(cere.foldchange.n6$q01_genes$Default), length(cere.foldchange.n6$q01_genes$Cont), 
                      length(cere.foldchange.n6$q01_genes$Limbic),  length(cere.foldchange.n6$q01_genes$VentAttn), 
                      length(cere.foldchange.n6$q01_genes$DorsAttn), length(cere.foldchange.n6$q01_genes$SomMot),
                      length(cere.foldchange.n6$q01_genes$Vis), 
                      sum(length(cere.foldchange.n6$q01_genes$Default), length(cere.foldchange.n6$q01_genes$Cont), 
                          length(cere.foldchange.n6$q01_genes$Limbic),  length(cere.foldchange.n6$q01_genes$VentAttn), 
                          length(cere.foldchange.n6$q01_genes$DorsAttn), length(cere.foldchange.n6$q01_genes$SomMot),
                          length(cere.foldchange.n6$q01_genes$Vis))), 
                    nrow = 1, ncol = 8)
  colnames(Sig_cur) <- c('Default', 'Cont', 'Limbic', 'VentAttn', 'DorsAttn', 'SomMot', 'Vis', "Total")
  Sig <- rbind(Sig, Sig_cur)
  print(i)
  print("Total_Sig")
  Total_Sig <- Sig_cur[,8]
  print(Total_Sig)
  ## Calculate the gene-fc correlation 
  if (Total_Sig > 1) {  # >3,5,7
    for ( donor in donor.nums[1:2]){ # donor <- "10021"
      cerebellar.atlas <- 'BucknerMNI152_17'
      cerebellar.num   <- 17
      cere.use.regions <- twohemi.net17.regions
      # Mean Normalized cerebellar expression values
      all_data[[donor]]$cere_expr_17         <- averageCereExpr(cere.use.regions, all_data[[donor]], cerebellar.num, 
                                                                all_data[[donor]][[cerebellar.atlas]], 'cere_meanNorm')
      # non-Mean Normalized cerebellar expression values
      # all_data[[donor]]$cere_expr_nonorm_17  <- averageCereExpr(cere.use.regions, all_data[[donor]], cerebellar.num, 
      # all_data[[donor]][[cerebellar.atlas]], 'all_cere_micros')
    }
    # Make the mRNA cerebellum-cerebellum correlation matrix
    gene.idxs <- which(rownames(all_data[['10021']]$cere_expr) %in% unique(cere.n6.sig.genes))
    out.mats.17  <- NULL
    # gene.idxs <- c("AIDA" , "AIF1" )
    # gene.idxs <- "AIDA" 
    # Calculate each gene co-experssion matrix seperately for each subject
    for ( donor in donor.nums[1:2] ){
      # donor <- donor.nums[[1]]
      cur.mat           <- all_data[[donor]]$cere_expr_17[gene.idxs, ] # the difference between cort_expr_17 and cere_expr_17
      # cur.mat           <- all_data[[donor]]$cere_expr_17[gene.idxs[1], ] # the difference between cort_expr_17 and cere_expr_17
      untransformed     <- cor(cur.mat, method = 'spearman')
      # untransformed     <- cor(cur.mat, method = 'pearson')
      # res[[donor]]      <- rcorr(as.matrix(cur.mat), type = "spearman" )
      scaled  <- untransformed # z-tranfsorm spearman correlations
      
      z.corrs <- fisherz(untransformed[upper.tri(untransformed)]) # Upper portion of corr grid
      scaled[upper.tri(untransformed)] <- z.corrs
      
      z.corrs <- fisherz(untransformed[lower.tri(untransformed)])  # Lower portion of corr grid
      scaled[lower.tri(untransformed)] <- z.corrs
      
      # Label the rows/columns
      out.mats.17[[donor]] <- scaled # scaled sig genes expression
      rownames(out.mats.17[[donor]]) <- plot.order.17[twohemi.net17.regions]
      colnames(out.mats.17[[donor]]) <- plot.order.17[twohemi.net17.regions]
    }
    # Average the correlation matrices
    avg.mat   <- (out.mats.17[[1]] + out.mats.17[[2]])/2                            # average the z-transformed correlations of each donors
    rownames(avg.mat) <- plot.order.17[twohemi.net17.regions]
    colnames(avg.mat) <- plot.order.17[twohemi.net17.regions]
    
    # Select the rows that correspond to the mRNA data, oj
    cere2cere.fcmri.17  <- read.csv(paste0(base.dir,'/Data/FC/cere17-cor17-all_fisher/gsr_cerecere.csv'), header = FALSE)
    colnames(cere2cere.fcmri.17) <-  c("VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB", "SalVentAttnA", "SalVentAttnB", "Limbic_Pole", "Limbic_OFC", "ContC", "ContA", "ContB", "TempPar", "DefaultC", "DefaultA", "DefaultB")
    rownames(cere2cere.fcmri.17) <-  c("VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB", "SalVentAttnA", "SalVentAttnB", "Limbic_Pole", "Limbic_OFC", "ContC", "ContA", "ContB", "TempPar", "DefaultC", "DefaultA", "DefaultB")
    
    use_idxs        <- colnames(cere2cere.fcmri.17) %in% rownames(avg.mat)  # avg.mat is the mRNA coorelation martix of two donors which both have left and right hemisphere 49*49
    plot.func       <- cere2cere.fcmri.17[use_idxs, use_idxs]
    
    # Global correspondence of mrna/fcmri, here just use the bi-hemisphere donors: 9861, 10021
    func.1d.arr <- plot.func[upper.tri(plot.func)]
    mrna.1d.arr <- avg.mat[upper.tri(avg.mat)]
    test <- cor.test(func.1d.arr, mrna.1d.arr)
    R_cur <- rbind(test$estimate, test$p.value)
    rownames(R_cur) <- c("R", "p")
    R_random <- cbind(R_random, R_cur)
    print(R_cur)
  }
  else {
    R_random <- cbind(R_random, NA)
  }
}
write.csv(x=R_random, file=paste0(base.dir, '/Output/SupplementaryData/Permutation/Threshold/Permutation_R1.csv')) # same as generated in the cere script
write.csv(x=Sig, file=paste0(base.dir, '/Output/SupplementaryData/Permutation/Threshold/Permutation_SigGene1.csv')) # same as generated in the cere script
# -------------------------------
