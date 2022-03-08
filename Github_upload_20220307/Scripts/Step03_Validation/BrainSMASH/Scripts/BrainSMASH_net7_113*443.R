# -------------------------------
# -------------------------------
# Code to run the BrainSMASH test for Gene-FC correlation analyses in:
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

# Calculate cerebellar differential expression of genes in all 6 donors 
cere.use.regions.7 <- getCereRegions(all_data, donor.nums, atlas.num='7', thresh=2)  # 7 parcel ID (each parcel have some samples) of left hemisphere donor used to define the differential gene 
cere.use.regions.7 <- c('Default', 'Cont', 'Limbic', 'VentAttn', 'DorsAttn', 'SomMot', 'Vis') # because for all cere networks, they have more than 2 donors samples, which get from before line

# Get Average cerebellar expression for each subject, in 6 subjects. 
region.names         <- c('Vis','VentAttn','DorsAttn','SomMot','Limbic','Cont','Default')
out <- averageWithinCereNetworks(all_data=all_data, 
                                 donor.nums=donor.nums,  
                                 use.cere.networks=cere.use.regions.7, 
                                 type='all_cere_micros',
                                 net_names=buckner.names$sev,
                                 atlas_field='BucknerMNI152_7')

expr    <- out[[1]]  
regions <- out[[2]]  
donors  <- out[[3]]  
colnames(expr) <- paste(donors, regions, sep = '_')

# Use limma to calculate differential expression for each network, relative to all others
fac              <- as.factor(regions)                
design           <- model.matrix(~0 + fac)            
colnames(design) <- gsub('fac', '', colnames(design)) 
corfit           <- duplicateCorrelation(expr, design, block=donors) 

# Calculate Differential expression for each region
sig.genes       <- NULL
cere.foldchange <- NULL
cere.foldchange[['q05_genes']] <- list()
for ( net in cere.use.regions.7 ) { 
  print(paste('Getting preferential expression for: ', net, sep = ''))
  
  # Negative weight for the contrast matrix, depends on number of comparison networks
  mult.term    <- round(1/(length(region.names)-1),6)
  o.nets       <- region.names[region.names != net]
  cur.contrast <- paste('1*', net, '-', mult.term, '*', paste(o.nets, collapse = paste('-', mult.term, '*', sep = '')), sep = '')
  
  # Make the contrast matrix
  cmtx <- NULL
  cmtx <- makeContrasts(contrasts=cur.contrast, levels=colnames(design)) 
  tmplm   <- lmFit(expr, design, block=donors, correlation=corfit$consensus.correlation )
  fit     <- eBayes(contrasts.fit( tmplm, cmtx ) ) 
  cere.foldchange[['fit_df']][[net]] <- fit        
  tmp     <- topTable(fit, number=Inf)            
  cere.foldchange[['stats']][[net]] <- tmp[order(rownames(tmp)),] 
  
  # Positive fold change, FDR corrected p<0.01
  pos.idxs       <- which(cere.foldchange$stats[[net]]$logFC > 0)         
  adjusted.ps    <- which(cere.foldchange$stats[[net]]$adj.P.Val <= .05) 
  genes.tmp      <- rownames(cere.foldchange$stats[[net]])[intersect(adjusted.ps, pos.idxs)]
  cere.foldchange[['q05_genes']][[net]] <- genes.tmp
  print(length(genes.tmp)) 
  sig.genes <- c(sig.genes, genes.tmp)
}
cere.foldchange.n6 <- cere.foldchange
cere.n6.sig.genes  <- sig.genes

plot.order.7 <- read.csv(paste0(base.dir, '/Reference_files/7network_names.csv'), header = FALSE)
plot.order.7 <- plot.order.7[-1,]
twohemi.net7.regions <- getCereRegions(all_data, donor.nums[1:2], '7', 2) # 6 region ID 
# Average expression within each parcel that contains data from each bi-hemispheric subject. 
cere.use.regions.7 <- twohemi.net7.regions 
# ----------------


## 1. Write out the sample * sig gene to "SA_control_permutation.py"
# ----------------
cere.use.regions.7 <- twohemi.net7.regions 
Lab <- NULL
all.expr <- NULL
for ( donor in donor.nums[1:2]){ # donor = "9861"
  print(paste('Averaging cere data for: ', donor, sep = ''))
  
  # Set the modifiable atlas parameters, for instance you can run cerebellum 17/cortical 17
  cerebellar.atlas <- 'BucknerMNI152_7'
  cerebellar.num   <- 7
  cere.use.regions <- twohemi.net7.regions
  bucknernames = as.character(read.csv(paste0(base.dir, '/Reference_files/7network_names.csv'), header = FALSE)$V1)
  # Mean Normalized cerebellar expression values
  all_data[[donor]]$cere_expr_7         <- averageCereExpr(cere.use.regions, all_data[[donor]], bucknernames, 
                                                            all_data[[donor]][[cerebellar.atlas]], 'cere_meanNorm')
  
  # Write out the sample annotation
  cere.regions <- cere.use.regions
  data.struct = all_data[[donor]]
  
  samp.labels = all_data[[donor]][[cerebellar.atlas]]
  data.type = 'cere_meanNorm'
  
  
  for (reg in cere.regions) { # reg = 3
    # name array to add as column headers at the finish
    col_idx  <- samp.labels == reg                     # indices for this dth subregion
    reg_data <- data.struct[[data.type]][, col_idx]
    
    ann_cur <- all_data[[donor]]$cere_samples[col_idx, ]
    cname    <- bucknernames[as.integer(reg) + 1]     # so here the bucknernames can contains the 1 None
    if (sum(col_idx)==1) {
      
      lab <- data.frame(header = colnames(data.struct[[data.type]])[col_idx],
                        name = cname,
                        label = reg,
                        donor = donor)
      lab <- cbind(lab, ann_cur)
    } else{
      count_cur <- sum(col_idx)
      lab <- data.frame(header = colnames(reg_data),
                        name = rep(cname, count_cur),
                        label = rep(reg, count_cur),
                        donor = rep(donor, count_cur))
      lab <- cbind(lab, ann_cur)
    }
    
    Lab <- rbind(Lab, lab)
    all.expr <- cbind(all.expr, reg_data)
    
  }
  
  #Write out the 443
  gene.idxs <- which(rownames(all_data[['10021']]$cere_expr) %in% unique(cere.n6.sig.genes))
  mat <- all.expr[gene.idxs, ]
} 
write.csv(x = Lab, file = paste0(base.dir, '/Scripts/Step03_Validation/BrainSMASH/Results/7_network/Annotation.csv'))
write.csv(x = mat, file = paste0(base.dir, '/Scripts/Step03_Validation/BrainSMASH/Results/7_network/Mat.csv'))
# ----------------


## 2. Run the "SA_control_permutation.py"


## 3. Calculate the Gene-FC correlation using the permutated genetic correlation matrix
# ----------------
datalist <- NULL
for ( donor in donor.nums[1:2] ){
  setwd(paste0('/n02dat01/users/lchai/BrainsMASH_yp/Net_7/output_dir/Permutationmap_10000/', donor, '/'))
  filelist <- list.files(pattern=".*.csv")
  m <-length(filelist)
  m
  datalist[[donor]] <- lapply(filelist, function(x) read.csv(x,header=F,stringsAsFactors=F))
  # datalist[[donor]] <- lapply(filelist[505], function(x) read.csv(x,header=F,stringsAsFactors=F))
}


# BrainMesh Permutation: Read in data
datalist2 <- NULL
R_random <- NULL
for (i in 1:m) {
  out.mats.7 <- NULL
  for ( donor in donor.nums[1:2] ){ # donor <- donor.nums[[1]]
    datalist2[[donor]][[i]] <- data.frame(datalist[[donor]][[i]])
    cur.mat           <- t(data.frame(datalist2[[donor]][[i]])) # the difference between cort_expr_17 and cere_expr_17
    untransformed     <- cor(cur.mat, method = 'spearman')
    scaled  <- untransformed # z-tranfsorm spearman correlations
    
    z.corrs <- fisherz(untransformed[upper.tri(untransformed)]) # Upper portion of corr grid
    scaled[upper.tri(untransformed)] <- z.corrs
    
    z.corrs <- fisherz(untransformed[lower.tri(untransformed)])  # Lower portion of corr grid
    scaled[lower.tri(untransformed)] <- z.corrs
    
    # Label the rows/columns
    out.mats.7[[donor]] <- scaled # scaled sig genes expression
    rownames(out.mats.7[[donor]]) <- plot.order.7[twohemi.net7.regions]
    colnames(out.mats.7[[donor]]) <- plot.order.7[twohemi.net7.regions]
  }
  avg.mat   <- (out.mats.7[[1]] + out.mats.7[[2]])/2                            # average the z-transformed correlations of each donors
  rownames(avg.mat) <- plot.order.7[twohemi.net7.regions]
  colnames(avg.mat) <- plot.order.7[twohemi.net7.regions]
  
  cere2cere.fcmri.7 <- read.csv(paste0(base.dir,'/Data/FC/cere7-cor17-all_fisher/gsr_cerecere.csv'), header = FALSE)
  rownames(cere2cere.fcmri.7) <- plot.order.7
  colnames(cere2cere.fcmri.7) <- plot.order.7
  use_idxs        <- colnames(cere2cere.fcmri.7) %in% rownames(avg.mat)  # avg.mat is the mRNA coorelation martix of two donors which both have left and right hemisphere 49*49
  plot.func       <- cere2cere.fcmri.7[use_idxs, use_idxs]
  
  func.1d.arr <- plot.func[upper.tri(plot.func)]
  mrna.1d.arr <- avg.mat[upper.tri(avg.mat)]
  R <- cor.test(func.1d.arr, mrna.1d.arr)
  R_cur <- rbind(R$estimate, R$p.value)
  R_random <- cbind(R_random, R_cur)
}
length(R_random[R_random[1,]>0.7645752])/length(R_random[1,])  # 0.0128
sum(R_random[1,]>0.7645752)/length(R_random[1,])               # 0.0064
# ----------------
