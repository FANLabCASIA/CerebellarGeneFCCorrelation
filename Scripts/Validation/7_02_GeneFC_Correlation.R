# -------------------------------
# -------------------------------
# Code to run Gene-FC correlation analyses in:
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

# Prepare before 
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

# Calculate cerebellar differential expression of genes in all 6 donors 
# ----------------
cere.use.regions.7 <- getCereRegions(all_data, donor.nums, atlas.num='7', thresh=2)  # 7 parcel ID (each parcel have some samples) of left hemisphere donor used to define the differential gene 
# cere.use.regions.7 <- rev(buckner.names$sev[cere.use.regions.7+1]) 
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
cere.n6.sig.genes <- unique(cere.n6.sig.genes)
# --------------------------


## Calculate the gene-fc correlation in 7 network
# --------------------------
plot.order.7 <- read.csv(paste0(base.dir, '/Reference_files/7network_names.csv'), header = FALSE)
plot.order.7 <- plot.order.7[-1,]
twohemi.net7.regions <- getCereRegions(all_data, donor.nums[1:2], '7', 2) # 6 region ID 
# Average expression within each parcel that contains data from each bi-hemispheric subject. 
cere.use.regions.7 <- twohemi.net7.regions # 59 parcel ID from both bi-hemispheric AHBA donors, same as the twohemi.net17.regions
for ( donor in donor.nums[1:2]){
  print(paste('Averaging cere data for: ', donor, sep = ''))
  
  # Set the modifiable atlas parameters, for instance you can run cerebellum 17/cortical 17
  cerebellar.atlas <- 'BucknerMNI152_7'
  cerebellar.num   <- 7
  cort.use.regions <- twohemi.net7.regions
  # Mean Normalized cerebellar expression values
  all_data[[donor]]$cere_expr_7         <- averageCereExpr(cort.use.regions, all_data[[donor]], cerebellar.num, 
                                                           all_data[[donor]][[cerebellar.atlas]], 'cere_meanNorm')
  # non-Mean Normalized cerebellar expression values
  all_data[[donor]]$cere_expr_nonorm_7  <- averageCereExpr(cort.use.regions, all_data[[donor]], cerebellar.num, 
                                                           all_data[[donor]][[cerebellar.atlas]], 'all_cere_micros')
} # here for the cere_expr_17, we just get 10 column, which correspond to the 10 region ID derived from the both bi-hemisphere donor
# Network-associated genes used for the cortico-cortical correlations



# Make the mRNA cerebellum-cerebellum correlation matrix
gene.idxs <- which(rownames(all_data[['10021']]$cere_expr_7) %in% unique(cere.n6.sig.genes))
out.mats.7  <- NULL
res <- NULL
# Calculate each gene co-experssion matrix seperately for each subject
for ( donor in donor.nums[1:2] ){
  
  cur.mat           <- all_data[[donor]]$cere_expr_7 [gene.idxs,] # cur.mat  <- all_data[[donor]]$cere_expr_7 value the all 20738 genes
  untransformed     <- cor(cur.mat, method = 'spearman')
  # res[[donor]]      <- rcorr(as.matrix(cur.mat), type = "spearman" )
  
  scaled  <- untransformed # z-tranfsorm spearman correlations
  
  z.corrs <- fisherz(untransformed[upper.tri(untransformed)]) # Upper portion of corr grid
  scaled[upper.tri(untransformed)] <- z.corrs
  
  z.corrs <- fisherz(untransformed[lower.tri(untransformed)])  # Lower portion of corr grid
  scaled[lower.tri(untransformed)] <- z.corrs
  
  # Label the rows/columns
  out.mats.7[[donor]] <- scaled
  rownames(out.mats.7[[donor]]) <- plot.order.7[twohemi.net7.regions]
  colnames(out.mats.7[[donor]]) <- plot.order.7[twohemi.net7.regions]
}
# Average the correlation matrices
avg.mat   <- (out.mats.7[[1]] + out.mats.7[[2]])/2                            # average the z-transformed correlations of each donors
# reg.names <- as.character(atlas.key.17$Name[2:length(atlas.key.17$Name)]) # add region names
rownames(avg.mat) <- plot.order.7[twohemi.net7.regions]
colnames(avg.mat) <- plot.order.7[twohemi.net7.regions]


cere2cere.fcmri.7 <- read.csv(paste0(base.dir,'/Data/FC/cere7-cor17-all_fisher/gsr_cerecere.csv'), header = FALSE)
rownames(cere2cere.fcmri.7) <- plot.order.7
colnames(cere2cere.fcmri.7) <- plot.order.7
use_idxs        <- colnames(cere2cere.fcmri.7) %in% rownames(avg.mat)  # avg.mat is the mRNA coorelation martix of two donors which both have left and right hemisphere 49*49
plot.func       <- cere2cere.fcmri.7[use_idxs, use_idxs]
# plot by corrplot

# Global correspondence of mrna/fcmri, here just use the bi-hemisphere donors: 9861, 10021
func.1d.arr <- plot.func[upper.tri(plot.func)]
mrna.1d.arr <- avg.mat[upper.tri(avg.mat)]
cor.test(func.1d.arr, mrna.1d.arr)

write.csv(x=avg.mat, file=paste0(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet18_Gene_correlation_net7.csv')) # same as generated in the cere script
write.csv(x=plot.func, file=paste0(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet19_FC_WithinCere_net7.csv')) # same as generated in the cere script
# --------------------------