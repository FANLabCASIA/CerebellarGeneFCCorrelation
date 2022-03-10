# -------------------------------
# -------------------------------
# Code to validate the stability of network-specific genes in:
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
base.dir <- ''     # Enter your path 
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


# 1. Segregation: Validate using 4 left donors
# Calculate cerebellar differential expression of genes in all 6 donors 
# ----------------
cere.use.regions.7 <- getCereRegions(all_data, donor.nums[3:6], atlas.num='7', thresh=2)  # 7 parcel ID (each parcel have some samples) of left hemisphere donor used to define the differential gene 
cere.use.regions.7 <- c('Default', 'Cont', 'Limbic', 'VentAttn', 'DorsAttn', 'SomMot', 'Vis') # because for all cere networks, they have more than 2 donors samples, which get from before line

# Get Average cerebellar expression for each subject, in 6 subjects. 
region.names         <- c('Vis','VentAttn','DorsAttn','SomMot','Limbic','Cont','Default')
out <- averageWithinCereNetworks(all_data=all_data, 
                                 donor.nums=donor.nums[3:6], 
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
cere.foldchange.n4 <- cere.foldchange
cere.n4.sig.genes  <- sig.genes
cere.n4.sig.genes <- unique(cere.n4.sig.genes)
# ----------------

# 2. Overlap with main strategy
# ----------------
sig.All <- read.csv(paste(base.dir,"/Output/SupplementaryData/SUPPDATA_sheet3_Gene_cere_GeneList.csv", sep = ""), header = FALSE)
sig.overlap <- intersect(cere.n4.sig.genes, as.matrix(sig.All)) 
# write.csv(file = paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet5_Gene_cere_Network_4donor_Overlap.csv', sep=''), x = sig.overlap )

# name.out     <- paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet5_Gene_cere_Network_4donor_geneList.csv', sep = '') # same as generated in the cere script
# write.table(x=cere.n4.sig.genes, file=name.out, quote=TRUE, col.names=FALSE, row.names=FALSE)

gene.net <- read.csv(paste(base.dir,"/Output/SupplementaryData/SUPPDATA_sheet2_Gene_cere_Network.csv", sep = ""), header = TRUE)
colnames(gene.net) <- as.matrix(gene.net[1,])
gene.net <- as.matrix(gene.net)
gene.net <- gene.net[-c(1,2),]
length(intersect(cere.n4.sig.genes, as.matrix(sig.All)))               # 141
length(intersect(cere.foldchange.n4$q05_genes$Vis, gene.net[,"Vis"]))       # 97, 
length(intersect(cere.foldchange.n4$q05_genes$Limbic, gene.net[,"Limbic"]))    # 44, 
length(intersect(cere.foldchange.n4$q05_genes$DorsAttn, gene.net[,"DorsAttn"]))  # 0, 
length(intersect(cere.foldchange.n4$q05_genes$SomMot, gene.net[,"SomMot"]))  # 0

# SUPPLEMENTAL DATA 
# format the differential expression information into one big table
out.matrix  <- NULL
col.headers <- NULL
for ( reg in region.names ) {
  out.dat     <- cbind(cere.foldchange.n4$stats[[reg]]$logFC, cere.foldchange.n4$stats[[reg]]$P.Value, cere.foldchange.n4$stats[[reg]]$adj.P.Val)
  out.matrix  <- cbind(out.matrix, out.dat)
  col.headers <- c(col.headers, c(paste(reg, 'logFC', sep='_'), paste(reg, 'pval', sep='_'), paste(reg, 'adj.pval', sep='_')))
}
colnames(out.matrix) <- col.headers
rownames(out.matrix) <- rownames(cere.foldchange.n4$stats[[reg]])
# write.csv(x=out.matrix, file=paste0 (base.dir, '/Output/SupplementaryData/SUPPDATA_sheet5_Gene_cere_Network_4donor_diffExpr.csv')) 

# Genes that are significantly differentially expressed in one of the cerebellar networks, n=6
write.me    <- NULL
write.genes <- vector()
write.nums  <- NULL
for (buckner in region.names){
  print(buckner)
  genes.tmp   <- cere.foldchange.n4[['q05_genes']][[buckner]]
  genes.in    <- paste(genes.tmp, collapse = ',')
  write.genes <- rbind.fill.matrix(write.genes, t(genes.tmp)) # t(): change the rows and columns; rbind.fill.matrix: Bind matrices by row, and fill missing columns with NA.
  write.nums  <- c(write.nums, length(genes.tmp))
  write.me    <- rbind(write.me, c(buckner, length(genes.tmp), paste(genes.tmp, collapse = ',')))
}
# Write the significant genes in readable format
write.genes[is.na(write.genes)] <- ''
t.write.genes <- t(write.genes)
write.me      <- rbind(region.names, rbind(write.nums, t.write.genes))
write.csv(x=write.me, file=paste0(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet5_Gene_cere_Network_4donor.csv'),row.names=FALSE)
# ----------------

