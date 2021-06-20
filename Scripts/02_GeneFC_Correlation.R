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


# 1. Segregation: Network-specific genes
# 1.1  Segregation: Network-specific genes in Cere, n6
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
name.out     <- paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet3_Gene_cere_geneList.csv', sep = '') # same as generated in the cere script
write.table(x=cere.n6.sig.genes, file=name.out, quote=TRUE, col.names=FALSE, row.names=FALSE)



# SUPPLEMENTAL DATA 
# format the differential expression information into one big table
out.matrix  <- NULL
col.headers <- NULL
for ( reg in region.names ) {
  out.dat     <- cbind(cere.foldchange.n6$stats[[reg]]$logFC, cere.foldchange.n6$stats[[reg]]$P.Value, cere.foldchange.n6$stats[[reg]]$adj.P.Val)
  out.matrix  <- cbind(out.matrix, out.dat)
  col.headers <- c(col.headers, c(paste(reg, 'logFC', sep='_'), paste(reg, 'pval', sep='_'), paste(reg, 'adj.pval', sep='_')))
}
colnames(out.matrix) <- col.headers
rownames(out.matrix) <- rownames(cere.foldchange.n6$stats[[reg]])
write.csv(x=out.matrix, file=paste0 (base.dir, '/Output/SupplementaryData/SUPPDATA_sheet1_Gene_cere_diffExpr.csv')) 

# Genes that are significantly differentially expressed in one of the cerebellar networks, n=6
write.me    <- NULL
write.genes <- vector()
write.nums  <- NULL
for (buckner in region.names){
  print(buckner)
  genes.tmp   <- cere.foldchange.n6[['q05_genes']][[buckner]]
  genes.in    <- paste(genes.tmp, collapse = ',')
  write.genes <- rbind.fill.matrix(write.genes, t(genes.tmp)) # t(): change the rows and columns; rbind.fill.matrix: Bind matrices by row, and fill missing columns with NA.
  write.nums  <- c(write.nums, length(genes.tmp))
  write.me    <- rbind(write.me, c(buckner, length(genes.tmp), paste(genes.tmp, collapse = ',')))
}
# Write the significant genes in readable format
write.genes[is.na(write.genes)] <- ''
t.write.genes <- t(write.genes)
write.me      <- rbind(region.names, rbind(write.nums, t.write.genes))
write.csv(x=write.me, file=paste0(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet2_Gene_cere_Network.csv'),row.names=FALSE)
# --------------------------

# 1.2  Segregation: Network-specific genes in Cort, n6
# Preapare befor differentiall expressed analysis
# ----------------
# Read info about split label names for cortex, IDs, and color labels
# this splits the Yeo atlas into 57/114 spatially contiguous regions, depending on whether it's the 7- or 17-network atlas
atlas.key.7  <- read.table(paste0(base.dir, '/Data/Yeo_JNeurophysiol11_SplitLabels/Yeo_JNeurophysiol11_SplitLabels/MNI152/7Networks_ColorLUT_freeview.txt'), col.names=c('ID','Name','R','G','B','A'))
atlas.key.17 <- read.table(paste0(base.dir, '/Data/Yeo_JNeurophysiol11_SplitLabels/Yeo_JNeurophysiol11_SplitLabels/MNI152/17Networks_ColorLUT_freeview.txt'), col.names=c('ID','Name','R','G','B','A'))

# Find cortical parcels that contain a sample in the at least 2 of 6 subjects
use.regions.7 <- getCortRegions(all_data, donor.nums, atlas.num='7', thresh=2)

all_data <- avg_parcel_expression(all_data, 
                                   cerebellar.atlas='BucknerMNI152_7', # We opt for the 7 network parcellation in cerebellum to minimize functional regions with sparse sampling. 
                                   cerebellar.num='sev',
                                   cortical.atlas='splitLabel_7',
                                   cortical.num=51,
                                   cort.use.regions=use.regions.7)

# Array to match each of the 51 spatially contiguous parcels to a Yeo cortical network name (e.g. Parcel 23=Default)
reg.names.27 <- as.character(atlas.key.7$Name[2:27]) # Left hemispheres
reg.names.51 <- as.character(atlas.key.7$Name[2:52]) # both hemispheres
reg2yeo.27   <- matrix(NA,26)
reg2yeo.51   <- matrix(NA,51)
for(buckner in buckner.names[['sev']]){
  cur.idxs     <- grep(buckner, reg.names.27)
  cur.idxs.51  <- grep(buckner, reg.names.51)
  # Treat TempPar as Default, like in original Yeo paper
  if (buckner == 'Default'){
    cur.idxs     <- c(cur.idxs, grep('TempPar', reg.names.27))
    cur.idxs.51  <- c(cur.idxs.51, grep('TempPar', reg.names.51))
  }
  reg2yeo.27[cur.idxs] <- buckner     # give us the parcel ID for each networks only for left hemisphere
  reg2yeo.51[cur.idxs.51] <- buckner
}

# Calculate cortical differential expression across all 6 subjects. 
region.names         <- c('Vis','VentAttn','DorsAttn','SomMot','Limbic','Cont','Default')
use.n6.net7.regions  <- getCortRegions(all_data, donor.nums, '7', 2)
diff.expr.list <- ahba_diff_expression(all_data, 
                                        use.donors=donor.nums, 
                                        use.regions=use.n6.net7.regions,
                                        reg2yeo=reg2yeo.51,
                                        dat2avg='cort_expr_nonorm',
                                        rest.networks=region.names)
cort.foldchanges.n6  <- diff.expr.list[[1]]
cort.n6.sig.genes    <- diff.expr.list[[2]]    # 7689 gene totally, Vis: 3586; VentAttn: 103; DorsAttn: 80; SomMot: 960; Limbic: 2920; Cont: 33; Default: 25.
cort.n6.sig.genes    <- unique(cort.n6.sig.genes ) # 6987 unique genes

write.me    <- NULL
write.genes <- vector()
write.nums  <- NULL
for (buckner in region.names){
  print(buckner)
  genes.tmp   <- cort.foldchanges.n6[['q05_genes']][[buckner]]
  genes.in    <- paste(genes.tmp, collapse = ',')
  write.genes <- rbind.fill.matrix(write.genes, t(genes.tmp)) # t(): change the rows and columns; rbind.fill.matrix: Bind matrices by row, and fill missing columns with NA.
  write.nums  <- c(write.nums, length(genes.tmp))
  write.me    <- rbind(write.me, c(buckner, length(genes.tmp), paste(genes.tmp, collapse = ',')))
}
# Write the significant genes in readable format
write.genes[is.na(write.genes)] <- ''
t.write.genes <- t(write.genes)
write.me      <- rbind(region.names, rbind(write.nums, t.write.genes))
write.csv(x=write.me, file=paste0(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet5_Gene_cort_diffExpr_Network.csv'), row.names=FALSE)

# Write the genes that passed the sig threshold
cort.n6.sig.genes  <- unique(cort.n6.sig.genes )
write.table(x=cort.n6.sig.genes , file=paste0(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet6_Gene_cort_diffExpr_geneList.csv'), row.names=FALSE, col.names=FALSE, quote=TRUE)

# format the differential expression information into one big table
out.matrix  <- NULL
col.headers <- NULL
for ( reg in region.names ) {
  out.dat     <- cbind(cort.foldchanges.n6$stats[[reg]]$logFC, cort.foldchanges.n6$stats[[reg]]$P.Value, cort.foldchanges.n6$stats[[reg]]$adj.P.Val)
  out.matrix  <- cbind(out.matrix, out.dat)
  col.headers <- c(col.headers, c(paste(reg, 'logFC', sep='_'), paste(reg, 'pval', sep='_'), paste(reg, 'adj.pval', sep='_')))
}
colnames(out.matrix) <- col.headers
rownames(out.matrix) <- rownames(cort.foldchanges.n6$stats[[reg]])
write.csv(x=out.matrix, file=paste0(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet4_Gene_cort_diffExpr.csv')) 
# --------------------------

# 1.3 Segregation: Network-specific genes Overlapped
# --------------------------

# Run hypergeometric tests - significant values indicate that a greater than expected number of genes overlap between cortex/cerebellum for a network
# Variable names are meant to mimick the popular example of drawing black/white balls from an urn. 
use.cere.networks  <- c('Default','Cont','Limbic','VentAttn','DorsAttn','SomMot','Vis')
hyper.sig.array <- NULL
write.genes <- vector()
write.nums  <- NULL
for (buckner in use.cere.networks){
  sig.cere.genes.tmp <- cere.foldchange.n6$q05_genes[[buckner]]
  sig.cortical.genes <- cort.foldchanges.n6$q05_genes[[buckner]]
  
  num.overlap        <- length(intersect(sig.cere.genes.tmp, sig.cortical.genes))     # the number of white balls needed every sample time, x
  num.genes          <- length(cort.foldchanges.n6$fit_df$Limbic$coefficients)   
  white.genes        <- length(sig.cere.genes.tmp)                                    # the total number of white balls in the whole urn, a
  black.genes        <- num.genes-length(sig.cere.genes.tmp)                          # N-a
  genes.drawn        <- length(sig.cortical.genes)                                # the number of balls drawn from the urn,n
  print(buckner)
  print(num.overlap)
  
  p.val <- phyper(num.overlap-1, white.genes, black.genes, genes.drawn, lower.tail = F) 
  hyper.sig.array <- c(hyper.sig.array, p.val) 
  
  overlap.genes.tmp <- intersect(sig.cere.genes.tmp, sig.cortical.genes)
  if ( is.null(overlap.genes.tmp) ){
    overlap.genes.tmp <- ''
  }
  write.genes   <- rbind.fill.matrix(write.genes, t(overlap.genes.tmp))
  write.nums    <- c(write.nums, length(overlap.genes.tmp))
}
write.genes[is.na(write.genes)] <- ''
t.write.genes <- t(write.genes)
write.me  <- rbind(use.cere.networks, rbind(write.nums, t.write.genes))
base      <- paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet7_Gene_CerebroCerebellar_Overlap_Network.csv', sep = '')
write.csv(x = write.me, file = base, row.names=FALSE)

# Discrete overlapped genes
# limbic
limbic.genes  <- intersect(cere.foldchange.n6$q05_genes$Limbic, cort.foldchanges.n6$q05_genes$Limbic)
limbic.idxs   <- which(all_data$`9861`$probes_collapse$gene_symbol %in% limbic.genes)
limbic.entrez <- all_data$`9861`$probes_collapse$entrez_id[limbic.idxs]                                        # 56

# Somato/motor loop
sommot.genes  <- intersect(cere.foldchange.n6$q05_genes$SomMot, cort.foldchanges.n6$q05_genes$SomMot)       # 2, "TESC", C6orf167
sommot.idxs   <- which(all_data$`9861`$probes_collapse$gene_symbol %in% sommot.genes)
sommot.entrez <- all_data$`9861`$probes_collapse$entrez_id[sommot.idxs]
# Vis loop
vis.genes  <- intersect(cere.foldchange.n6$q05_genes$Vis, cort.foldchanges.n6$q05_genes$Vis)
vis.idxs   <- which(all_data$`9861`$probes_collapse$gene_symbol %in% vis.genes)
vis.entrez <- all_data$`9861`$probes_collapse$entrez_id[vis.idxs]                                            # 33
# DorsAttn loop
DorsAttn.genes  <- intersect(cere.foldchange.n6$q05_genes$DorsAttn, cort.foldchanges.n6$q05_genes$DorsAttn)
DorsAttn.idxs   <- which(all_data$`9861`$probes_collapse$gene_symbol %in% DorsAttn.genes)
DorsAttn.entrez <- all_data$`9861`$probes_collapse$entrez_id[DorsAttn.idxs]                                  # 0

overlap.genes <- c(limbic.genes, DorsAttn.genes, sommot.genes, vis.genes)  # 91
overlap.genes <- unique(overlap.genes)                                     # 90
write.csv(x = overlap.genes, 
          file = paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet8_Gene_CortCere_geneList.csv', sep = ''), 
          row.names=FALSE)
# Calculate the correlation between the cerebellar and cortical network-specific genes distribution 
pat_cere <- c(0,0,170,0,51,3,221)
pat_cort <- c(25,33,2920,103,80,960,3586)
cor.test(pat_cere, pat_cort)
# ------------------------------------


# 2. Integration
# 2.1 Integration: Intra-cerebellum
# 2.1.1 Gene-FC correlation
# --------------------------

# Order of cerebellar networks for plotting in a correlation grid. 
# 17-network parcels are used To increase resolution of cortico-cortical relationships. 
plot.order.17 <- read.csv(paste0(base.dir, '/Reference_files/17network_names.csv'), header = FALSE)
plot.order.17 <- plot.order.17[-1,]
twohemi.net17.regions <- getCereRegions(all_data, donor.nums[1:2], '17', 2) # 10 parcel ID 

# Average expression within each parcel that contains data from each bi-hemispheric subject. 
cere.use.regions.17 <- twohemi.net17.regions 
for ( donor in donor.nums[1:2]){
  print(paste('Averaging cere data for: ', donor, sep = ''))
  
  # Set the modifiable atlas parameters, for instance you can run cerebellum 17/cortical 17
  cerebellar.atlas <- 'BucknerMNI152_17'
  cerebellar.num   <- 17
  cere.use.regions <- twohemi.net17.regions
  # Mean Normalized cerebellar expression values
  all_data[[donor]]$cere_expr_17         <- averageCereExpr(cere.use.regions, all_data[[donor]], buckner.names[['sevteen']], 
                                                            all_data[[donor]][[cerebellar.atlas]], 'cere_meanNorm')
  # non-Mean Normalized cerebellar expression values
  all_data[[donor]]$cere_expr_nonorm_17  <- averageCereExpr(cere.use.regions, all_data[[donor]], buckner.names[['sevteen']], 
                                                            all_data[[donor]][[cerebellar.atlas]], 'all_cere_micros')
} 

# Network-associated genes used for the cortico-cortical correlations
cere.n6.sig.genes <- unique(cere.n6.sig.genes)

# Make the mRNA cerebellum-cerebellum correlation matrix
gene.idxs <- which(rownames(all_data[['10021']]$cere_expr) %in% unique(cere.n6.sig.genes))
out.mats.17  <- NULL
res <- NULL
# Calculate each gene co-experssion matrix seperately for each subject
for ( donor in donor.nums[1:2] ){
  # donor <- donor.nums[[1]]
  cur.mat           <- all_data[[donor]]$cere_expr_17[gene.idxs, ] # the difference between cort_expr_17 and cere_expr_17
  untransformed     <- cor(cur.mat, method = 'spearman')
  # res[[donor]]      <- rcorr(as.matrix(cur.mat),type ="spearman")  # Significance level
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
#  The correlation significance level of the gene co-expression was evaluated using the overlap 
#  between the correlation matrix for these two individuals and adjusted by Bonferroni correction.
p.adj <- 0.05/45
test1 <- res$`9861`$P> p.adj
test2 <- res$`10021`$P> p.adj
p.mat <- test1 + test2 

# Plot the Cor in Gene and FC
col2 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                           "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                           "#D6604D", "#B2182B", "#67001F"))
diag(avg.mat) <- 0

pdf(file = paste(base.dir, '/Output/SupplementaryData/Fig/corrplot_443genes.pdf', sep = ''),
    width = 10, height = 10)

corrplot(corr = as.matrix(avg.mat), type="upper", method = "color",  # diag=FALSE,
         col = col2(100), tl.pos = "lt", tl.cex = 1.2, tl.col = "black",
         cl.lim = c(-0.75, 0.75))
corrplot(corr = as.matrix(avg.mat), add=TRUE, type="upper", method = "color", diag=FALSE, cl.pos = "n",
         tl.pos = "n", tl.cex = 0.9, tl.col = "black", col = col2(100),
         addCoef.col = "black", number.cex = 1)
corrplot(corr = as.matrix(avg.mat), add=TRUE, type="lower", method="color",diag=FALSE,tl.pos="n", cl.pos="n", 
         col = col2(100), 
         p.mat = as.matrix(p.mat), insig = "label_sig", sig.level = .05, pch.cex = 1, pch.col = "black")
dev.off()
                      
diag(avg.mat) <- 1
write.csv(file = paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet9_Gene_correlations_net17.csv', sep=''), x = avg.mat)

cere2cere.fcmri.17  <- read.csv(paste0(base.dir,'/Data/FC/cere17-cor17-all_fisher/gsr_cerecere.csv'), header = FALSE)

colnames(cere2cere.fcmri.17) <-  c("VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB", "SalVentAttnA", "SalVentAttnB", "Limbic_Pole", "Limbic_OFC", "ContC", "ContA", "ContB", "TempPar", "DefaultC", "DefaultA", "DefaultB")
rownames(cere2cere.fcmri.17) <-  c("VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB", "SalVentAttnA", "SalVentAttnB", "Limbic_Pole", "Limbic_OFC", "ContC", "ContA", "ContB", "TempPar", "DefaultC", "DefaultA", "DefaultB")

use_idxs        <- colnames(cere2cere.fcmri.17) %in% rownames(avg.mat)  # avg.mat is the mRNA coorelation martix of two donors which both have left and right hemisphere 49*49
plot.func       <- cere2cere.fcmri.17[use_idxs, use_idxs]

diag(plot.func) <- 0 
pdf(file = paste(base.dir, "/Output/SupplementaryData/Fig/corrplot_fc_17_fisher.pdf", sep=''),
    width = 10, height = 10)
corrplot(corr = as.matrix(plot.func), type="full", method = "color",  # diag=FALSE,
         col = col2(100), tl.pos = "lt", tl.cex = 1.2, tl.col = "black",
         cl.lim = c(0,0.7))
corrplot(corr = as.matrix(plot.func), add=TRUE, type="upper", method = "color", diag=FALSE, cl.pos = "n",
         tl.pos = "n", tl.cex = 0.9, tl.col = "black", col = col2(100), cl.lim = c(0,0.7),
         addCoef.col = "black", number.cex = 1)
dev.off()

diag(plot.func) <- 1
write.csv(file = paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet10_FC_WithinCere_net17.csv', sep=''), x = plot.func)


# Global correspondence of mrna/fcmri, here just use the bi-hemisphere donors: 9861, 10021
func.1d.arr <- plot.func[upper.tri(plot.func)]
mrna.1d.arr <- avg.mat[upper.tri(avg.mat)]
cor.test(func.1d.arr, mrna.1d.arr)  # 0.4787607
# cor.test(func.1d.arr, mrna.1d.arr, method = 'spearman')

# --------------------------

# 2.1 Integration: Intra-cerebellum
# 2.1.2 Virtual Knock out
# --------------------------

## Rank the sig genes based on the virtual knock-out
cere.n6.sig.genes <- unique(cere.n6.sig.genes)
cor_all <- NULL
cor_dif <- NULL
for (gene_out in cere.n6.sig.genes){
  # Make the mRNA cerebellum-cerebellum correlation matrix
  id_out  <- which(cere.n6.sig.genes == gene_out)
  gene_in <- cere.n6.sig.genes[-id_out]
  gene.idxs <- which(rownames(all_data[['10021']]$cere_expr) %in% gene_in)
  out.mats.17  <- NULL
  # Calculate each gene co-experssion matrix seperately for each subject
  for ( donor in donor.nums[1:2] ){
    # donor <- donor.nums[[1]]
    cur.mat           <- all_data[[donor]]$cere_expr_17[gene.idxs, ] # the difference between cort_expr_17 and cere_expr_17
    untransformed     <- cor(cur.mat, method = 'spearman')
    # untransformed     <- cor(cur.mat, method = 'pearson')
    
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
  use_idxs        <- colnames(cere2cere.fcmri.17) %in% rownames(avg.mat)  # avg.mat is the mRNA coorelation martix of two donors which both have left and right hemisphere 49*49
  plot.func       <- cere2cere.fcmri.17[use_idxs, use_idxs]
  
  # Global correspondence of mrna/fcmri, here just use the bi-hemisphere donors: 9861, 10021
  func.1d.arr <- plot.func[upper.tri(plot.func)]
  mrna.1d.arr <- avg.mat[upper.tri(avg.mat)]
  cor_cur    <- cor.test(func.1d.arr, mrna.1d.arr) 
  # cor_cur    <- cor.test(func.1d.arr, mrna.1d.arr, method = 'spearman') 
  # cor_rp     <- c(cor_cur$estimate, cor_cur$p.value)
  cor_rp_dif  <- c(cor_cur$estimate, cor_cur$estimate-0.4787607, cor_cur$p.value)
  # cor_all <- cbind(cor_all, cor_rp)
  cor_dif <- cbind(cor_dif, cor_rp_dif)
}

colnames(cor_dif) <- cere.n6.sig.genes
rownames(cor_dif) <- c("Rx","GCI","p")
cor_dif<- t(cor_dif)
Entrez_id <- matrix(0,443,1)
cor_dif <- cbind(cor_dif, Entrez_id)
# GCI <- cor_dif[,2]/mean(cor_dif[,2])
# cor_dif<- cbind(cor_dif, GCI)
for (sym in row.names(cor_dif)){
  # sym <- "ABCG2"
  cor_dif[sym,4] <- all_data[["9861"]][["probes_collapse"]]$entrez_id[all_data[["9861"]][["probes_collapse"]]$gene_symbol==sym]
}
colnames(cor_dif) <- c("Rx","GCI","p","Entrez_id")
    
rank_gene <- NULL
rank_gene$In <- cor_dif[which(cor_dif[,2]>0),]
rank_gene$De <- cor_dif[which(cor_dif[,2]<0),]

rank_gene$In <- rank_gene$In[order(-rank_gene$In[,1]),]
rank_gene$De <- rank_gene$De[order(rank_gene$De[,1]),]

write.csv(file = paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet13_Fisher_246_GCI+.csv', sep=''), x = rank_gene$In )
write.csv(file = paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet14_Fisher_197_GCI-.csv', sep=''), x = rank_gene$De )

# --------------------------

# 2.2 Integration: FC-Gene correlation across cerebro-cerebellum circuit
# 2.2 90 overlap genes
# --------------------------

# Determine which cortical parcels are present across in at least 2 donor
n17.use.regions  <- getCortRegions(all_data, donor.nums, '17', 2)
n17.use.regions  <- sort(n17.use.regions)

# Get the average expression of each 17-network parcel, across all six subjects. 
lh.use.regions.17 <- getCortRegions(all_data, donor.nums, atlas.num='17', thresh=2)
for ( donor in donor.nums){
  print(paste('Averaging cort data for: ', donor, sep = ''))
  
  cortical.atlas <- 'splitLabel_17'
  cortical.num   <- 114
  cort.use.regions <- lh.use.regions.17
  # Mean Normalized cortical expression values
  all_data[[donor]]$cort_expr_17_73parcels <- averageCortExpr(cort.use.regions, all_data[[donor]], cortical.num, 
                                                              all_data[[donor]][[cortical.atlas]], 'cort_meanNorm')
}


# Genetic correlation between each cortical and cerebellar region, LH/RH, seperately for each donor
targ.idxs <- which(rownames(all_data[[donor]]$cort_expr) %in% overlap.genes)

cort.use.regions <- getCortRegions(all_data, donor.nums[1:2], atlas.num='17', thresh=2)
plot.order.7 <- read.csv(paste0(base.dir, '/Reference_files/7network_names.csv'), header = FALSE)
plot.order.7 <- plot.order.7[-1,]
twohemi.net7.regions <- getCereRegions(all_data, donor.nums[1:2], '7', 2) # 6 parcel ID 
for ( donor in donor.nums[1:2]){
  print(paste('Averaging cere data for: ', donor, sep = ''))
  
  # Set the modifiable atlas parameters, for instance you can run cerebellum 17/cortical 17
  cerebellar.atlas <- 'BucknerMNI152_7'
  cerebellar.num   <- 7
  cere.use.regions <- twohemi.net7.regions
  # Mean Normalized cerebellar expression values
  all_data[[donor]]$cere_expr_7         <- averageCereExpr(cere.use.regions, all_data[[donor]], buckner.names[['seven']], 
                                                            all_data[[donor]][[cerebellar.atlas]], 'cere_meanNorm')
  # non-Mean Normalized cerebellar expression values
  all_data[[donor]]$cere_expr_nonorm_7  <- averageCereExpr(cere.use.regions, all_data[[donor]], buckner.names[['seven']], 
                                                            all_data[[donor]][[cerebellar.atlas]], 'all_cere_micros')
} 
for (donor in donor.nums[1:2]){ # donor <- "10021"
  n.corts            <- 114 #
  n.ceres            <- dim(all_data[[donor]]$cere_expr_7)[2] # this will vary by subj
  
  # use_regions refers to the cortical areas that have enough data from each subject
  corticocere.cors             <- matrix(NA, n.corts, n.ceres)
  colnames(corticocere.cors)   <- plot.order.7[twohemi.net7.regions]
  
  donor.regions <- getCortRegions(all_data, donor, atlas.num='17', thresh=1)
  cur.expr      <- averageCortExpr(cort.use.regions, all_data[[donor]], cortical.num, all_data[[donor]][[cortical.atlas]], 'cort_meanNorm')
  for ( cort in donor.regions ) { # cort <- 2
    cur.cort <- cur.expr[targ.idxs, cort]
    if ( length(cur.cort[is.na(cur.cort)]) == 0 ) { # make sure the cortical area exists for this subject
      # Iterate over every cerebellum region
      for ( cere in 1:n.ceres ){ # cere <- 1
        cur.cere                       <- all_data[[donor]]$cere_expr_7[targ.idxs, cere]
        corticocere.cors[cort, cere]   <- cor(cur.cort, cur.cere, method = 'spearman')
      }
      
    } else {
      for ( cere in 1:n.ceres ){
        corticocere.cors[cort, cere] <- NA
      }
    }
  }
  all_data[[donor]]$corticocere.cors           <- as.data.frame(corticocere.cors)
  rownames(all_data[[donor]]$corticocere.cors) <- atlas.key.17$Name[2:115]
  rownames(corticocere.cors) <- atlas.key.17$Name[2:115]
  
  # z-transform correlation values
  for (idx in 1:n.ceres){
    corrs   <- corticocere.cors[,idx]
    z.corrs <- fisherz(corrs[!is.na(corrs)])
    corrs[!is.na(corrs)] <- z.corrs
    corticocere.cors[,idx] <- corrs
  }
}

# Average all of the cortico-cerebellum correlations across donors (after z-transform), for plotting
# and comparison with rs-fcMRI data
avg.cort.cere.express <- NULL
min.ceres <- plot.order.7[twohemi.net7.regions]
for ( reg in min.ceres ){ # reg = "SomMot"
  expr.arr <- NULL
  cur.name <- reg
  for ( donor in donor.nums[1:2] ){ # donor <- "10021"
    if (cur.name %in% colnames(all_data[[donor]]$corticocere.cors) ){
      corrs   <- as.matrix(all_data[[donor]]$corticocere.cors[[cur.name]])
      z.corrs <- fisherz(corrs[!is.na(corrs)])
      corrs[!is.na(corrs)] <- z.corrs
      expr.arr <- cbind(expr.arr, as.matrix(all_data[[donor]]$corticocere.cors[[cur.name]]))
    }
  }
  avg.cort.cere.express <- cbind(avg.cort.cere.express, rowMeans(expr.arr, na.rm = TRUE))
}
# Add column names and convert to data.frame
colnames(avg.cort.cere.express) <- min.ceres
avg.cort.cere.express.df        <- as.data.frame(avg.cort.cere.express)
rownames(avg.cort.cere.express.df) <- atlas.key.17$Name[2:115]

cere2cort.fcmri     <- read.csv(paste0(base.dir,'/Data/FC/cere7-cor17-all_fisher/gsr_corcere.csv'))
rownames(cere2cort.fcmri)    <- atlas.key.17$Name[2:115]

# Get rid of missing rows in fc and gene
avg.cort.cere.express.df.nonan  <- avg.cort.cere.express.df[which(!is.na(avg.cort.cere.express.df$SomMot)),]
# colnames(cere2cort.fcmri)   <- rownames(cere2cort.fcmri)
cere2cort.fcmri.nonan           <- cere2cort.fcmri[which(!is.na(avg.cort.cere.express.df$SomMot)),]
cere2cort.fcmri.nonan           <- cere2cort.fcmri.nonan[, colnames(cere2cort.fcmri.nonan)%in%colnames(avg.cort.cere.express.df.nonan)]
rownames(cere2cort.fcmri.nonan) <- rownames(avg.cort.cere.express.df.nonan)

# Get the Correlation for each pair of network in Gene, and Gene-FC coupling
Cor <- NULL
GFS <- cor.test(avg.cort.cere.express.df.nonan$SomMot, cere2cort.fcmri.nonan$SomMot)
GFD <- cor.test(avg.cort.cere.express.df.nonan$DorsAttn, cere2cort.fcmri.nonan$DorsAttn)        
GFV <- cor.test(avg.cort.cere.express.df.nonan$VentAttn, cere2cort.fcmri.nonan$VentAttn) 
GFL <- cor.test(avg.cort.cere.express.df.nonan$Limbic, cere2cort.fcmri.nonan$Limbic)    # R = 0.3601294, p =  0.005083
GFC <- cor.test(avg.cort.cere.express.df.nonan$Cont , cere2cort.fcmri.nonan$Cont)       # R = -0.3259672, p = 0.01175
GFDe <- cor.test(avg.cort.cere.express.df.nonan$Default, cere2cort.fcmri.nonan$Default)
# FDR
Cor$Gene_FC <- cbind(rbind(GFS$estimate, GFS$p.value), rbind(GFD$estimate, GFD$p.value),
                     rbind(GFV$estimate, GFV$p.value), rbind(GFL$estimate, GFL$p.value),
                     rbind(GFC$estimate, GFC$p.value), rbind(GFDe$estimate, GFDe$p.value))
colnames(Cor$Gene_FC ) <- plot.order.7[twohemi.net7.regions]
rownames(Cor$Gene_FC)  <- c("R", "p")
Cor$Gene_FC <- t(Cor$Gene_FC)  
Cor$Gene_FC <- cbind(Cor$Gene_FC, 
                     p.adjust(Cor$Gene_FC[,2] , method = "BH", n = length(Cor$Gene_FC[,2])))
colnames(Cor$Gene_FC)  <- c("R", "p", "adj_p")


g <- NULL
g$Limbic_SomMot   <- cor.test(avg.cort.cere.express.df.nonan$Limbic, avg.cort.cere.express.df.nonan$SomMot)
g$Limbic_VentAttn <- cor.test(avg.cort.cere.express.df.nonan$Limbic, avg.cort.cere.express.df.nonan$VentAttn)
g$Limbic_Default <- cor.test(avg.cort.cere.express.df.nonan$Limbic, avg.cort.cere.express.df.nonan$Default)
g$Limbic_DorsAttn <- cor.test(avg.cort.cere.express.df.nonan$Limbic, avg.cort.cere.express.df.nonan$DorsAttn)
g$Limbic_Cont <- cor.test(avg.cort.cere.express.df.nonan$Limbic, avg.cort.cere.express.df.nonan$Cont)

g$Cont_SomMot <- cor.test(avg.cort.cere.express.df.nonan$Cont, avg.cort.cere.express.df.nonan$SomMot)
g$Cont_VentAttn <- cor.test(avg.cort.cere.express.df.nonan$Cont, avg.cort.cere.express.df.nonan$VentAttn)
g$Cont_Default <- cor.test(avg.cort.cere.express.df.nonan$Cont, avg.cort.cere.express.df.nonan$Default)
g$Cont_DorsAttn <- cor.test(avg.cort.cere.express.df.nonan$Cont, avg.cort.cere.express.df.nonan$DorsAttn)


Cor$Gene <- cbind(rbind(g$Limbic_SomMot$estimate, g$Limbic_SomMot$p.value), 
                  rbind(g$Limbic_VentAttn$estimate, g$Limbic_VentAttn$p.value),
                  rbind(g$Limbic_Default$estimate, g$Limbic_Default$p.value), 
                  rbind(g$Limbic_DorsAttn$estimate, g$Limbic_DorsAttn$p.value),
                  rbind(g$Limbic_Cont$estimate, g$Limbic_Cont$p.value), 
                  rbind(g$Cont_SomMot$estimate, g$Cont_SomMot$p.value),
                  rbind(g$Cont_VentAttn$estimate, g$Cont_VentAttn$p.value),
                  rbind(g$Cont_Default$estimate, g$Cont_Default$p.value),
                  rbind(g$Cont_DorsAttn$estimate, g$Cont_DorsAttn$p.value)
                  )
colnames(Cor$Gene ) <- names(g)
rownames(Cor$Gene)  <- c("R", "p")
Cor$Gene <- t(Cor$Gene)  
Cor$Gene <- cbind(Cor$Gene, 
                     p.adjust(Cor$Gene[,2] , method = "BH", n = length(Cor$Gene[,2])))
colnames(Cor$Gene)  <- c("R", "p", "adj_p")


f <- NULL
f$Limbic_SomMot   <- cor.test(cere2cort.fcmri.nonan$Limbic, cere2cort.fcmri.nonan$SomMot)
f$Limbic_VentAttn <- cor.test(cere2cort.fcmri.nonan$Limbic, cere2cort.fcmri.nonan$VentAttn)
f$Limbic_Default <- cor.test(cere2cort.fcmri.nonan$Limbic, cere2cort.fcmri.nonan$Default)
f$Limbic_DorsAttn <- cor.test(cere2cort.fcmri.nonan$Limbic, cere2cort.fcmri.nonan$DorsAttn)
f$Limbic_Cont <- cor.test(cere2cort.fcmri.nonan$Limbic, cere2cort.fcmri.nonan$Cont)

f$Cont_SomMot <- cor.test(cere2cort.fcmri.nonan$Cont, cere2cort.fcmri.nonan$SomMot)
f$Cont_VentAttn <- cor.test(cere2cort.fcmri.nonan$Cont, cere2cort.fcmri.nonan$VentAttn)
f$Cont_Default <- cor.test(cere2cort.fcmri.nonan$Cont, cere2cort.fcmri.nonan$Default)
f$Cont_DorsAttn <- cor.test(cere2cort.fcmri.nonan$Cont, cere2cort.fcmri.nonan$DorsAttn)


Cor$FC <- cbind(rbind(f$Limbic_SomMot$estimate, f$Limbic_SomMot$p.value), 
                  rbind(f$Limbic_VentAttn$estimate, f$Limbic_VentAttn$p.value),
                  rbind(f$Limbic_Default$estimate, f$Limbic_Default$p.value), 
                  rbind(f$Limbic_DorsAttn$estimate, f$Limbic_DorsAttn$p.value),
                  rbind(f$Limbic_Cont$estimate, f$Limbic_Cont$p.value), 
                  rbind(f$Cont_SomMot$estimate, f$Cont_SomMot$p.value),
                  rbind(f$Cont_VentAttn$estimate, f$Cont_VentAttn$p.value),
                  rbind(f$Cont_Default$estimate, f$Cont_Default$p.value),
                  rbind(f$Cont_DorsAttn$estimate, f$Cont_DorsAttn$p.value)
)
colnames(Cor$FC ) <- names(f)
rownames(Cor$FC)  <- c("R", "p")
Cor$FC <- t(Cor$FC)  
Cor$FC <- cbind(Cor$FC, 
                  p.adjust(Cor$FC[,2] , method = "BH", n = length(Cor$FC[,2])))
colnames(Cor$FC)  <- c("R", "p", "adj_p")

filename = paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet11_Gene_correlations_CortCere_90genes_cere7_cort17.csv', sep = '')
write.csv(x=avg.cort.cere.express.df.nonan, file = filename)
filename = paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet12_FC_CortCere_cere7_cort17.csv', sep = '')
write.csv(x = cere2cort.fcmri.nonan, file = filename)


filename = paste(base.dir, '/Output/SupplementaryData/Self/cere_cort/cere7_cort17_GeneFC.csv', sep = '')
write.csv(x = Cor$Gene_FC, file = filename)
filename = paste(base.dir, '/Output/SupplementaryData/Self/cere_cort/cere7_cort17_Gene.csv', sep = '')
write.csv(x = Cor$Gene, file = filename)
filename = paste(base.dir, '/Output/SupplementaryData/Self/cere_cort/cere7_cort17_FC.csv', sep = '')
write.csv(x = Cor$FC, file = filename)

# Partial correlation between genetic and functional correlation between limbic and control
# BiocManager::install("ggm")
library(ggm)
test <- cbind(avg.cort.cere.express.df.nonan$Cont, cere2cort.fcmri.nonan$Cont, avg.cort.cere.express.df.nonan$Limbic)
colnames(test)<- c("Cont_Gene", "Cont_FC", "Limbic_Gene")
pc <- pcor(c("Cont_Gene", "Cont_FC", "Limbic_Gene"),cov(test))
pcor.test(pc, 1, 59)
# -------------------------- 
