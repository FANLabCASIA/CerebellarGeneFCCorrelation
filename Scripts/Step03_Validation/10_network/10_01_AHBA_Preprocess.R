# -------------------------------
# -------------------------------
# Code to run AHBA pre-process in MDTB atlas:
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

## R package needed
# ----------------
# Use the BiocManager::install() command if they haven't been installed
library(data.table) 
library(WGCNA) 
library(limma) 
# ----------------


## Modify these filepaths for your local directory structure
# ----------------
afni.dir <- ' '                        # enter the abin directory for local install of AFNI
base.dir <- ' '          # enter the directory containing the script repository
source(paste(base.dir, '/Scripts/function_library.R', sep = ''))   # Source function library for this project
data_path  <- paste(base.dir, '/Data/AHBA/AHBA_original_data', sep = '') # path to AHBA microarray data
donor.nums <- c('9861', '10021', '12876', '14380', '15496', '15697')
filenames  <- c('donor_10021', 'donor_9861', 
                'donor_12876', 'donor_15697',
                'donor_14380', 'donor_15496')
# ----------------


## Read/process AHBA microarray data from each donor
# ----------------
all_data <- NULL
for ( donor in donor.nums ) {
  file <- paste('donor_', donor, sep='')
  print(paste('Reading and collapsing data for donor: ', donor, sep = ''))
  
  # Initiate empty field for this donor's data
  all_data[[donor]] <- NULL
  
  
  # Sample Information
  saname    <- paste(data_path, file, 'SampleAnnot.csv', sep='/')
  samp_info <- read.csv(saname)
  all_data[[donor]]$raw_samp <- samp_info
  
  
  # Read Microexpression Data
  fname                <- paste(data_path, file, 'MicroarrayExpression.csv', sep='/')
  microdata            <- fread(fname, header = F, sep = ',')
  micro_arr            <- as.matrix(microdata)
  micro_temp           <- micro_arr[,2:dim(micro_arr)[2]] # First column contains probe IDs
  rownames(micro_temp) <- micro_arr[,1]
  micro_df             <- as.data.frame(micro_temp)
  all_data[[donor]]$raw_micros <- micro_df
  
  
  # Read PA-Call (signal present vs absent)
  paname      <- paste(data_path, file, 'PACall.csv', sep='/')
  pacall      <- fread(paname, header = F, sep = ',')
  pacall_arr  <- as.matrix(pacall)
  pa_dat      <- pacall_arr[,2:dim(pacall_arr)[2]]
  rownames(pa_dat) <- pacall_arr[,1]
  all_data[[donor]]$raw_pas <- pa_dat
  
  
  # Information about each Gene Probe (~50,000 probes for ~20,000 genes)
  pname     <- paste(data_path, file, 'Probes.csv', sep='/')
  probes    <- read.csv(pname)
  all_data[[donor]]$raw_probes <- probes
  
  
  # Gene Ontology Info
  oname    <- paste(data_path, file, 'Ontology.csv', sep='/')
  ont_data <- read.csv(oname)
  all_data[[donor]]$raw_ont <- ont_data
  
  
  # output filenames for collapseRows
  write_micro_name   <- paste(data_path, '/', file, '/collapsed_micro_', donor, '.csv', sep="")
  write_select_name  <- paste(data_path, '/', file, '/selectedRows_', donor, '.csv', sep="")
  write_grp2row_name <- paste(data_path, '/', file, '/grp2row_', donor,'.csv', sep="")
  
  # discard probes without an entrez-id
  num_samples   <- dim(all_data[[donor]]$raw_micros)[2]
  trash_me      <- is.na(all_data[[donor]]$raw_probes$entrez_id) # & pasums > cutoff
  good_probes   <- all_data[[donor]]$raw_probes[trash_me == FALSE,]  
  
  # Select the probes with valid Entrez IDs
  all_data[[donor]]$probes_filter <- good_probes
  all_data[[donor]]$micro_filter  <- all_data[[donor]]$raw_micros[rownames(all_data[[donor]]$raw_micros) %in% good_probes$probe_id,]
  all_data[[donor]]$pas_filter    <- all_data[[donor]]$raw_pas[rownames(all_data[[donor]]$raw_pas) %in% good_probes$probe_id,]
  
  
  # Collapse Rows Function selects the "best"/maxMean probe for each gene, WCGNA toolbox
  out     <- collapseRows(all_data[[donor]]$micro_filter, 
                          all_data[[donor]]$probes_filter$gene_symbol, 
                          rownames(all_data[[donor]]$micro_filter), 
                          method='MaxMean', 
                          connectivityBasedCollapsing = TRUE)
  
  # Write collapsed data
  all_data[[donor]]$micro_collapsed <- out$datETcollapsed
  write.csv(x = all_data[[donor]]$micro_collapsed, file = write_micro_name)
  
  # Write group2row
  all_data[[donor]]$group2row  <- out$group2row
  write.csv(x = all_data[[donor]]$group2row, file = write_grp2row_name)
  
  # Write selected Rows
  all_data[[donor]]$selectedRow  <- out$selectedRow
  write.csv(x = all_data[[donor]]$selectedRow, file = write_select_name)
  
  # Write collapsed probes
  write_probe_name                  <- paste(data_path, '/', file, '/collapsed_probes_', donor, '.csv', sep="")
  all_data[[donor]]$probes_collapse <- all_data[[donor]]$probes_filter[out$selectedRow,]
  write.csv(x = all_data[[donor]]$probes_collapse, file = write_probe_name)
  
  # Write collapsed PA calls
  write_PAcall_name                  <- paste(data_path, '/', file, '/collapsed_PAcall_', donor, '.csv', sep="")
  all_data[[donor]]$pas_collapse     <- all_data[[donor]]$pas_filter[out$selectedRow,]
  write.csv(x = all_data[[donor]]$pas_collapse, file = write_PAcall_name)
}
# save(file=paste(base.dir, '/Output/RData/all_data.Rdata', sep=''), x=all_data)
# load(file=paste(base.dir, '/Output/RData/all_data.Rdata', sep=''))
# --------------------- 


## Calculate the mean normalized expression for Cerebellum and Cortex seperetedly
# ID numbers for Cortex + Cerebellum structures
# ----------------
cort_in      <- read.csv(paste(base.dir, '/Reference_files/cort_regions.csv', sep=''), header=FALSE)
cortex       <- as.numeric(cort_in$V1)

# Get info about cortical/cerebellum samples
for ( donor in donor.nums ){
  print(paste('Retrieving cortical and cerebellum information for: ', donor, sep = ''))
  
  # Get expression and sample information about cortical samples
  cort_idxs                          <- all_data[[donor]]$raw_samp$structure_id %in% cortex
  all_data[[donor]]$cort_samples     <- all_data[[donor]]$raw_samp[cort_idxs,]                        # select info from cerebral cortex samples
  all_data[[donor]]$cort_acros       <- factor(all_data[[donor]]$cort_samples$structure_acronym)      # list of corresponding struct acronyms
  all_data[[donor]]$all_cort_micros  <- as.data.frame(all_data[[donor]]$micro_collapsed)[, cort_idxs]
  
  
  # Cere - ID numbers for cerebellum structures
  test           <- all_data[[donor]]$raw_samp
  cere_in        <- test[which(test[,"slab_type"]=="CB"),] 
  cere_in_wn     <- cere_in[!grepl("nucleus", cere_in$structure_name),]
  reference_path <- paste(base.dir, '/Reference_files/', sep = '')
  write.csv(x = cere_in_wn,file=paste( reference_path, 'cere_regions_', donor, '.csv', sep=''))
  cere           <- as.numeric(cere_in_wn$structure_id)
  cere_idxs                              <- all_data[[donor]]$raw_samp$structure_id %in% cere
  all_data[[donor]]$cere_samples         <- all_data[[donor]]$raw_samp[cere_idxs,]                    # select info from striatal samples
  all_data[[donor]]$all_cere_micros      <- as.data.frame(all_data[[donor]]$micro_collapsed)[, cere_idxs]
  
  
  # Mean Normalize expression values
  all_data[[donor]]$ceres_meanNorm     <- t(apply(all_data[[donor]]$all_cere_micros, 1, function(x) x-mean(x) ))
  all_data[[donor]]$corts_meanNorm     <- t(apply(all_data[[donor]]$all_cort_micros, 1, function(x) x-mean(x) ))
}
# ----------------


## Identify which samples overlap with cere 10 network Atlases, using AHBA provided MNI locations, 
# this step will take a long time: 50 min per donor, so 50*6=300min/60min/h=5h in 200
# ----------------
atlas_dir <- paste(base.dir, '/Output/Atlas_overlap/', sep = '')
for ( donor in donor.nums ){
  atlas_names       <- 'MDTB_10Regions_MNI_2MM.nii'  
  atlas             <- paste(base.dir, '/Reference_files/MDTB_10Regions_2MM.nii.gz', sep = '')
  net10.assignments <- cere_query_10(data_struct=all_data[[donor]], atlas=atlas, MNI_coords=all_data[[donor]]$raw_samp, rad=0, afni.dir)
  write.csv(x=net10.assignments, file=paste(atlas_dir, 'MDTB_10cere_', donor, '_10net.csv', sep='') )
}
# ----------------


## Read the Atlas Overlap Information that was just created above
# ----------------
for ( donor in donor.nums ){
  atlas_1   <- read.csv(paste(atlas_dir,'MDTB_10cere_', donor, '_10net.csv', sep=''))
  atlas_2   <- read.csv(paste(atlas_dir,'/BucknerMNI152_cere_', donor, '_7net.csv', sep=''))
  
  all_data[[donor]][["MDTB_10"]] <- atlas_1$x # store atlas assignments in the all_data structure
  all_data[[donor]][["BucknerMNI152_7"]] <- atlas_2$x # store atlas assignments in the all_data structure
}
# ----------------


## Info about split label names and color schemes
# ----------------
network_names     <- NULL
mdtb_names_in  <- read.csv(paste(base.dir, '/Reference_files/MDTB_10network_names.CSV', sep = ''), header=FALSE)
network_names[['mdtb']] <- as.character(mdtb_names_in$V1)

buckner_names_in  <- read.csv(paste(base.dir, '/Reference_files/7network_names.csv', sep = ''), header=FALSE)
network_names[['buckner7']] <- as.character(buckner_names_in$V1)
# ----------------


## Calculate and write frequency information about for each atlas and each donor, output to csv file
# ----------------
atlas_vector <- c('BucknerMNI152_7', 'MDTB_10')
ref_vector   <- c('buckner7', 'mdtb')

for (idx in 1:length(atlas_vector)){ # idx <- 1
  atlas_name <- atlas_vector[idx]
  all_tables <- NULL
  for ( donor in donor.nums ){
    cur_ref   <- network_names[[ref_vector[idx]]]
    cur_atlas <- all_data[[donor]][[atlas_name]]
    cur_names    <- cur_ref[cur_atlas+1]
    cur_table    <- table(cur_names)
    unique_ref   <- unique(cur_ref)
    add_these    <- unique_ref[!unique_ref %in% cur_names]
    cur_table[add_these] <- 0
    sort_table   <- cur_table[sort(names(cur_table))]
    all_tables   <- cbind(all_tables, sort_table)
  }
  colnames(all_tables) <- donor.nums
  out_path <- paste(atlas_dir, '/', atlas_name, '_freq.csv', sep = '')
  write.csv(x = all_tables, file = out_path)
}
save(file=paste(base.dir, '/Output/RData/all_data_cere_net10.Rdata', sep=''), x=all_data)
# load(file=paste(data_path, '/all_data_cere_10.Rdata', sep=''))
# ----------------


