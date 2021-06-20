# -------------------------------
# -------------------------------
# Code to run AHBA pre-process in:
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

# R package needed
# ----------------
# use the BiocManager::install() command if they haven't been installed
library(data.table) 
library(WGCNA) 
library(limma) 
# ----------------

# Modify these filepaths for your local directory structure
# ----------------
afni.dir <- '/share/soft/linux_centos_7_64'                       # enter the abin directory for local install of AFNI
base.dir <- '/n01dat01/ypwang/AHBA/CerebellarGeneFCCorrelation'   # enter the directory containing the script repository

# Source function library
source(paste(base.dir, '/Scripts/function_library.R', sep = ''))
# ----------------

# SUBJECT LIST
# ----------------
data_path  <- paste(base.dir, '/Data/AHBA/AHBA_original_data', sep = '') # path to AHBA microarray data
donor.nums <- c('9861', '10021', '12876', '14380', '15496', '15697')
filenames  <- c('donor_10021', 'donor_9861', 
                'donor_12876', 'donor_15697',
                'donor_14380', 'donor_15496')
# ----------------

# Read/process AHBA microarray data for each donor
# ----------------
all_data <- NULL
for ( donor in donor.nums ) {
  file <- paste('donor_', donor, sep='')
  print(paste('Reading and collapsing data for donor: ', donor, sep = ''))
  
  # Initiate empty field for this donor's data
  all_data[[donor]] <- NULL
  
  
  # Read Sample Information
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
  
  
  # Read Information about each Gene Probe (~50,000 probes for ~20,000 genes)
  pname     <- paste(data_path, file, 'Probes.csv', sep='/')
  probes    <- read.csv(pname)
  all_data[[donor]]$raw_probes <- probes
  
  
  # Read Gene Ontology Info
  oname    <- paste(data_path, file, 'Ontology.csv', sep='/')
  ont_data <- read.csv(oname)
  all_data[[donor]]$raw_ont <- ont_data
  
  
  # Output filenames for collapseRows
  write_micro_name   <- paste(data_path, '/', file, '/collapsed_micro_', donor, '.csv', sep="")
  write_select_name  <- paste(data_path, '/', file, '/selectedRows_', donor, '.csv', sep="")
  write_grp2row_name <- paste(data_path, '/', file, '/grp2row_', donor,'.csv', sep="")
  
  # Discard probes without an entrez-id
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
# save(file=paste(data_path, '/all_data.Rdata', sep=''), x=all_data)
# load(file=paste(data_path, '/all_data.Rdata', sep=''))
# ----------------


## Calculate the mean normalized expression for Cerebellum and Cortex seperetedly
# ----------------
# ID numbers for Cortex structures
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
  cere           <- as.numeric(cere_in_wn$structure_id)
  cere_idxs                              <- all_data[[donor]]$raw_samp$structure_id %in% cere
  all_data[[donor]]$cere_samples         <- all_data[[donor]]$raw_samp[cere_idxs,]                    # select info from cerebellar samples
  all_data[[donor]]$all_cere_micros      <- as.data.frame(all_data[[donor]]$micro_collapsed)[, cere_idxs]
  
  # Mean Normalize expression values
  all_data[[donor]]$ceres_meanNorm     <- t(apply(all_data[[donor]]$all_cere_micros, 1, function(x) x-mean(x) ))
  all_data[[donor]]$corts_meanNorm     <- t(apply(all_data[[donor]]$all_cort_micros, 1, function(x) x-mean(x) ))
}
# ----------------



# Identify which samples overlap with Yeo/Buckner Atlases, using samples' MNI locations provided by AHBA, 
# this step will take a long time: 50 min per donor, so 50*6=300min/60min/h=5h 
# -------------------------------------------
atlas_dir <- paste(base.dir, '/Output/Atlas_overlap/', sep = '')
for ( donor in donor.nums ){
  base_path <- paste(base.dir, '/Data/', sep = '')
  path      <- paste(base_path, 'Yeo_JNeurophysiol11_SplitLabels/Yeo_JNeurophysiol11_SplitLabels/MNI152/', sep = '')
  
  # 17-network split component - Cortex
  atlas_names       <- 'Yeo2011_17Networks_N1000.split_components.FSL_MNI152_FreeSurferConformed_1mm.nii'
  atlas             <- paste(path, atlas_names, sep = '')
  net17.assignments <- cort_query(data_struct=all_data[[donor]], atlas=atlas, MNI_coords=all_data[[donor]]$raw_samp, rad=0, afni.dir)
  write.csv(x=net17.assignments, file=paste(atlas_dir, 'splitLabel_cort_', donor, '_17net.csv', sep='') )
  
  # 7-network split component - Cortex
  atlas_names       <- 'Yeo2011_7Networks_N1000.split_components.FSL_MNI152_FreeSurferConformed_1mm.nii'
  atlas             <- paste(path, atlas_names, sep = '')
  net7.assignments  <- cort_query(data_struct=all_data[[donor]], atlas=atlas, MNI_coords=all_data[[donor]]$raw_samp, 0, afni.dir)
  write.csv(x=net7.assignments, file=paste(atlas_dir, '/splitLabel_cort_', donor, '_7net.csv', sep=''))
  
  
  # 17-network - Cerebellum
  cerebellum.path <- paste(base_path, 'Buckner_JNeurophysiol11_MNI152/', sep = '')
  atlas           <- paste(cerebellum.path, 'Buckner2011_17Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz', sep = '')
  atlas_conf      <- paste(cerebellum.path, 'Buckner2011_17NetworksConfidence_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz', sep = '')
  net17.assignments.cerebellum <- cere_query(data_struct=all_data[[donor]], atlas=atlas, atlas_conf=atlas_conf, MNI_coords=all_data[[donor]]$raw_samp, afni.dir)
  write.csv(x=net17.assignments.cerebellum[[1]], file=paste(atlas_dir, '/BucknerMNI152_cere_', donor, '_17net.csv', sep=''))
  
  # 7-network - Cerebellum
  path        <- paste(base_path, 'Buckner_JNeurophysiol11_MNI152/', sep = '')
  atlas       <- paste(cerebellum.path, 'Buckner2011_7Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz', sep = '')
  atlas_conf  <- paste(cerebellum.path, 'Buckner2011_7NetworksConfidence_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz', sep = '')
  net7.assignments.cerebellum <- cere_query(data_struct=all_data[[donor]], atlas=atlas, atlas_conf=atlas_conf, MNI_coords=all_data[[donor]]$raw_samp, afni.dir)
  write.csv(x=net7.assignments.cerebellum[[1]], file=paste(atlas_dir, '/BucknerMNI152_cere_', donor, '_7net.csv', sep=''))
}
# ----------------

# Read the Atlas Overlap Information that was just created above
# -------------------------------------------
for ( donor in donor.nums ){ # donor <- "9861"
  types  <- c('splitLabel_cort_', 'BucknerMNI152_cere_')
  n_regs <- c('7','17')
  
  for (type in types){ # for each atlas # type <- 'BucknerMNI152_cere_'
    for (nreg in n_regs){ # for each atlas subtype # nreg <- '7'
      atlas_name <- paste(type, donor, '_', nreg, 'net.csv', sep = '')
      atlas_in   <- read.csv(file = paste(atlas_dir, atlas_name, sep=''))
      
      cur_n    <- strsplit(atlas_name, '_')[[1]][1]
      use_name <- paste(gsub('/', '', cur_n), nreg, sep = '_')
      all_data[[donor]][[use_name]] <- atlas_in$x # store atlas assignments in the all_data structure
    }
  }
}
# -------------------------------------------

# Info about split label names and color schemes
# ---------------------------------------
buckner_names     <- NULL
buckner_names_in  <- read.csv(paste(base.dir, '/Reference_files/7network_names.csv', sep = ''), header=FALSE)
buckner_names[['sev']] <- as.character(buckner_names_in$V1)

buckner_names_in  <- read.csv(paste(base.dir, '/Reference_files/17network_names.csv', sep = ''), header=FALSE)
buckner_names[['sevteen']] <- as.character(buckner_names_in$V1)

atlas_key           <- read.table(paste(base_path, '/Yeo_JNeurophysiol11_SplitLabels/Yeo_JNeurophysiol11_SplitLabels/MNI152/7Networks_ColorLUT_freeview.txt', sep=''))
colnames(atlas_key) <- c('ID','Name','R','G','B','A')
split.names.7        <- atlas_key$Name[2:length(atlas_key$Name)]
buckner_names[['split.sev']]  <- split.names.7

atlas_key           <- read.table(paste(base_path, '/Yeo_JNeurophysiol11_SplitLabels/Yeo_JNeurophysiol11_SplitLabels/MNI152/17Networks_ColorLUT_freeview.txt', sep=''))
colnames(atlas_key) <- c('ID','Name','R','G','B','A')
split.names.17        <- atlas_key$Name[2:length(atlas_key$Name)]
buckner_names[['split.sevteen']]  <- split.names.17
# -------------------------------------------


# Calculate and write frequency information about for each atlas and each donor, output to csv file
# -------------------------------------------
atlas_vector <- c('BucknerMNI152_7', 'BucknerMNI152_17', 'splitLabel_7', 'splitLabel_17')
ref_vector   <- c('sev', 'sevteen', 'split.sev', 'split.sevteen')

for (idx in 1:length(atlas_vector)){ # idx <- 3
  atlas_name <- atlas_vector[idx]
  all_tables <- NULL
  
  for ( donor in donor.nums ){ # donor <- "9861"
    cur_ref   <- buckner_names[[ref_vector[idx]]]
    cur_atlas <- all_data[[donor]][[atlas_name]]
    if ( length(grep('splitLabel', atlas_name)) > 0){
      category.splits <- strsplit2(x=cur_ref, split='_')
      tmp             <- paste(category.splits[,3][grep('Limbic', category.splits[,3])], category.splits[,4][grep('Limbic', category.splits[,3])], sep = '_') # limbic names are a special case
      category.splits[,3][grep('Limbic', category.splits[,3])] <- tmp
      cur_ref      <- c('NONE', category.splits[,3])
      cur_names    <- cur_ref[cur_atlas+1]
    } else {
      cur_names    <- cur_ref[cur_atlas+1]
    }
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
# -------------------------------------------
save(file=paste(data_path, '/all_data_cere_net17.Rdata', sep=''), x=all_data)
# load(file=paste(data_path, '/all_data_cere_net17.Rdata', sep=''))



