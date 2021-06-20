# Function library for cerebellar Gene-FC analyses

ahba_diff_expression <- function(all_data, use.donors, use.regions, reg2yeo, dat2avg, rest.networks){
  
  out <- averageWithinCortNetworks(all_data=all_data, 
                                   use.regions=use.regions, 
                                   reg2yeo=reg2yeo, 
                                   use.donors=use.donors, # LH donors
                                   dat2avg=dat2avg, 
                                   regs=rest.networks)
  expr    <- out[[1]]  #  all.expr = num [1:20738, 1:28], rows means the 20738 genes, colums means the 4 donor * 7 networks name
  regions <- out[[2]]  #  region_arr = donor_arr = chr[1: 28], 28 = 7 network * 4 donor, 28 is the region ID, 4 donor * 7 networks name
  donors  <- out[[3]]  #  So here we get the averaged expression value for each gene(20738) in each networks(7) in each donors(4)
  colnames(expr) <- paste(donors, regions, sep = '_')
  
  # Use limma to calculate differential expression for each network, relative to all others
  # ------------------------
  fac              <- as.factor(regions)    # the network type of each column in 'expr', get 7 levels which corresponde to 7 networks name
  design           <- model.matrix(~0 + fac) # design matrix, 0 means no intercept
  colnames(design) <- gsub('fac', '', colnames(design)) # colnames(design) = fac+networks name, eg., facCont, after this step, we get the colname euqal to the networks name
  corfit           <- duplicateCorrelation(expr, design, block=donors) 
  # Use duplicate correlation to account for random effect of subject. The alternative is to use donor as an explicit fixed effect
  # in the model for a slight boost in statistical power. Random effects approach is used here for consistency with later GTEx/Brainspan 
  # analyses, where the random effects approach allows for more data retention (i.e. keep donors without complete pairwise data)
  # corfit give us 3 outcomes: 
  #   consensus.correlation:	the average estimated inter-duplicate correlation. The average is the trimmed mean of the individual correlations on the atanh-transformed scale.
  #   cor:	same as consensus.correlation, for compatibility with earlier versions of the software. 与consensus.correlation相同，以与该软件的早期版本兼容
  #   atanh.correlations:	numeric vector of length nrow(object)/ndups giving the individual genewise atanh-transformed correlations.
  
  # Calculate Differential expression for each region
  # ------------------------
  sig.genes       <- NULL
  cort.foldchange <- NULL
  cort.foldchange[['q05_genes']] <- list()
  for ( net in rest.networks ) { # rest.networks = regioon.names, is the array of the 7 networks name
    print(paste('Getting preferential expression for: ', net, sep = ''))
    
    # negative weight for the contrast matrix, depends on number of comparison networks
    mult.term    <- round(1/(length(region.names)-1),6) # round(a,b),四舍五入 , a 是四舍五入的对象， b是保留的小数位
    o.nets       <- region.names[region.names != net]   # name of the other networks, 
    cur.contrast <- paste('1*', net, '-', mult.term, '*', paste(o.nets, collapse = paste('-', mult.term, '*', sep = '')), sep = '')
    
    # Make the contrast matrix
    cmtx <- NULL
    cmtx <- makeContrasts(contrasts=cur.contrast, levels=colnames(design)) # cmtx is our contrast matrix, eg., for vis, is1*Vis-0.166667*Default-0.166667*Cont-0.166667*Limbic-0.166667*VentAttn-0.166667*DorsAttn-0.166667*SomMot
    tmplm   <- lmFit(expr, design, block=donors, correlation=corfit$consensus.correlation ) # Fit the linear model to the data
    # If we only wanted to compared diseased to normal, we could do an ordinary two group comparison.
    # Since we need to make comparisons both within and between subjects, it is necessary to treat Patient
    # as a random effect. This can be done in limma using the duplicateCorrelation function
    fit     <- eBayes(contrasts.fit( tmplm, cmtx ) ) # contrast.fit: fit the linear model to estimate a set of contrast
    # eBtaes(fit): Given a microarray linear model fit, compute moderated t statixtics
    cort.foldchange[['fit_df']][[net]] <- fit        # cort.foldchange has three list: fit_df = fit; stats= ordered topTable(fit, number=Inf); q05_genes = row names of selected stats
    tmp     <- topTable(fit, number=Inf)             # A number of summary statistics are presented by topTable() for the top genes and the selected contrast.
    cort.foldchange[['stats']][[net]] <- tmp[order(rownames(tmp)),] 
    
    # Positive fold change, FDR corrected p<0.01
    pos.idxs       <- which(cort.foldchange$stats[[net]]$logFC > 0)         # The logFC column gives the value of the contrast. Usually this represents a log2-fold change  between two or more experimental 
    # conditions although sometimes it represents a log2-expression level.
    adjusted.ps    <- which(cort.foldchange$stats[[net]]$adj.P.Val <= .05)  # adj.P.Value is the p-value adjusted for multiple testing
    # reasons: As gene expression is relatively homogenous across the cortex, an arbitrary fold change threshold could
    # mask subtle, yet biologically meaningful, variations in expression
    # Accordingly, we adopted a statistical threshold to account for multiple comparisons 
    # (false discovery rate (FDR) corrected q ≤ 0.01) and independant data to replicate this results, assess population-level stability.
    genes.tmp      <- rownames(cort.foldchange$stats[[net]])[intersect(adjusted.ps, pos.idxs)]
    cort.foldchange[['q05_genes']][[net]] <- genes.tmp
    print(length(genes.tmp)) # the number of diff genes of each networks
    sig.genes <- c(sig.genes, genes.tmp)
  }
  return(list(cort.foldchange, sig.genes))
}




avg_parcel_expression <- function(all_data, cerebellar.atlas, cerebellar.num, cortical.atlas, cortical.num, cort.use.regions){
  
  for ( donor in names(all_data)){
    # print(paste('Averaging cort and cerebellun data for: ', donor, sep = ''))
    
    # Get the cerebellar parcels that contain samples for this subject
    # ----------------------
    cere.regions.tmp <- unique(all_data[[donor]][[cerebellar.atlas]])
    cere.regions     <- sort(cere.regions.tmp[cere.regions.tmp != 0]) # Don't count areas that weren't assigned (i.e not 0)
    
    # Mean Normalized cerebellar expression values, error in the difference between ceres_meanNorm and cere_meanNorm
    all_data[[donor]]$cere_expr <- averageCereExpr(cere.region=cere.regions, 
                                                   data.struct=all_data[[donor]], 
                                                   buckner.names=buckner.names[[ cerebellar.num]], 
                                                   samp.labels=all_data[[donor]][[cerebellar.atlas]], 
                                                   data.type='ceres_meanNorm')
    # non-mean Normalized cerebellar expression values
    all_data[[donor]]$cere_expr_nonorm <- averageCereExpr(cere.region=cere.regions, 
                                                          data.struct=all_data[[donor]], 
                                                          buckner.names=buckner.names[[cerebellar.num]], 
                                                          samp.labels=all_data[[donor]][[cerebellar.atlas]], 
                                                          data.type='all_cere_micros')
    
    # Mean Normalized cortical expression values
    all_data[[donor]]$cort_expr         <- averageCortExpr(use_regions=cort.use.regions, 
                                                           data_struct=all_data[[donor]], 
                                                           max_reg=cortical.num, 
                                                           sample_labels=all_data[[donor]][[cortical.atlas]], 
                                                           data_type='corts_meanNorm')
    # non-Mean Normalized cortical expression values
    all_data[[donor]]$cort_expr_nonorm  <- averageCortExpr(use_regions=cort.use.regions, 
                                                           data_struct=all_data[[donor]], 
                                                           max_reg=cortical.num, 
                                                           sample_labels=all_data[[donor]][[cortical.atlas]], 
                                                           data_type='all_cort_micros')
  }
  return(all_data)
}



# Cross-reference MNI coordinates with group atlas
# ----------------------------
cort_query <- function(data_struct, atlas, MNI_coords, rad, afni.dir) {
  
  # Initialize output arrays
  net_arr     <- NULL
  in_network  <- NULL
  out_network <- NULL
  # Iterate of all cortical sample locations
  # ---------------------
  for ( iter in 1:length(data_struct$cort_samples[,1]) ){
    use_coords       <- MNI_coords[MNI_coords$well_id %in% data_struct$cort_sample$well_id[iter],]
    coords           <- c(use_coords$mni_x, use_coords$mni_y, use_coords$mni_z)
    network          <- query_atlas(atlas, coords, 0, afni.dir)
    print(network)
    
    # Check neighboring voxels if needed
    if ( rad == 1 ){
      neigh_arr <- query_atlas(atlas, coords, 1, afni.dir)
      if ( length(neigh_arr[neigh_arr %in% 0]) == 27 ){
        out_network <- rbind(out_network, data_struct$cort_samples[iter,])
      } else if ( length(unique(neigh_arr[neigh_arr > 0])) == 1 ){
        if ( network == 0 ){
          use     <- neigh_arr[neigh_arr>0]
          network <- as.integer(names(sort(summary(as.factor(use))))[1])
        }
        in_network  <- rbind(in_network, data_struct$cort_samples[iter,]) 
      } else {
        network <- 999
        out_network <- rbind(out_network, data_struct$cort_samples[iter,])
      }
      net_arr   <- c(net_arr, network)
    } else {
      net_arr   <- c(net_arr, network)
    }
  }
  return(net_arr)
}



cere_query <- function(data_struct, atlas, atlas_conf, MNI_coords, afni.dir) {
  net_arr      <- NULL
  net_arr_conf <- NULL
  in_network   <- NULL
  conf_network <- NULL
  out_network  <- NULL
  
  for ( iter in 1:dim(data_struct$cere_samples)[1] ){
    print(iter)
    #coords       <- c(data_struct$striat_samples[iter,]$mni_x, data_struct$striat_samples[iter,]$mni_y, data_struct$striat_samples[iter,]$mni_z)
    use_coords <- MNI_coords[MNI_coords$well_id %in% data_struct$cere_samples$well_id[iter],]
    coords     <- c(use_coords$mni_x, use_coords$mni_y, use_coords$mni_z)
    
    network      <- query_atlas(atlas, coords, 0, afni.dir)
    network_conf <- query_atlas(atlas_conf, coords, 0, afni.dir)
    
    # Check neighboring voxels
    if ( network == 0 ){
      neigh_arr <- query_atlas(atlas, coords, 1, afni.dir)
      num_zeros <- length(neigh_arr[neigh_arr %in% 0])
      if ( num_zeros == 27 ){
        neigh_arr <- query_atlas(atlas, coords, 2, afni.dir)
        num_zeros <- length(neigh_arr[neigh_arr %in% 0])
        if ( num_zeros == 125 ){
          network <- 0 
        } else {
          use       <- neigh_arr[neigh_arr>0]
          use_table <- summary(as.factor(use))
          network   <- as.integer(names(sort(use_table)))
          network_conf <- 0
        }
      } else {
        use       <- neigh_arr[neigh_arr>0]
        use_table <- summary(as.factor(use))
        network   <- as.integer(names(sort(use_table)))
        network_conf <- 0
      }
    }
    print(paste(iter,':', network))
    net_arr      <- c(net_arr, network[1])
    net_arr_conf <- c(net_arr_conf, network_conf)
    
  }
  return(list(net_arr, net_arr_conf))
}

cere_query_10 <- function(data_struct, atlas, MNI_coords, rad, afni.dir) {
  
  # Initialize output arrays
  net_arr     <- NULL
  in_network  <- NULL
  out_network <- NULL
  # Iterate of all cere sample locations
  # ---------------------
  for ( iter in 1:length(data_struct$cere_samples[,1]) ){
    use_coords       <- MNI_coords[MNI_coords$well_id %in% data_struct$cere_sample$well_id[iter],]
    coords           <- c(use_coords$mni_x, use_coords$mni_y, use_coords$mni_z)
    network          <- query_atlas(atlas, coords, 0, afni.dir)
    print(network)
    
    # Check neighboring voxels if needed
    if ( rad == 1 ){
      neigh_arr <- query_atlas(atlas, coords, 1, afni.dir)
      if ( length(neigh_arr[neigh_arr %in% 0]) == 27 ){
        out_network <- rbind(out_network, data_struct$cere_samples[iter,])
      } else if ( length(unique(neigh_arr[neigh_arr > 0])) == 1 ){
        if ( network == 0 ){
          use     <- neigh_arr[neigh_arr>0]
          network <- as.integer(names(sort(summary(as.factor(use))))[1])
        }
        in_network  <- rbind(in_network, data_struct$cere_samples[iter,]) 
      } else {
        network <- 999
        out_network <- rbind(out_network, data_struct$cere_samples[iter,])
      }
      net_arr   <- c(net_arr, network)
    } else {
      net_arr   <- c(net_arr, network)
    }
  }
  return(net_arr)
}


readAtlasOverlap <- function(all_data, filenames, atlas.dir){
  # Read in the atlas overlap information that aligns the MNI coordinates of each sample to the 
  # Buckner/Yeo striatal and cortical atlases. These were calculated ahead of time in 'fyi_preprocess_data.R'
  # 
  # Args:
  #   all_data:  data structure containing microarray expr, as well as sample/ontology/probe info
  #   filenames: donor samples, corresponds to field names in 'all_data'
  #   atlas.dir: directory with all the atlas information
  #
  # Returns:
  #   Atlas overlap information for each subject
  
  
  for ( donor in filenames ){
    types  <- c('splitLabel_cort_', 'BucknerMNI152_cere_')
    n.regs <- c('7','17')
    
    for (type in types){
      for (nreg in n.regs){
        atlas.name <- paste0('/', type, donor, '_', nreg, 'net.csv')
        atlas.in   <- read.csv(file = paste(atlas.dir, atlas.name, sep=''))
        atlas.out  <- atlas.in$x
        cur.n      <- strsplit(atlas.name, '_')[[1]][1]
        use_name   <- paste(gsub('/', '', cur.n), nreg, sep = '_')
        print(use_name)
        all_data[[donor]][[use_name]] <- atlas.out
      }
    }
  }
  return(all_data)
}


get_region_info_cort <- function(all_data, filenames, reg_IDs, name){
  # Given a list of AHBA region ID numbers, subset the data for each subject
  
  # Args:
  #   all_data:  data structure containing microarray expression, as well as sample/ontology/probe info, 
  #   filenames: list of subject numbers, which correspond to field names of 'all_data'
  #   reg_IDs:   AHBA ontology IDs for the region being analyzed
  #   name:      string for distinguishing the output data
  #
  # Returns:
  #   Region specific versions of microarray/ ontology/ probe/ and sample information
  
  
  for ( donor in filenames ){
    print(paste('Getting ', name, ' data for: ', donor, sep = ''))
    
    # Get expression and sample information about cortical samples
    # -----------------
    idxs      <- all_data[[donor]]$raw_samp$structure_id %in% reg_IDs
    
    samp_name <- paste(name, '_samples', sep='')
    all_data[[donor]][[samp_name]]  <- all_data[[donor]]$raw_samp[idxs, ] # sample information
    
    acro_name <- paste(name, '_acros', sep='')
    all_data[[donor]][[acro_name]]  <- factor(all_data[[donor]][[samp_name]]$structure_acronym) # list of corresponding struct acronyms
    
    micro_name <- paste('all_', name, '_micros', sep='')
    all_data[[donor]][[micro_name]] <- as.data.frame(all_data[[donor]]$micro_collapsed)[, idxs]
    
    mean_name <- paste(name, '_meanNorm', sep='')
    all_data[[donor]][[mean_name]]  <- t(apply(all_data[[donor]][[micro_name]], 1, function(x) x-mean(x) ))
    
    pa_name <- paste('all_', name, '_pas', sep='')
    all_data[[donor]][[pa_name]] <- as.data.frame(all_data[[donor]]$pas_collapse)[, idxs]
    
  }
  return (all_data)
}


get_region_info_cere <- function(all_data, filenames, name){
  # Given a list of AHBA region ID numbers, subset the data for each subject
  
  # Args:
  #   all_data:  data structure containing microarray expression, as well as sample/ontology/probe info, 
  #   filenames: list of subject numbers, which correspond to field names of 'all_data'
  #   reg_IDs:   AHBA ontology IDs for the region being analyzed
  #   name:      string for distinguishing the output data
  #
  # Returns:
  #   Region specific versions of microarray/ ontology/ probe/ and sample information
  for ( donor in donor.nums ){
    print(paste('Retrieving cerebellum information for: ', donor, sep = ''))
    # Cere - ID numbers for cerebellum structures
    # -----------------
    test           <- all_data[[donor]]$raw_samp
    cere_in        <- test[which(test[,"slab_type"]=="CB"),] 
    cere_in_wn     <- cere_in[!grepl("nucleus", cere_in$structure_name),]
    cere           <- as.numeric(cere_in_wn$structure_id)
    reg_IDs        <- cere
    
    print(paste('Getting ', name, ' data for: ', donor, sep = ''))
    
    # Get expression and sample information about cortical samples
    # -----------------
    idxs      <- all_data[[donor]]$raw_samp$structure_id %in% reg_IDs
    
    samp_name <- paste(name, '_samples', sep='')
    all_data[[donor]][[samp_name]]  <- all_data[[donor]]$raw_samp[idxs, ] # sample information
    
    acro_name <- paste(name, '_acros', sep='')
    all_data[[donor]][[acro_name]]  <- factor(all_data[[donor]][[samp_name]]$structure_acronym) # list of corresponding struct acronyms
    
    micro_name <- paste('all_', name, '_micros', sep='')
    all_data[[donor]][[micro_name]] <- as.data.frame(all_data[[donor]]$micro_collapsed)[, idxs]
    
    mean_name <- paste(name, '_meanNorm', sep='')
    all_data[[donor]][[mean_name]]  <- t(apply(all_data[[donor]][[micro_name]], 1, function(x) x-mean(x) ))
    
    pa_name <- paste('all_', name, '_pas', sep='')
    all_data[[donor]][[pa_name]] <- as.data.frame(all_data[[donor]]$pas_collapse)[, idxs]
    
  }
  return (all_data)
}


averageWithinCortNetworks <- function(all_data, use.regions, reg2yeo, use.donors, dat2avg, regs){
  # For each subject, get average cortical expression across all cortical regions that fall 
  #  within the same Yeo functional network
  # --------------
  # Args:
  #   all_data:    overall data struct
  #   use.regions: pre-defined array of indices for regions that meet some criteria, e.g. samples in at least 2 subs
  #   reg2yeo:     maps between region id and yeo network label
  #   use.donors:  donors you want to analyse
  #   dat2avg:     name of the field containing the data we will be averaging
  #
  # Returns:
  #     Matrix that is "# of genes" X "(# of subjects*# of regions)"
  
  ct <- 0 
  all.expr   <- NULL
  region_arr <- NULL
  donor_arr  <- NULL
  for (donor in use.donors){
    cur.cort.data <- all_data[[donor]][[dat2avg]]
    
    donor_data <- NULL
    for (reg in regs) {
      idxs_for_cur_reg   <- intersect(grep(reg, reg2yeo), use.regions)
      # use.regions <- lh.use.regions.7 : pre-defined array of indices for regions that meet some criteria, e.g. samples in 4 left-hemisphere donor,
      #                                   21 PARCEL ID of left hemisphere donor 
      # reg2yeo.27:  27 parcel ID for each networks within 7 networks only for left hemisphere, maps between region id and yeo network label
      ct <- ct + length(idxs_for_cur_reg)
      if (is.null(dim(cur.cort.data[,idxs_for_cur_reg])) == TRUE){ # if this subject has only one region in the current network
        avgdat_for_cur_reg <- cur.cort.data[,idxs_for_cur_reg]
      } else { # more than one region
        avgdat_for_cur_reg <- rowMeans(cur.cort.data[,idxs_for_cur_reg], na.rm = TRUE)
      }
      if ( is.nan(avgdat_for_cur_reg[1]) == FALSE ){               # NaN：表示非数值，是“Not a Number”的缩写
        donor_data <- cbind(donor_data, avgdat_for_cur_reg)
        region_arr <- c(region_arr, reg)
        donor_arr  <- c(donor_arr, donor)
      }
    }
    all.expr <- cbind(all.expr, donor_data)
  }
  out <- list(all.expr, region_arr, donor_arr) # region_arr = donor_arr = chr[1: 28], all.expr = num [1:20738, 1:28]
  return(out)
}


averageWithinCereNetworks <- function(all_data, donor.nums, use.cere.networks, type, net_names, atlas_field){
  # For each subject, get average cerebellar expression across all cortical regions that fall within the same Buckner functional network
  # --------------
  # Args:
  #   all_data:             overall data struct
  #   use.cere.networks:    pre-defined array of indices for regions that meet some criteria, e.g. samples in at least 2 subs
  #   use.donors:           donors you want to analyse
  #   dat2avg:              name of the field containing the data we will be averaging
  #
  # Returns:
  #     Matrix that is "# of genes" X "(# of subjects*# of regions)"
  use.idxs    <- which(net_names %in% use.cere.networks)-1
  cere.idxs   <- which(net_names %in% use.cere.networks)-1
  
  # Get Average cerebellar expression for each subject, in 6 
  # ------------------------------------------------------
  out.mat     <- NULL
  regions     <- NULL
  donor.array <- NULL
  for ( donor in donor.nums ) {
    cur.cere       <- all_data[[donor]][[type]]
    cur.cere.atlas <- all_data[[donor]][[atlas_field]]
    for (idx in cere.idxs){
      cdat <- cur.cere[,cur.cere.atlas == idx]
      if (is.null(dim(cdat)) == TRUE) {
        use.arr <- cdat
      } else {
        use.arr <- rowMeans(cdat)
      }
      if ( is.nan(sum(use.arr)) == FALSE ){
        out.mat     <- cbind(out.mat, use.arr)
        regions     <- c(regions, net_names[idx+1])
        donor.array <- c(donor.array, donor)
      }
    }
  }
  out <- list(out.mat,  regions, donor.array )
  return(out)
}


getCortRegions <- function(all_data, filenames, atlas.num, thresh){
  all_use <- NULL
  for ( donor in filenames ){
    split.name        <- paste('splitLabel_', atlas.num, sep = '') # Each individual cortical parcel (values = 1-114), eg., 466 samples for donor 9861
    split.regions     <- unique(all_data[[donor]][[split.name]])   # split_region index (corresponding to the functional atlas) for each sample
    use.split.regions <- NULL
    all_use <- c(all_use, unique(split.regions[split.regions != 0]))
  }
  # regions that are present in at least X number of subjects (set with thresh variable)
  reg_counts   <- sort(table(all_use))
  regions_min  <- reg_counts[reg_counts >= thresh]
  regions.out  <- as.numeric(rownames(as.matrix(regions_min)))
  regions.out  <- sort(regions.out)
  return(regions.out)
}

getCereRegions <- function(all_data, filenames, atlas.num, thresh){
  all_use <- NULL
  for ( donor in filenames ){
    split.name        <- paste('BucknerMNI152_', atlas.num, sep = '') # Each individual cortical parcel (values = 1-114), eg., 466 samples for donor 9861
    split.regions     <- unique(all_data[[donor]][[split.name]])   # split_region index (corresponding to the functional atlas) for each sample
    use.split.regions <- NULL
    all_use <- c(all_use, unique(split.regions[split.regions != 0]))
  }
  # regions that are present in at least X number of subjects (set with thresh variable)
  reg_counts   <- sort(table(all_use))
  regions_min  <- reg_counts[reg_counts >= thresh]
  regions.out  <- as.numeric(rownames(as.matrix(regions_min)))
  regions.out  <- sort(regions.out)
  return(regions.out)
}



query_atlas <- function(atlas, coords, rad, afni.dir) {
  cmd       <- paste(afni.dir, '/3dmaskdump -nbox ', coords[1]+rad, ':', coords[1]-rad, ' ', coords[2]+rad, ':', coords[2]-rad, ' ', coords[3]+rad, ':', coords[3]-rad, ' ', atlas, sep='')
  output    <- system(cmd, intern=TRUE)
  # output    <- shell(cmd, intern=TRUE)
  neigh_arr <- NULL
  for ( idx in 1:length(output) ) {
    cur_out   <- output[[idx]]
    split_out <- strsplit(cur_out, split=' ')
    network   <- as.numeric(split_out[[1]][4])
    neigh_arr <- c(neigh_arr, network)
  }
  return(neigh_arr)
}



averageCortExpr <- function(use_regions, data_struct, max_reg, sample_labels, data_type){
  # 
  # Args:
  #   use_regions: array w/ info about the network assignment of each cortical sample
  #   data.struct:    data for the current subject 
  #   max_reg:        sometimes you just want to run on the left-hemisphere, which is equivalent to atlas indices < 57
  #   sample_labels:    numeric network assignment of each sample
  #   data.type:      the name of the cortical field in 'data.struct' that we want to average
  #
  # Returns:
  #   subject-wise cortical data averaged for each network
  n_genes  <- dim(data_struct[[data_type]])[1]
  cort_arr <- matrix(NA, n_genes, max_reg)
  
  for ( reg in use_regions ){
    cur_dat <- data_struct[[data_type]][, sample_labels == reg ] # all_data[["9861"]][["cort_meanNorm"]] rows means the 20738 genes, colums means the 466 samples lable, the label corresponds to the parcel ID
    if ( is.null(dim(cur_dat)) == TRUE ) {
      use_dat <- cur_dat
    } else {
      use_dat <- rowMeans(cur_dat)
    }
    cort_arr[,reg] <- use_dat
  }
  out           <- cort_arr
  rownames(out) <- rownames(data_struct[[data_type]])   # 20738 genes * 114 parcel ID, some colume of them is NA, because for each donors, don't have all 114 percel ID 
  return(out)
}


averageCereExpr <- function(cere.regions, data.struct, buckner.names, samp.labels, data.type) {
  # 
  # Args:
  #   cere_regions:   array w/ info about the network assignment of each striatal sample
  #   data.struct:    data for the current subject 
  #   buckner.names:  list of all possible network assignments
  #   samp.labels:    numeric network assignment of each sample
  #   data.type:      the name of the cerebellar field in 'data.struct' that we want to average
  #
  # Returns:
  #   subject-wise cerebellar data averaged for each network
  all.expr <- NULL
  cols.out <- NULL
  for (reg in cere.regions) {
    # name array to add as column headers at the finish
    col_idx  <- samp.labels == reg                     # indices for this dth subregion
    reg_data <- data.struct[[data.type]][, col_idx]
    cname    <- buckner.names[as.integer(reg) + 1]     # so here the buckner.names can contains the 1 None
    
    # if there is only one column for this iteration
    if ( is.null(dim(reg_data)) == TRUE ){
      cur_reg <- reg_data
    } else{
      cur_reg <- rowMeans(reg_data) # Get the average across all samples
    }
    # append the data to an output array
    cols.out <- c(cols.out, cname)
    all.expr <- cbind(all.expr, cur_reg)
  }
  colnames(all.expr) <- cols.out
  cere_expr          <- as.data.frame(all.expr)
  return(cere_expr)
}


