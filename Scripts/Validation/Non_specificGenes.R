# -------------------------------
# -------------------------------
# Code to test Gene-FC correlation analyses using non-networ-specific genes:
# 
# Uncovering the Genetic Profiles Underlying the Intrinsic Organization of the Human Cerebellum
# Yaping Wang, Lin Chai, Deying Li, Chaohong Gao, Congying Chu, Zhengyi Yang, 
# Yu Zhang, Junhai Xu, Jens Randel Nyengaard1, Bing Liu, Kristoffer Hougaard Madsen, Tianzi Jiang, Lingzhong Fan
#
#
# Written by: Yaping Wang
# Contact:    wangyaping19@mails.ucas.ac.cn
# Noted: this code is need to run after 02_Gene-FC_Correlation
# -------------------------------
# -------------------------------

# Test non network-specific genes
# -----------------------------------
writ <- NULL
for (x in 8:11) { # x=2
  print(x)
  Sig <- NULL
  R_random <- NULL
  for (i in 1:10000){ # i=1
    print(i)
    gene.idxs <- which(rownames(all_data[['10021']]$cere_expr) %in% unique(cere.n6.sig.genes))
    non_id <- which(!(rownames(all_data[['10021']]$cere_expr) %in% unique(cere.n6.sig.genes)))
    non <- as.character(rownames(all_data[['10021']]$cere_expr[ non_id,]))
    test <- as.character(rownames(all_data[['10021']]$cere_expr)[(rowMeans(all_data$`10021`$cere_expr_nonorm)>x)]) # n =1821
    # test <- as.character(rownames(all_data[['10021']]$cere_expr)[(rowMeans(all_data$`10021`$cere_expr_nonorm)>12)]) # n =1821
    
    # intersect(non, cere.n6.sig.genes)  #65
    cur.nm <- intersect(non,test)
    # non.specific <- which(rownames(all_data[['10021']]$cere_expr) %in% non)  
    non.specific <- which(rownames(all_data[['10021']]$cere_expr) %in% cur.nm) 
    non.specific.idxs <- sample(non.specific, 443, replace = F)
    out.mats.17  <- NULL
    p_correct <- NULL
    p <- NULL
    # Calculate each gene co-experssion matrix seperately for each subject
    for ( donor in donor.nums[1:2] ){
      # donor <- donor.nums[[1]]
      # cur.mat           <- all_data[[donor]]$cere_expr_17[non.specific, ] # the difference between cort_expr_17 and cere_expr_17
      cur.mat           <- all_data[[donor]]$cere_expr_17[non.specific.idxs, ] # the difference between cort_expr_17 and cere_expr_17
      untransformed     <- cor(cur.mat, method = 'spearman')
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
    test <- (all_data$`9861`$cere_expr_17[gene.idxs, ]
             + all_data$`10021`$cere_expr_17[gene.idxs, ])/2
    avg.mat   <- (out.mats.17[[1]] + out.mats.17[[2]])/2                            # average the z-transformed correlations of each donors
    rownames(avg.mat) <- plot.order.17[twohemi.net17.regions]
    colnames(avg.mat) <- plot.order.17[twohemi.net17.regions]
    
    cere2cere.fcmri.17  <- read.csv(paste0(base.dir,'/Data/FC/cere17-cor17-all_fisher/gsr_cerecere.csv'), header = FALSE)
    colnames(cere2cere.fcmri.17) <-  c("VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB", "SalVentAttnA", "SalVentAttnB", "Limbic_Pole", "Limbic_OFC", "ContC", "ContA", "ContB", "TempPar", "DefaultC", "DefaultA", "DefaultB")
    rownames(cere2cere.fcmri.17) <-  c("VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB", "SalVentAttnA", "SalVentAttnB", "Limbic_Pole", "Limbic_OFC", "ContC", "ContA", "ContB", "TempPar", "DefaultC", "DefaultA", "DefaultB")
    use_idxs        <- colnames(cere2cere.fcmri.17) %in% rownames(avg.mat)  # avg.mat is the mRNA coorelation martix of two donors which both have left and right hemisphere 49*49
    plot.func       <- cere2cere.fcmri.17[use_idxs, use_idxs]
    
    func.1d.arr <- plot.func[upper.tri(plot.func)]
    mrna.1d.arr <- avg.mat[upper.tri(avg.mat)]
    test <- cor.test(func.1d.arr, mrna.1d.arr)
    # R_cur <- test$p.value
    # R_random <- cbind(R_random, test$estimate)
    R_random <- cbind(R_random, test$p.value)
    # print(test$estimate)
    print(test$p.value)
  }
  # length R_random[,R_random<0.05]
  # wt <- NULL
  pg <- length(R_random[,R_random<0.05])/length(R_random)
  writ <- rbind(writ, cbind(x, length(cur.nm), min(R_random),max(R_random),rowMeans(R_random), rowMedians(R_random), pg ))
  colnames(writ) <- c("Threshold", "GeneCount", "Min", "Max", "Means", "Medians", "Percentage")
  print(writ)
}
# -----------------------------------

# Non threshold
# -----------------------------------
Sig <- NULL
R_random <- NULL
for (i in 1:10000){ # i=1
  print(i)
  gene.idxs <- which(rownames(all_data[['10021']]$cere_expr) %in% unique(cere.n6.sig.genes))
  non_id <- which(!(rownames(all_data[['10021']]$cere_expr) %in% unique(cere.n6.sig.genes)))
  non <- as.character(rownames(all_data[['10021']]$cere_expr[ non_id,]))
  # test <- as.character(rownames(all_data[['10021']]$cere_expr)[(rowMeans(all_data$`10021`$cere_expr_nonorm)>x)]) # n =1821
  # intersect(non, cere.n6.sig.genes)  #65
  # cur.nm <- intersect(non,test)
  non.specific <- which(rownames(all_data[['10021']]$cere_expr) %in% non)
  # non.specific <- which(rownames(all_data[['10021']]$cere_expr) %in% cur.nm) 
  non.specific.idxs <- sample(non.specific, 443, replace = F)
  out.mats.17  <- NULL
  p_correct <- NULL
  p <- NULL
  # Calculate each gene co-experssion matrix seperately for each subject
  for ( donor in donor.nums[1:2] ){
    # donor <- donor.nums[[1]]
    # cur.mat           <- all_data[[donor]]$cere_expr_17[non.specific, ] # the difference between cort_expr_17 and cere_expr_17
    cur.mat           <- all_data[[donor]]$cere_expr_17[non.specific.idxs, ] # the difference between cort_expr_17 and cere_expr_17
    untransformed     <- cor(cur.mat, method = 'spearman')
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
  test <- (all_data$`9861`$cere_expr_17[gene.idxs, ]
           + all_data$`10021`$cere_expr_17[gene.idxs, ])/2
  avg.mat   <- (out.mats.17[[1]] + out.mats.17[[2]])/2                            # average the z-transformed correlations of each donors
  rownames(avg.mat) <- plot.order.17[twohemi.net17.regions]
  colnames(avg.mat) <- plot.order.17[twohemi.net17.regions]
  
  cere2cere.fcmri.17  <- read.csv(paste0(base.dir,'/Data/FC/cere17-cor17-all_fisher/gsr_cerecere.csv'), header = FALSE)
  colnames(cere2cere.fcmri.17) <-  c("VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB", "SalVentAttnA", "SalVentAttnB", "Limbic_Pole", "Limbic_OFC", "ContC", "ContA", "ContB", "TempPar", "DefaultC", "DefaultA", "DefaultB")
  rownames(cere2cere.fcmri.17) <-  c("VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB", "SalVentAttnA", "SalVentAttnB", "Limbic_Pole", "Limbic_OFC", "ContC", "ContA", "ContB", "TempPar", "DefaultC", "DefaultA", "DefaultB")
  use_idxs        <- colnames(cere2cere.fcmri.17) %in% rownames(avg.mat)  # avg.mat is the mRNA coorelation martix of two donors which both have left and right hemisphere 49*49
  plot.func       <- cere2cere.fcmri.17[use_idxs, use_idxs]
  
  func.1d.arr <- plot.func[upper.tri(plot.func)]
  mrna.1d.arr <- avg.mat[upper.tri(avg.mat)]
  test <- cor.test(func.1d.arr, mrna.1d.arr)
  # R_cur <- test$p.value
  # R_random <- cbind(R_random, test$estimate)
  R_random <- cbind(R_random, test$p.value)
  # print(test$estimate)
  print(test$p.value)
}
# length R_random[,R_random<0.05]

pg <- length(R_random[,R_random<0.05])/length(R_random)
writ <- rbind(writ, cbind("None", length(non), min(R_random),max(R_random),rowMeans(R_random), rowMedians(R_random), pg ))
colnames(writ) <- c("Threshold", "GeneCount", "Min", "Max", "Means", "Medians", "Percentage")
print(writ)
# -----------------------------------
write.csv(x=writ, file=paste0(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet28_ControlTest.csv')) # same as generated in the cere script

