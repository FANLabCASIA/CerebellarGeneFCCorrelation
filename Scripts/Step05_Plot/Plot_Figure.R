# -------------------------------
# -------------------------------
# Code to plot in:
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


## R packages need
# ------------
library(cowplot)
library(ggpubr)
library(pheatmap)
library(corrplot)
library(cairo)
library(plotly)
library(ggplot2)
library(fmsb)
library(reshape2)
library(ggplot2)
library(ggpmisc)
library(RColorBrewer)
library(ggpubr)
library(ggsignif)
library(ggsci)
library(scales)
library(basicTrendline)
library(Rmisc)
library(stringr)
base.dir <- ''  # Modify these filepaths for your local directory structure
# ------------


## Validation: Scatterplot
# ------------
# Read in gene correlation and FC within cere
mrna_net7 <- read.csv(paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet13_Gene_cere_Corr_net7.csv', sep=''),
                      header = TRUE, row.names = 1 )
fc_net7    <- read.csv(paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet14_FC_cere_Corr_net7.csv', sep=''),
                       header = TRUE, row.names = 1)

mrna_net17 <- read.csv(paste(base.dir, "/Output/SupplementaryData/SUPPDATA_sheet11_Gene_cere_Corr_net17.csv",sep = ""), 
                       header = TRUE, row.names = 1 )
fc_net17 <- read.csv(paste(base.dir, "/Output/SupplementaryData/SUPPDATA_sheet12_FC_cere_Corr_net17.csv",sep = ""), 
                     header = TRUE, row.names = 1 )

mrna_net10 <- read.csv(paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet15_Gene_cere_Corr_net10.csv', sep=''),
                       header = TRUE, row.names = 1 )
fc_net10    <- read.csv(paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet16_FC_cere_Corr_net10.csv', sep=''),
                        header = TRUE, row.names = 1)

data <- as.data.frame((rbind(cbind(mrna_net7[upper.tri(mrna_net7)], fc_net7[upper.tri(fc_net7)]),
                             cbind(mrna_net17[upper.tri(mrna_net17)], fc_net17[upper.tri(fc_net17)]),
                             cbind(mrna_net10[upper.tri(mrna_net10)], fc_net10[upper.tri(fc_net10)]))))
data$Atlas <- c(rep("1", 15), rep("2", 45),  rep("4", 45))
names(data) <- c("Gene_correlation", "FC_correlation", "Atlas")

pdf(file = paste(base.dir, '/Output/Figure/Fig2d_GeneFC_Scatterplot.pdf', sep = ''),
    width = 10, height = 5)
par(mar = c(2, 2, 2, 2))
ggplot(data, aes(x= Gene_correlation , y= FC_correlation, color = Atlas, shape = Atlas, color_palette())) +
  geom_point(aes(color = Atlas), size = 2.0) +    # Use hollow circles, shape = 1, , shape = Atlas , ggpubr::show_point_shapes()   
  scale_color_manual(name = "" ,
                     values = c( "#f4977a","#bc1d2c" , "#779af7"), 
                     labels=c("7-network", "17-network",  "MDTB")) + #label
  scale_shape_manual(name = "",
                     values = c(15,17,16),
                     labels=c("7-network", "17-network", "MDTB")) +
  # labs(color  = "Guide name", shape = "Guide name")
  scale_fill_discrete(guide = FALSE) +
  geom_smooth(aes(color = Atlas, fill = Atlas), se =FALSE, # se =FALSE, # fill = Atlas, 
              method=lm, formula = y ~ x,   fullrange=TRUE) +
  guides()+
  stat_cor(method = "pearson",label.x = -0.45, show.legend = FALSE) +
  coord_quickmap() + 
  scale_y_continuous('Functional connectivity') + scale_x_continuous('Gene correlation') + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 18,  face ="plain"),
        axis.title.y = element_text(size = 18,  face ="plain"),
        axis.text.x  = element_text(size = 15,  face ="plain"),
        axis.text.y  = element_text(size = 15,  face ="plain"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
        plot.title = element_text(size = 20,  face ="bold", hjust = 0.5),
        axis.line.x  = element_line(size = 0.8),
        axis.line.y  = element_line(size = 0.8),
        legend.background = element_blank(), legend.key = element_blank()
  ) 

dev.off()
# ------------


## GO-plot
# ------------
# Read in the ToppGene result for GCI+ and GCI-
Enrich <- NULL
Enrich$GO$In <- read.csv(paste(base.dir, "/Reference_files/GO/Plot_GO_GCI+_fisher.CSV", sep = ""), header = TRUE)
Enrich$GO$De <- read.csv(paste(base.dir, "/Reference_files/GO/Plot_GO_GCI-_fisher.CSV", sep = ""), header = TRUE)

Enrich$Dis$In <- read.csv(paste(base.dir, "/Reference_files/GO/Plot_Dis_GCI+_fisher.CSV", sep = ""), header = TRUE)
Enrich$Dis$De <- read.csv(paste(base.dir, "/Reference_files/GO/Plot_Dis_GCI-_fisher.CSV", sep = ""), header = TRUE)


# type <- c("GO", "Path", "Dis")
type <- c("GO", "Dis")
gci <- c("In", "De")

for ( typ in type) {
  for (gc in gci) {
    Enrich[[typ]][[gc]] <- as.data.frame(Enrich[[typ]][[gc]])
  }
}


# GO
pdf(file = paste(base.dir, '/Output/Figure/Fig5a_GO_GCI+.pdf', sep = ''),
    width = 7.5, height = 10)
ggplot(Enrich$GO$In[1:20,], aes(x= GeneRatio, y = reorder(Term, GeneRatio)))+
  geom_point(aes(size = Count, color = Pvalue)) +
  scale_size_area(max_size = 8) +
  scale_color_gradient("Adjusted P", low='#a11729',high='#f4c5af') +
  # scale_color_gradient("Adjusted P", low='#b81029',high='#f4c5af') +
  theme_bw() +
  scale_y_discrete(labels = function(Term) str_wrap(Term, width = 30) ) + 
  scale_x_continuous(limits = c(0, 0.082),breaks = seq(0.00,0.08,0.02), expand = c(0,0)) +
  ggtitle("GCI+") +
  guides(color = guide_colorbar(order = 0), size = guide_legend(order = 1)) +
  theme(axis.title.x = element_text(size = 13,  face ="plain"),
        axis.title.y = element_blank(),
        axis.text.x  = element_text(size = 12,  face ="plain"),
        axis.text.y  = element_text(size = 13,  face ="plain"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
        plot.title = element_text(size = 15,  face ="plain", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(), axis.line = element_line(size = 0.8)
        # axis.ticks.x = element_blank()
        # legend.position = c(0.8,0.2), legend.text = element_text( size = 10,  face ="plain"),
        # legend.title = element_text( size = 12,  face ="plain")
  ) 
dev.off()


pdf(file = paste(base.dir, '/Output/Figure/Fig5b_GO_GCI-.pdf', sep = ''),
    width = 7.5, height = 10)
ggplot(Enrich$GO$De[1:20,], aes(x= GeneRatio, y = reorder(Term, GeneRatio)))+
  geom_point(aes(size = Count, color = Pvalue)) +
  scale_size_area(max_size = 8) +
  # scale_color_gradient("Adjusted P", low='#3b4ec1',high='#89acfa') +
  scale_color_gradient("Adjusted P", low='#18407b',high='#89acfa') + 
  # scale_color_gradient("Adjusted P", low='#E43535',high='#00438E') +
  theme_bw() +
  scale_y_discrete(labels = function(Term) str_wrap(Term, width = 30) ) + 
  scale_x_continuous(limits = c(0, 0.15), breaks = seq(0,0.15, 0.05), expand = c(0,0)) +
  ggtitle("GCI-") +
  guides(color = guide_colorbar(order = 0), size = guide_legend(order = 1)) +
  theme(axis.title.x = element_text(size = 13,  face ="plain"),
        axis.title.y = element_blank(),
        axis.text.x  = element_text(size = 12,  face ="plain"),
        axis.text.y  = element_text(size = 13,  face ="plain"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
        plot.title = element_text(size = 15,  face ="plain", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(), axis.line = element_line(size = 0.8)
        # legend.position = c(0.8,0.2), legend.text = element_text( size = 10,  face ="plain"),
        # legend.title = element_text( size = 12,  face ="plain")
  )
dev.off()



# Disease
pdf(file = paste(base.dir, '/Output/Figure/Fig5c_Dis_GCI+.pdf', sep = ''),
    width = 8, height = 7)
ggplot(Enrich$Dis$In, aes(x= GeneRatio, y = reorder(Term, GeneRatio)))+
  geom_bar(aes( fill = Pvalue), stat="identity") +
  # scale_fill_manual(values = c("#053061", "#2166AC", "#4393C3", "#92C5DE",
  #                                                          "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
  #                                                          "#D6604D", "#B2182B", "#67001F")) +
  # scale_fill_gradient("Adjusted P", low='#a11729',high='#f4c5af') +
  scale_fill_gradient("Adjusted P", low='#cf3840',high='#f4c5af') +
  theme_bw() +
  scale_y_discrete(labels = function(Term) str_wrap(Term, width = 36) ) + 
  # scale_x_continuous(expand = c(0,0))
  
  ggtitle("GCI+") +
  theme(axis.title.x = element_text(size = 14,  face ="plain"),
        axis.title.y = element_blank(),
        axis.text.x  = element_text(size = 13,  face ="plain"),
        axis.text.y  = element_text(size = 14,  face ="plain"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
        plot.title = element_text(size = 17,  face ="plain", hjust = 0.2),
        legend.text = element_text( size = 12,  face ="plain"), # legend.position = c(0.8,0.2), 
        legend.title = element_text( size = 12,  face ="plain"),
        panel.grid = element_blank(),
        panel.border = element_blank(), axis.line = element_line(size = 0.8)
  )
dev.off()

pdf(file = paste(base.dir, '//Output/Figure/Fig5d_Dis_GCI-.pdf', sep = ''),
    width = 8, height = 7)
ggplot(Enrich$Dis$De, aes(x= GeneRatio, y = reorder(Term, GeneRatio)))+
  geom_bar(aes( fill = Pvalue), stat="identity") +
  # scale_fill_manual(values = c("#053061", "#2166AC", "#4393C3", "#92C5DE",
  #                                                          "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
  #                                                          "#D6604D", "#B2182B", "#67001F")) +
  scale_fill_gradient("Adjusted P", low='#18407b',high='#89acfa')+ 
  # scale_fill_gradient("Adjusted P", low='#3b4ec1',high='#89acfa') +
  theme_bw() +
  scale_y_discrete(labels = function(Term) str_wrap(Term, width = 38) ) + 
  # scale_x_continuous(expand = c(0,0))
  # scale_x_continuous(limits = c(0, 0.12) ,breaks = c(0, 0.04, 0.08,  0.12)) +
  ggtitle("GCI-") +
  theme(axis.title.x = element_text(size = 14,  face ="plain"),
        axis.title.y = element_blank(),
        axis.text.x  = element_text(size = 13,  face ="plain"),
        axis.text.y  = element_text(size = 14,  face ="plain"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
        plot.title = element_text(size = 17,  face ="plain", hjust = 0.2),
        legend.text = element_text( size = 12,  face ="plain"), # legend.position = c(0.8,0.2), 
        legend.title = element_text( size = 12,  face ="plain"),
        panel.grid = element_blank(),
        panel.border = element_blank(), axis.line = element_line(size = 0.8)
  )
dev.off()

# ------------


## BrainSpan
# ------------
# Read in the results p value from CSEA website
all_p <- read.csv(paste(base.dir, "/Reference_files/SpatialTemporal/allp_trans.CSV", sep=""), header = TRUE)
x <- c(all_p$In0.05, all_p$In0.01, all_p$In0.001, all_p$In0.0001,
       all_p$De0.05, all_p$De0.01, all_p$De0.001, all_p$De0.0001)
# BH-FDR
p_adj <- p.adjust(x , method = "BH", n = length(x))

q_mat <- matrix(p_adj,60,8)
colnames(q_mat) <- c(  "In0.05",  "In0.01",  "In0.001",  "In0.0001",
                       "De0.05",  'De0.01',  "De0.001",  "De0.0001")
rownames(q_mat) <- all_p$Brain_Regions_Development
write.csv(q_mat, file = paste(base.dir, '/Reference_files/SpatialTemporal/all_q_trans.csv', sep = ''))
# Select the cerebellar results manually
Q <- NULL
Q$All <- read.csv(paste(base.dir, "/Reference_files/SpatialTemporal/GCI_All_adj_p_fisher.CSV", sep=""), header = TRUE, row.names = 1)
Q$In <- read.csv(paste(base.dir, "/Reference_files/SpatialTemporal/GCI_In_adj_p_fisher.CSV", sep=""), header = TRUE, row.names = 1)
Q$De <- read.csv(paste(base.dir, "/Reference_files/SpatialTemporal/GCI_De_adj_p_fisher.CSV", sep=""), header = TRUE, row.names = 1)
name <- c("All", "In", "De")
FDR.q <-NULL

for (reg in c( "In", "De") ) {
  FDR.q[[reg]] <- data.frame(period = c(rep(rownames(Q$All),4)),
                             fdr=c(-log10(Q[[reg]][,1]), -log10(Q[[reg]][,2]), 
                                   -log10(Q[[reg]][,3]), -log10(Q[[reg]][,4])), 
                             type = c(rep("1", 10), rep("2", 10),rep("3", 10), rep("4", 10)))
  FDR.q[[reg]][FDR.q[[reg]]$fdr>2.5,]$fdr <- 2.5
}

FDR.q$All <- data.frame(period = c(rep(rownames(Q$All),8)),
                        fdr=c(Q$All[,1],Q$All[,2], Q$All[,3], Q$All[,4], Q$All[,5], Q$All[,6], Q$All[,7], Q$All[,8]), 
                        type = c(rep("1", 10), rep("2", 10),rep("3", 10), rep("4", 10)
                                 , rep("1", 10), rep("2", 10), rep("3", 10), rep("4", 10)),
                        GCI = c(rep("1",40), rep("2", 40)))

mytheme = theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                axis.title.x = element_text(size = 18,  face ="plain"),
                axis.title.y = element_text(size = 18,  face ="plain"),
                axis.text.x  = element_text(size = 15,  face ="plain", angle = 45, hjust = 1),
                axis.text.y  = element_text(size = 15,  face ="plain"),
                plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
                plot.title = element_text(size = 20,  face ="bold", hjust = 0.5),
                axis.line.x  = element_line(size = 0.8),
                axis.line.y  = element_line(size = 0.8),
                legend.position = 'right', legend.text = element_text(size = 15,  face ="plain"),
                legend.title = element_blank()) 

# Horizontal
pdf(file = paste(base.dir, '/Output/Figure/Fig6b_GCI-.pdf', sep = ''),
    width = 10, height = 10)
par(mar = c(2, 2, 2, 2))
ggplot(data = FDR.q$De, aes(x=period, y=fdr, color = type) ) + 
  geom_bar(stat = "identity", position='dodge',alpha=0.8, width = 0.8,
           aes(x=period, color = type, fill = type)) +
  # scale_y_continuous('-log10(q)') + 
  scale_color_manual(values = c( "#0218A8","#475fd0" , "#7699f6", "#b6cefa"), 
                     labels =c("pSI=0.05", "pSI=0.01", "pSI=0.001", "pSI=0.001")) +
  scale_fill_manual(values = c( "#0218A8","#475fd0" , "#7699f6", "#b6cefa"),
                    labels =c("pSI=0.05", "pSI=0.01", "pSI=0.001", "pSI=0.001"))+
  # scale_x_discrete("Developmental period",limits=rev(rownames(Q$All)))+
  # scale_fill_discrete(guide = FALSE) +
  coord_flip() +
  scale_y_continuous("-log(q)",expand = c(0,0)) +
  geom_hline(yintercept = -log10(0.05),lty="dashed", show.legend = TRUE) +
  scale_x_discrete("Developmental period",labels = function(x) str_wrap(x, width = 13), limits=rev(rownames(Q$All) )) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 18,  face ="plain"),
        axis.title.y = element_text(size = 18,  face ="plain"),
        axis.text.x  = element_text(size = 15,  face ="plain",  hjust = 1),
        axis.text.y  = element_text(size = 15,  face ="plain"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
        plot.title = element_text(size = 20,  face ="bold", hjust = 0.5),
        axis.line.x  = element_line(size = 0.8),
        axis.line.y  = element_line(size = 0.8),
        legend.text = element_text(size = 15,  face ="plain"),
        legend.title = element_blank())
dev.off()



pdf(file = paste(base.dir, '/Output/Figure/Fig6a_GCI+.pdf', sep = ''),
    width = 10, height = 10)
par(mar = c(5, 5, 5, 5))
ggplot(data = FDR.q$In, aes(x=period, y=fdr, color = type) ) +
  geom_bar(stat = "identity", position='dodge',alpha=0.8, width = 0.8,
           aes(x=period, color = type, fill = type)) +
  scale_color_manual(values = c( "#b81029","#d95646" , "#ee8568", "#f7b093"), 
                     labels =c("pSI=0.05", "pSI=0.01", "pSI=0.001", "pSI=0.001")) +
  scale_fill_manual(values = c( "#b81029","#d95646" , "#ee8568", "#f7b093"),
                    labels =c("pSI=0.05", "pSI=0.01", "pSI=0.001", "pSI=0.001"))+
  # scale_x_discrete("Developmental period",limits=rev(rownames(Q$All)))+
  # scale_fill_discrete(guide = FALSE) +
  coord_flip() +
  scale_y_continuous("-log(q)",expand = c(0,0)) +
  geom_hline(yintercept = -log10(0.05),lty="dashed", show.legend = TRUE) +
  scale_x_discrete("Developmental period",labels = function(x) str_wrap(x, width = 13), limits=rev(rownames(Q$All) )) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 18,  face ="plain"),
        axis.title.y = element_text(size = 18,  face ="plain"),
        axis.text.x  = element_text(size = 15,  face ="plain",  hjust = 1),
        axis.text.y  = element_text(size = 15,  face ="plain"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
        plot.title = element_text(size = 20,  face ="bold", hjust = 0.5),
        axis.line.x  = element_line(size = 0.8),
        axis.line.y  = element_line(size = 0.8),
        legend.text = element_text(size = 15,  face ="plain"),
        legend.title = element_blank())
dev.off()
write.csv(file = paste(base.dir, '/Output/SupplementaryData/SUPPDATA_sheet28_TemporalSpecificity.csv', sep=''), x = Q$All)
# ------------
