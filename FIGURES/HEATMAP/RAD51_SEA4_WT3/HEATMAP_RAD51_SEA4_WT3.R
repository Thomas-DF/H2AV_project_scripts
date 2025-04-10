#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

library(GenomicRanges)
library(ggplot2)
library(ggpubr) 
require(BiocGenerics)
library(rtracklayer)
library(gridBase)

#####################################################################################-
#         FUNCTIONS  ----
#####################################################################################-

workdir = "~/Bureau/tdefreitas_genobioinfo/PROJET_H2AV_2025/"
source(paste0(workdir,"functionR/Script_HEATMAP_profile.R"))

outfig=paste0(workdir,"FIGURES/HEATMAP/RAD51_SEA4_WT3/")

#####################################################################################-
#          LOAD DATA  ----
#####################################################################################-

GNref = readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))

LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
Q_H3K36me3_2C4_GB_f=LIST_QUANTIF_K36$Q_H3K36me3_2C4_GB_f


####### PROFMAT #######

SEA4_Rad51_L1 = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/SEA4_Rad51_L1_RPGC_profmat.RDS"))
SEA4_Rad51_L2 = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/SEA4_Rad51_L2_RPGC_profmat.RDS"))
SEA4_Rad51_N1 = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/SEA4_Rad51_N1_RPGC_profmat.RDS"))
SEA4_Rad51_N2 = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/SEA4_Rad51_N2_RPGC_profmat.RDS"))


WT3_Rad51_L1 = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/WT3_Rad51_L1_RPGC_profmat.RDS"))
WT3_Rad51_L2 = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/WT3_Rad51_L2_RPGC_profmat.RDS"))
WT3_Rad51_N1 = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/WT3_Rad51_N1_RPGC_profmat.RDS"))
WT3_Rad51_N2 = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/WT3_Rad51_N2_RPGC_profmat.RDS"))



#####################################################################################-
#         PLOT  ----
#####################################################################################-

rangeheatmap = c(1:1000)

##### RAD51 #####

#### PROFILE MATRIX SEA4

pdf(paste0(outfig,"SEA4_Rad51_L1_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(SEA4_Rad51_L1[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="SEA4_Rad51_L1",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()

pdf(paste0(outfig,"SEA4_Rad51_L2_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(SEA4_Rad51_L2[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="SEA4_Rad51_L2",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()

pdf(paste0(outfig,"SEA4_Rad51_N1_bis_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(SEA4_Rad51_N1[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="SEA4_Rad51_N1",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()

pdf(paste0(outfig,"SEA4_Rad51_N2_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(SEA4_Rad51_N2[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="SEA4_Rad51_N2",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()




#### PROFILE MATRIX WT3

pdf(paste0(outfig,"WT3_Rad51_L1_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(WT3_Rad51_L1[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="WT3_Rad51_L1",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()

pdf(paste0(outfig,"WT3_Rad51_L2_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(WT3_Rad51_L2[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="WT3_Rad51_L2",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()

pdf(paste0(outfig,"WT3_Rad51_N1_bis_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(WT3_Rad51_N1[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="WT3_Rad51_N1",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()

pdf(paste0(outfig,"WT3_Rad51_N2_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(WT3_Rad51_N2[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="WT3_Rad51_N2",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()




