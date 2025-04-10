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

outfig=paste0(workdir,"FIGURES/HEATMAP/H2AV_SEA4_WT3/BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f/")

#####################################################################################-
#          LOAD DATA  ----
#####################################################################################-

GNref = readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))

LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f = LIST_QUANTIF_K36$ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f


####### ZSCORE #######
 
ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT_R3 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT_R3.RDS"))
ZSCORE_PROFMAT_GB_WT3_NELFKD_vs_WT_R3 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_WT3_NELFKD_vs_WT_R3.RDS"))


#####################################################################################-
#         PLOT  ----
#####################################################################################-

####### HEATMAP ZSCORE #######

rangeheatmap = c(1:1000)


# ZSCORE_SEA4_NELFKD_vs_WT_R3

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT_R3_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT_R3[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],order = FALSE,winsorize=c(5,95),main="ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT_R3",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()


# ZSCORE_WT3_NELFKD_vs_WT

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_WT3_NELFKD_vs_WT_R3_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_WT3_NELFKD_vs_WT_R3[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95),order = FALSE,main="ZSCORE_PROFMAT_GB_WT3_NELFKD_vs_WT_R3",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()






