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

outfig=paste0(workdir,"FIGURES/HEATMAP/H3K36me3/")

#####################################################################################-
#          LOAD DATA  ----
#####################################################################################-
PROFMAT_GB_H3K36m3_2C4_smth = readRDS("~/Bureau/tdefreitas_genobioinfo/PROJET_H2AV_2025/DATA/PROFMAT_SMOOTH/2C4K36me3_trimmed_filt_sort_RPGC_profmat.RDS")

GNref = readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))
LIST_QUANTIF=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF.RDS"))
ZSCORE_PROFMAT_H2AV_WT_vs_NELF = LIST_QUANTIF$ZSCORE_PROFMAT_H2AV_WT_vs_NELF
rownames(ZSCORE_PROFMAT_H2AV_WT_vs_NELF) <- paste0(rownames(ZSCORE_PROFMAT_H2AV_WT_vs_NELF), ".1")

LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
Q_H3K36me3_2C4_GB_f=LIST_QUANTIF_K36$Q_H3K36me3_2C4_GB_f

ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f = LIST_QUANTIF_K36$ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f
ZSCORE_PROFMAT_GB_H3K36m3_2C4_vs_2N4 = LIST_QUANTIF_K36$ZSCORE_PROFMAT_GB_H3K36m3_2C4_vs_2N4

#####################################################################################-
#         PLOT  ----
#####################################################################################-

rangeheatmap = c(1:1000)

#############  HLuc  #############

pdf(paste0(outfig, "ZSCORE_PROFMAT_GB_H3K36m3_2C4_vs_2N4","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_H3K36m3_2C4_vs_2N4[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_GB_H3K36m3_2C4_vs_2N4",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()


pdf(paste0(outfig, "PROFMAT_GB_H3K36m3_2C4_smth","_BY_","Q_H3K36me3_2C4_GB_f",".pdf"))
heatMatrixMat(PROFMAT_GB_H3K36m3_2C4_smth[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="PROFMAT_GB_H3K36m3_2C4_smth",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()

pdf(paste0(outfig, "ZSCORE_PROFMAT_H2AV_WT_vs_NELF","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(ZSCORE_PROFMAT_H2AV_WT_vs_NELF[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_H2AV_WT_vs_NELF",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()
