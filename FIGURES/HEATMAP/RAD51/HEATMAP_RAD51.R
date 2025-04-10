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

outfig=paste0(workdir,"FIGURES/HEATMAP/RAD51/")

#####################################################################################-
#          LOAD DATA  ----
#####################################################################################-

GNref = readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))

LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f = LIST_QUANTIF_K36$ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f


####### PROFMAT #######
PROFMAT_GB_Rad51_WT_bis = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_GB_Rad51_WT_bis_RPGC_smth.RDS"))
PROFMAT_GB_Rad51_N_bis = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_GB_Rad51_N_bis_RPGC_smth.RDS"))
PROFMAT_GB_Rad51_WT = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_GB_Rad51_WT_RPGC_smth.RDS"))
PROFMAT_GB_Rad51_N = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_GB_Rad51_N_RPGC_smth.RDS"))

########### RAD51 #############

ZSCORE_PROFMAT_GB_Rad51_WTH_vs_WT_R2 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_Rad51_WTH_vs_WT_R2.RDS"))
ZSCORE_PROFMAT_GB_Rad51_NH_vs_N_R2 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_Rad51_NH_vs_N_R2.RDS"))
ZSCORE_PROFMAT_GB_Rad51_N_vs_WT_R2 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_Rad51_N_vs_WT_R2.RDS"))
ZSCORE_PROFMAT_GB_Rad51_NH_vs_WTH_R2 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_Rad51_NH_vs_WTH_R2.RDS"))


#####################################################################################-
#         PLOT  ----
#####################################################################################-

rangeheatmap = c(1:1000)

##### RAD51 #####

#### PROFMAT_SMOOTH

pdf(paste0(outfig,"PROFMAT_GB_Rad51_WT_bis_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
heatMatrixMat(PROFMAT_GB_Rad51_WT_bis[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="PROFMAT_GB_Rad51_WT_bis",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()

pdf(paste0(outfig,"PROFMAT_GB_Rad51_WT_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
heatMatrixMat(PROFMAT_GB_Rad51_WT[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="PROFMAT_GB_Rad51_WT",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()

pdf(paste0(outfig,"PROFMAT_GB_Rad51_N_bis_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
heatMatrixMat(PROFMAT_GB_Rad51_N_bis[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="PROFMAT_GB_Rad51_N_bis",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()

pdf(paste0(outfig,"PROFMAT_GB_Rad51_N_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
heatMatrixMat(PROFMAT_GB_Rad51_N[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="PROFMAT_GB_Rad51_N",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()

#### ZSCORE

# ZSCORE WTH_vs_WT

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_Rad51_WTH_vs_WT_R2_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_Rad51_WTH_vs_WT_R2[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_GB_Rad51_WTH_vs_WT_R2",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()


# ZSCORE NH_vs_N

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_Rad51_NH_vs_N_R2_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_Rad51_NH_vs_N_R2[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_GB_Rad51_NH_vs_N_R2",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()


# ZSCORE N_vs_WT

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_Rad51_N_vs_WT_R2_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_Rad51_N_vs_WT_R2[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_GB_Rad51_N_vs_WT_R2",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()


# ZSCORE NH_vs_WTH

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_Rad51_NH_vs_WTH_R2_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_Rad51_NH_vs_WTH_R2[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_GB_Rad51_NH_vs_WTH_R2",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()








