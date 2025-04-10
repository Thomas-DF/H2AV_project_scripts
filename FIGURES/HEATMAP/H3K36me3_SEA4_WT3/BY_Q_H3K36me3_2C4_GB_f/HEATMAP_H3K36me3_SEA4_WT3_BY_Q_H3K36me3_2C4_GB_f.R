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

outfig=paste0(workdir,"FIGURES/HEATMAP/H3K36me3_SEA4_WT3/BY_Q_H3K36me3_2C4_GB_f/")

#####################################################################################-
#          LOAD DATA  ----
#####################################################################################-

GNref = readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))

LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
Q_H3K36me3_2C4_GB_f=LIST_QUANTIF_K36$Q_H3K36me3_2C4_GB_f


########################
######  CHIPSEQ  #######
########################
ChIP_H3K36me3_FSEA4_Luc_profmat = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/ChIP_H3K36me3_FSEA4_Luc_profmat.RDS"))
ChIP_H3K36me3_FSEA4_Nelf_profmat = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/ChIP_H3K36me3_FSEA4_Nelf_profmat.RDS"))

ChIP_H3K36me3_FWT3_Luc_profmat = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/ChIP_H3K36me3_FWT3_Luc_profmat.RDS"))
ChIP_H3K36me3_FWT3_Nelf_profmat = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/ChIP_H3K36me3_FWT3_Nelf_profmat.RDS"))


########################
######  DRIPSEQ  #######
########################

DRIP_HLuc_profmat = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/DRIP_HLuc_profmat.RDS"))
DRIP_HLuc_RNAse_profmat = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/DRIP_HLuc_RNAse_profmat.RDS"))

DRIP_Nelf_profmat = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/DRIP_Nelf_profmat.RDS"))
DRIP_Nelf_RNAse_profmat = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/DRIP_Nelf_RNAse_profmat.RDS"))

DRIP_HypB_profmat = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/DRIP_HypB_profmat.RDS"))
DRIP_HypB_RNAse_profmat = readRDS(paste0(workdir,"DATA/PROFILE_MATRIX/DRIP_HypB_RNAse_profmat.RDS"))


ZSCORE_PROFMAT_DRIP_HLuc_vs_HLuc_RNAse = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_DRIP_HLuc_vs_HLuc_RNAse.RDS"))
ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse =  readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse.RDS"))


#####################################################################################-
#         PLOT  ----
#####################################################################################-

########################
######  CHIPSEQ  #######
########################

rangeheatmap = c(1:1000)

#############  FSEA4  #############

pdf(paste0(outfig, "ChIP_H3K36me3_FSEA4_Luc_profmat","_BY_","Q_H3K36me3_2C4_GB_f",".pdf"))
heatMatrixMat(ChIP_H3K36me3_FSEA4_Luc_profmat[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ChIP_H3K36me3_FSEA4_Luc_profmat",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()

pdf(paste0(outfig, "ChIP_H3K36me3_FSEA4_Nelf_profmat","_BY_","Q_H3K36me3_2C4_GB_f",".pdf"))
heatMatrixMat(ChIP_H3K36me3_FSEA4_Nelf_profmat[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ChIP_H3K36me3_FSEA4_Nelf_profmat",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


#############  FWT3  #############

pdf(paste0(outfig, "ChIP_H3K36me3_FWT3_Luc_profmat","_BY_","Q_H3K36me3_2C4_GB_f",".pdf"))
heatMatrixMat(ChIP_H3K36me3_FWT3_Luc_profmat[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ChIP_H3K36me3_FWT3_Luc_profmat",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()

pdf(paste0(outfig, "ChIP_H3K36me3_FWT3_Nelf_profmat","_BY_","Q_H3K36me3_2C4_GB_f",".pdf"))
heatMatrixMat(ChIP_H3K36me3_FWT3_Nelf_profmat[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ChIP_H3K36me3_FWT3_Nelf_profmat",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()



########################
######  DRIPSEQ  #######
########################

rangeheatmap = c(1:1000)

#############  HLuc  #############

pdf(paste0(outfig, "DRIP_HLuc_profmat","_BY_","Q_H3K36me3_2C4_GB_f",".pdf"))
heatMatrixMat(DRIP_HLuc_profmat[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="DRIP_HLuc_profmat",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()

pdf(paste0(outfig, "DRIP_HLuc_RNAse_profmat","_BY_","Q_H3K36me3_2C4_GB_f",".pdf"))
heatMatrixMat(DRIP_HLuc_RNAse_profmat[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="DRIP_HLuc_RNAse_profmat",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


#############  NELF  #############

pdf(paste0(outfig, "DRIP_Nelf_profmat","_BY_","Q_H3K36me3_2C4_GB_f",".pdf"))
heatMatrixMat(DRIP_Nelf_profmat[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="DRIP_Nelf_profmat",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()

pdf(paste0(outfig, "ChIP_H3K36me3_FWT3_Nelf_profmat","_BY_","Q_H3K36me3_2C4_GB_f",".pdf"))
heatMatrixMat(ChIP_H3K36me3_FWT3_Nelf_profmat[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ChIP_H3K36me3_FWT3_Nelf_profmat",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


#############  HypB  #############

pdf(paste0(outfig, "DRIP_HypB_profmat","_BY_","Q_H3K36me3_2C4_GB_f",".pdf"))
heatMatrixMat(DRIP_HypB_profmat[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="DRIP_HypB_profmat",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()

pdf(paste0(outfig, "DRIP_HypB_RNAse_profmat","_BY_","Q_H3K36me3_2C4_GB_f",".pdf"))
heatMatrixMat(DRIP_HypB_RNAse_profmat[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="DRIP_HypB_RNAse_profmat",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


#############  ZSCORE  #############

pdf(paste0(outfig, "ZSCORE_PROFMAT_DRIP_HLuc_vs_HLuc_RNAse","_BY_","Q_H3K36me3_2C4_GB_f",".pdf"))
heatMatrixMat(ZSCORE_PROFMAT_DRIP_HLuc_vs_HLuc_RNAse[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_DRIP_HLuc_vs_HLuc_RNAse",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()

pdf(paste0(outfig, "ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse","_BY_","Q_H3K36me3_2C4_GB_f",".pdf"))
heatMatrixMat(ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse[names(Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()






