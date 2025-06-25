
#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

library(GenomicRanges)
library(ggplot2)
library(ggpubr) # magrittr
library(gplots)
'%ni%' = Negate('%in%')
require(BiocGenerics)
require(parallel)
library(gsubfn) # proto
library(Rsamtools)
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library(rtracklayer)
library(nucleR)
library(gridBase)

######## PATH 
workdir="/home/cperrois/work/"
outfig="/home/cperrois/work/PROJET_H2AV_2021/FIGURE/HEATMAPS/HEATMAP_H3K36me3/BY_H3K36/By_Q_H3K36m3/"

source(paste0(workdir, "functionR/Script_HEATMAP_profile.R"))

######## DATA 
GNref = readRDS(paste0(workdir, "PROJET_H2AV_2021/DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))


ZSCORE_PROFMAT_GB_H3K36m3_N_vs_L_bis=readRDS("/home/cperrois/work/PROJET_H2AV_2021/DATA/LIST_FEATURES/ZSCORE/ZSCORE_PROFMAT_GB_me3_N_vs_me3_L_bis.RDS")
ZSCORE_PROFMAT_GB_H3K36m3_2N4_vs_2C4=readRDS("/home/cperrois/work/PROJET_H2AV_2021/DATA/LIST_FEATURES/ZSCORE/ZSCORE_PROFMAT_GB_H3K36m3_2N4_vs_2C4.RDS")

######## HEATMAP 

rangeheatmap = c(1:1000)
pdf(paste0(outfig,"ZSCORE_PROMAT_GB_H3K36m3_WT_N_BY_Q_H3K36me3_2C4_GB_f_upd5K.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_H3K36m3_2N4_vs_2C4[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),main="ZSCORE_PROFMAT_GB_H3K36m3_2N4_vs_2C4",legend.name="Q_H3K36me3_2C4_GB_f")        
heatMatrixMat(ZSCORE_PROFMAT_GB_H3K36m3_N_vs_L_bis[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),main="ZSCORE_PROFMAT_GB_me3_N_vs_me3_L_bis",legend.name="Q_H3K36me3_2C4_GB_f")        

dev.off()





