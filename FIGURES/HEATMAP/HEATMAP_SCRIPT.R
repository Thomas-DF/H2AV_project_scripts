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

outfig=paste0(workdir,"FIGURES/HEATMAP/")

#####################################################################################-
#          LOAD DATA  ----
#####################################################################################-

GNref = readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))

LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
Q_H3K36me3_2C4_GB_f=LIST_QUANTIF_K36$Q_H3K36me3_2C4_GB_f

ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f = LIST_QUANTIF_K36$ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f


####### PROFMAT #######

PROFMAT_WT3_NELFKD_R3_RPGC_smth = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_WT3_NELFKD_R3_RPGC_smth.RDS"))
#PROFMAT_WT3_NELFKD_R3_RPGC_smth = readRDS("/work/user/tdefreitas/PROJET_H2AV_2025/DATA/PROFMAT_SMOOTH/PROFMAT_WT3_NELFKD_R3_RPGC_smth.RDS")

PROFMAT_WT3_NELFKD_R2_RPKM_smth = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_WT3_NELFKD_R2_RPKM_smth.RDS"))
PROFMAT_WT3_NELFKD_R3_RPKM_smth = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_WT3_NELFKD_R3_RPKM_smth.RDS"))

PROFMAT_WT3_WT_R3_RPGC_smth = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_WT3_WT_R3_RPGC_smth.RDS"))

PROFMAT_WT3_WT_R3_RPKM_smth = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_WT3_WT_R3_RPKM_smth.RDS"))
PROFMAT_WT3_WT_R2_RPKM_smth = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_WT3_WT_R2_RPKM_smth.RDS"))

PROFMAT_SEA4_NELFKD_R3_RPGC_smth = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_SEA4_NELFKD_R3_RPGC_smth.RDS"))

PROFMAT_SEA4_NELFKD_R3_RPKM_smth = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_SEA4_NELFKD_R3_RPKM_smth.RDS"))
PROFMAT_SEA4_NELFKD_R2_RPKM_smth = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_SEA4_NELFKD_R2_RPKM_smth.RDS"))

PROFMAT_SEA4_WT_R3_RPGC_smth = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_SEA4_WT_R3_RPGC_smth.RDS"))

PROFMAT_SEA4_WT_R3_RPKM_smth = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_SEA4_WT_R3_RPKM_smth.RDS"))
PROFMAT_SEA4_WT_R2_RPKM_smth = readRDS(paste0(workdir,"DATA/PROFMAT_SMOOTH/PROFMAT_SEA4_WT_R2_RPKM_smth.RDS"))



####### ZSCORE #######
 
ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT3_NELFKD_R3 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT3_NELFKD_R3.RDS"))

ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT_R3 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT_R3.RDS"))

ZSCORE_PROFMAT_GB_SEA4_WT_vs_NELFKD_R2 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_SEA4_WT_vs_NELFKD_R2.RDS"))
ZSCORE_PROFMAT_GB_SEA4_WT_vs_NELFKD_R3 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_SEA4_WT_vs_NELFKD_R3.RDS"))

ZSCORE_PROFMAT_GB_SEA4_WT_vs_WT3_WT_R3 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_SEA4_WT_vs_WT3_WT_R3.RDS"))

ZSCORE_PROFMAT_GB_WT3_NELFKD_vs_WT_R3 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_WT3_NELFKD_vs_WT_R3.RDS"))

ZSCORE_PROFMAT_GB_WT3_WT_vs_NELFKD_R2 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_WT3_WT_vs_NELFKD_R2.RDS"))
ZSCORE_PROFMAT_GB_WT3_WT_vs_NELFKD_R3 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_WT3_WT_vs_NELFKD_R3.RDS"))

ZSCORE_PROFMAT_GB_WTH_vs_WT_R3 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_WTH_vs_WT_R3.RDS"))


########### RAD51 #############

ZSCORE_PROFMAT_GB_Rad51_WTH_vs_WT_R2 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_Rad51_WTH_vs_WT_R2.RDS"))
ZSCORE_PROFMAT_GB_Rad51_NH_vs_N_R2 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_Rad51_NH_vs_N_R2.RDS"))
ZSCORE_PROFMAT_GB_Rad51_N_vs_WT_R2 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_Rad51_N_vs_WT_R2.RDS"))
ZSCORE_PROFMAT_GB_Rad51_NH_vs_WTH_R2 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_Rad51_NH_vs_WTH_R2.RDS"))


#####################################################################################-
#         PLOT  ----
#####################################################################################-

####### HEATMAP PROFMAT #######


rangeheatmap = c(1:1000)

# WT3_NELFKD_RPCG

pdf(paste0(outfig,"PROFMAT_WT3_NELFKD_R3_RPGC_smth_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(PROFMAT_WT3_NELFKD_R3_RPGC_smth[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="PROFMAT_WT3_NELFKD_R3_RPGC_smth",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


# WT3_NELFKD_RPKM

pdf(paste0(outfig,"PROFMAT_WT3_NELFKD_R2_R3_RPKM_smth_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(PROFMAT_WT3_NELFKD_R2_RPKM_smth[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="PROFMAT_WT3_NELFKD_R2_RPKM_smth",legend.name="Q_H3K36me3_2C4_GB_f")        
heatMatrixMat(PROFMAT_WT3_NELFKD_R3_RPKM_smth[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="PROFMAT_WT3_NELFKD_R3_RPKM_smth",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


# WT3_WT_RPGC

pdf(paste0(outfig,"PROFMAT_WT3_WT_R3_RPGC_smth_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(PROFMAT_WT3_WT_R3_RPGC_smth[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="PROFMAT_WT3_WT_R3_RPGC_smth",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


# WT3_WT_RPKM

pdf(paste0(outfig,"PROFMAT_WT3_WT_R2_R3_RPKM_smth_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(PROFMAT_WT3_WT_R2_RPKM_smth[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="PROFMAT_WT3_WT_R2_RPKM_smth",legend.name="Q_H3K36me3_2C4_GB_f")        
heatMatrixMat(PROFMAT_WT3_WT_R3_RPKM_smth[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="PROFMAT_WT3_WT_R3_RPKM_smth",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


# SEA4_NELFKD_RPGC

pdf(paste0(outfig,"PROFMAT_SEA4_NELFKD_R3_RPGC_smth_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(PROFMAT_SEA4_NELFKD_R3_RPGC_smth[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="PROFMAT_SEA4_NELFKD_R3_RPGC_smth",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


# SEA4_NELFKD_RPKM

pdf(paste0(outfig,"PROFMAT_SEA4_NELFKD_R2_R3_RPKM_smth_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(PROFMAT_SEA4_NELFKD_R2_RPKM_smth[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="PROFMAT_SEA4_NELFKD_R2_RPKM_smth",legend.name="Q_H3K36me3_2C4_GB_f")        
heatMatrixMat(PROFMAT_SEA4_NELFKD_R3_RPKM_smth[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="PROFMAT_SEA4_NELFKD_R3_RPKM_smth",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


# SEA4_WT_RPGC

pdf(paste0(outfig,"PROFMAT_SEA4_WT_R3_RPGC_smth_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(PROFMAT_SEA4_WT_R3_RPGC_smth[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="PROFMAT_SEA4_WT_R3_RPGC_smth",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


# SEA4_WT_RPKM

pdf(paste0(outfig,"PROFMAT_SEA4_WT_R2_R3_RPKM_smth_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(PROFMAT_SEA4_WT_R2_RPKM_smth[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="PROFMAT_SEA4_WT_R2_RPKM_smth",legend.name="Q_H3K36me3_2C4_GB_f")        
heatMatrixMat(PROFMAT_SEA4_WT_R3_RPKM_smth[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="PROFMAT_SEA4_WT_R3_RPKM_smth",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()



####### HEATMAP ZSCORE #######

rangeheatmap = c(1:1000)


# ZSCORE_SEA4_NELFKD_vs_WT3_NELFKD

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT3_NELFKD_R3_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT3_NELFKD_R3[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],order = TRUE,winsorize=c(5,95),main="ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT3_NELFKD_R3",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


# ZSCORE_SEA4_NELFKD_vs_WT

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT_R3_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT_R3[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT_R3",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


# ZSCORE_SEA4_WT_vs_NELFKD

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_SEA4_WT_vs_NELFKD_R3_BY_Q_H3K36me3_2C4_GB_f.pdf"))
#heatMatrixMat(ZSCORE_PROFMAT_GB_SEA4_WT_vs_NELFKD_R2[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),main="ZSCORE_PROFMAT_GB_SEA4_WT_vs_NELFKD_R2",legend.name="Q_H3K36me3_2C4_GB_f")        
heatMatrixMat(ZSCORE_PROFMAT_GB_SEA4_WT_vs_NELFKD_R3[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="ZSCORE_PROFMAT_GB_SEA4_WT_vs_NELFKD_R3",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


# ZSCORE_SEA4_WT_vs_WT3_WT

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_SEA4_WT_vs_WT3_WT_R3_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_SEA4_WT_vs_WT3_WT_R3[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="ZSCORE_PROFMAT_GB_SEA4_WT_vs_WT3_WT_R3",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


# ZSCORE_WT3_NELFKD_vs_WT

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_WT3_NELFKD_vs_WT_R3_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_WT3_NELFKD_vs_WT_R3[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="ZSCORE_PROFMAT_GB_WT3_NELFKD_vs_WT_R3",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


# ZSCORE_WT3_WT_vs_NELFKD

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_WT3_WT_vs_NELFKD_R3_BY_Q_H3K36me3_2C4_GB_f.pdf"))
#heatMatrixMat(ZSCORE_PROFMAT_GB_WT3_WT_vs_NELFKD_R2[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),main="ZSCORE_PROFMAT_GB_WT3_WT_vs_NELFKD_R2",legend.name="Q_H3K36me3_2C4_GB_f")        
heatMatrixMat(ZSCORE_PROFMAT_GB_WT3_WT_vs_NELFKD_R3[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="ZSCORE_PROFMAT_GB_WT3_WT_vs_NELFKD_R3",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


# ZSCORE_WTH_vs_WT

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_WTH_vs_WT_R3_BY_Q_H3K36me3_2C4_GB_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_WTH_vs_WT_R3[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),order = TRUE,main="ZSCORE_PROFMAT_GB_WTH_vs_WT_R3",legend.name="Q_H3K36me3_2C4_GB_f")        
dev.off()


########  ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f  ###########

rangeheatmap = c(1:1000)

# ZSCORE_SEA4_NELFKD_vs_WT

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT_R3_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))

heatMatrixMat(ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT_R3[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f),rangeheatmap],winsorize=c(5,95), order = TRUE, main="ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT_R3",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
heatMatrixMat(ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT_R3[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f),rangeheatmap],winsorize=c(5,95),RangeValue = c(-0.5,0.7),order = TRUE,main="ZSCORE_PROFMAT_GB_SEA4_NELFKD_vs_WT_R3",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()


# ZSCORE_WT3_NELFKD_vs_WT

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_WT3_NELFKD_vs_WT_R3_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))

heatMatrixMat(ZSCORE_PROFMAT_GB_WT3_NELFKD_vs_WT_R3[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f),rangeheatmap],winsorize=c(5,95), order = TRUE,main="ZSCORE_PROFMAT_GB_WT3_NELFKD_vs_WT_R3",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")   
heatMatrixMat(ZSCORE_PROFMAT_GB_WT3_NELFKD_vs_WT_R3[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f),rangeheatmap],winsorize=c(5,95),RangeValue = c(-0.5,0.7), order = TRUE,main="ZSCORE_PROFMAT_GB_WT3_NELFKD_vs_WT_R3",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")   

dev.off()


##### RAD51 #####


# ZSCORE WTH_vs_WT

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_Rad51_WTH_vs_WT_R2_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_Rad51_WTH_vs_WT_R2[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f),rangeheatmap],winsorize=c(5,95), order = TRUE, main="ZSCORE_PROFMAT_GB_Rad51_WTH_vs_WT_R2",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()


# ZSCORE NH_vs_N

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_Rad51_NH_vs_N_R2_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_Rad51_NH_vs_N_R2[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f),rangeheatmap],winsorize=c(5,95), order = TRUE, main="ZSCORE_PROFMAT_GB_Rad51_NH_vs_N_R2",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()



# ZSCORE N_vs_WT

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_Rad51_N_vs_WT_R2_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_Rad51_N_vs_WT_R2[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f),rangeheatmap],winsorize=c(5,95), order = TRUE, main="ZSCORE_PROFMAT_GB_Rad51_N_vs_WT_R2",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()


# ZSCORE NH_vs_WTH

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_Rad51_NH_vs_WTH_R2_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_Rad51_NH_vs_WTH_R2[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f),rangeheatmap],winsorize=c(5,95), order = TRUE, main="ZSCORE_PROFMAT_GB_Rad51_NH_vs_WTH_R2",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()








