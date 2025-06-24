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

outfig=paste0(workdir,"FIGURES/HEATMAP/H3K36me3_SEA4_WT3/BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f/")

#####################################################################################-
#          LOAD DATA  ----
#####################################################################################-

GNref = readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))

LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
LIST_QUANTIF=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF.RDS"))

ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f = LIST_QUANTIF_K36$ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f
ZSCORE_H2AV_GB_WT_N = LIST_QUANTIF$ZSCORE_H2AV_GB_WT_N



getNameList = function(Vec, topdown = "top", prct = 10){
  Vec = Vec[order(Vec, decreasing=T)]
  if(topdown %in% "top"){
    GN = names(Vec[Vec > quantile(Vec, (100-prct)/100)])
  }
  if(topdown %in% "down"){
    GN = names(Vec[Vec < quantile(Vec, (prct)/100)])
  }
  if(topdown %in% "mid"){
    tmp1 = names(Vec[Vec < quantile(Vec, (100/2-prct/2)/100)])
    tmp2 = names(Vec[Vec < quantile(Vec, (100/2-prct/2+prct)/100)])
    GN = tmp2[tmp2 %ni% tmp1]
  }
  return(GN)
}


TOP5_ZSCORE_H2AV_GB_WT_N = getNameList(ZSCORE_H2AV_GB_WT_N, topdown = "down", prct = 5)

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

pdf(paste0(outfig, "ChIP_H3K36me3_FSEA4_Luc_profmat","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(ChIP_H3K36me3_FSEA4_Luc_profmat[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ChIP_H3K36me3_FSEA4_Luc_profmat",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()

pdf(paste0(outfig, "ChIP_H3K36me3_FSEA4_Nelf_profmat","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(ChIP_H3K36me3_FSEA4_Nelf_profmat[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ChIP_H3K36me3_FSEA4_Nelf_profmat",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()


#############  FWT3  #############

pdf(paste0(outfig, "ChIP_H3K36me3_FWT3_Luc_profmat","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(ChIP_H3K36me3_FWT3_Luc_profmat[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ChIP_H3K36me3_FWT3_Luc_profmat",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()

pdf(paste0(outfig, "ChIP_H3K36me3_FWT3_Nelf_profmat","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(ChIP_H3K36me3_FWT3_Nelf_profmat[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ChIP_H3K36me3_FWT3_Nelf_profmat",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()



########################
######  DRIPSEQ  #######
########################

rangeheatmap = c(1:1000)

#############  HLuc  #############

pdf(paste0(outfig, "DRIP_HLuc_profmat","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(DRIP_HLuc_profmat[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="DRIP_HLuc_profmat",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()

pdf(paste0(outfig, "DRIP_HLuc_RNAse_profmat","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(DRIP_HLuc_RNAse_profmat[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="DRIP_HLuc_RNAse_profmat",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()


#############  NELF  #############

pdf(paste0(outfig, "DRIP_Nelf_profmat","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(DRIP_Nelf_profmat[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="DRIP_Nelf_profmat",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()

pdf(paste0(outfig, "ChIP_H3K36me3_FWT3_Nelf_profmat","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(ChIP_H3K36me3_FWT3_Nelf_profmat[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ChIP_H3K36me3_FWT3_Nelf_profmat",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()


#############  HypB  #############

pdf(paste0(outfig, "DRIP_HypB_profmat","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(DRIP_HypB_profmat[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="DRIP_HypB_profmat",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()

pdf(paste0(outfig, "DRIP_HypB_RNAse_profmat","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(DRIP_HypB_RNAse_profmat[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="DRIP_HypB_RNAse_profmat",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()


#############  ZSCORE  #############

pdf(paste0(outfig, "ZSCORE_PROFMAT_DRIP_HLuc_vs_HLuc_RNAse","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(ZSCORE_PROFMAT_DRIP_HLuc_vs_HLuc_RNAse[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_DRIP_HLuc_vs_HLuc_RNAse",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()

pdf(paste0(outfig, "ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()



pdf(paste0(outfig, "ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap], RangeValue = c(-7,5) ,winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()

pdf(paste0(outfig, "TEST_ZSCORE_PROFMAT_DRIP_HLuc_vs_HLuc_RNAse","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(ZSCORE_PROFMAT_DRIP_HLuc_vs_HLuc_RNAse[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap], RangeValue = c(-8,5) ,winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_DRIP_HLuc_vs_HLuc_RNAse",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()







pdf(paste0(outfig, "ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse","_BY_","ZSCORE_H2AV_GB_WT_N",".pdf"))
heatMatrixMat(ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse[names(ZSCORE_H2AV_GB_WT_N[order(ZSCORE_H2AV_GB_WT_N, decreasing=F)]),rangeheatmap], RangeValue = c(-7,5) ,winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse",legend.name="ZSCORE_H2AV_GB_WT_N")        
dev.off()





ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse = ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),]
l=length(rownames(ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse))
Q1 = ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse[1:(l*0.25),]
Q2 = ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse[(l*0.25):(l*0.5),]
Q3 = ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse[(l*0.5):(l*0.75),]
Q4 = ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse[(l*0.75):l,]


# Comptage du nombre de séquences TOP5 dans chaque quartile
sum(rownames(Q1) %in% TOP5_ZSCORE_H2AV_GB_WT_N)
sum(rownames(Q2) %in% TOP5_ZSCORE_H2AV_GB_WT_N)
sum(rownames(Q3) %in% TOP5_ZSCORE_H2AV_GB_WT_N)
sum(rownames(Q4) %in% TOP5_ZSCORE_H2AV_GB_WT_N)





ZSCORE_H2AV_GB_WT_N
pdf(paste0(outfig, "ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse_TOP5_ZSCORE_H2AV_GB_WT_N","_BY_","ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f",".pdf"))
heatMatrixMat(filtered_mat, winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse_TOP5_ZSCORE_H2AV_GB_WT_N",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()

#########
# Tri du ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse pour le TOP5_ZSCORE_H2AV_GB_WT_N selon ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f 
common_names = intersect(names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f),rownames(ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse_TOP5_ZSCORE_H2AV_GB_WT_N))
ordered_names = common_names[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[common_names], decreasing = TRUE)]
filtered_mat = ZSCORE_PROFMAT_DRIP_Nelf_vs_Nelf_RNAse_TOP5_ZSCORE_H2AV_GB_WT_N[ordered_names, rangeheatmap]
##########







len = length(ZSCORE_H2AV_GB_WT_N)
TOP5 = ZSCORE_H2AV_GB_WT_N[(len*0.95):len]
TOP5_H3 = ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[1:(len*0.05)]


genes_all = rownames(ZSCORE_H2AV_GB_WT_N)  # ou un vecteur de tous les gènes considérés

top5_H2AV = names(TOP5)
top5_H3K36me3 = names(TOP5_H3)

# Calcul des catégories
a = length(intersect(top5_H2AV, top5_H3K36me3))  # dans les deux
b = length(setdiff(top5_H2AV, top5_H3K36me3))    # dans H2AV seulement
c = length(setdiff(top5_H3K36me3, top5_H2AV))    # dans H3K36me3 seulement
d = length(setdiff(genes_all, union(top5_H2AV, top5_H3K36me3)))  # dans aucun

# Tableau de contingence
contingency_table = matrix(c(a, b, c, d), nrow = 2)
fisher.test(contingency_table)
contingency_table

# pval = 2.2e-16
# odds ratio = 0