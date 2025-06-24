#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(gplots)
'%ni%' = Negate('%in%')
require(BiocGenerics)
require(parallel)
library(gsubfn)
library(Rsamtools)
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library(rtracklayer)
library(nucleR)
#install.packages("ggforce")
library("ggforce")

#####################################################################################-
#         FUNCTIONS  ----
#####################################################################################-

workdir = "~/Bureau/tdefreitas_genobioinfo/PROJET_H2AV_2025/"
source(paste0(workdir,"functionR/SCATTER_2features.R"))
outfig = paste0(workdir,"FIGURES/SCATTER_PLOT/")

#####################################################################################-
#         DATA  ----
#####################################################################################-
Q_H2AV_R3_HypB_KD = readRDS(paste0(workdir,"DATA/QUANTIF/Q_H2AV_R3_HypB_KD_CPM_profmat_readsCounts_GB_SCALED.RDS"))
Q_H2AV_HypB_KD = readRDS(paste0(workdir,"DATA/QUANTIF/Q_H2AV_PHYPB_A_L1_RPGC_upstr500_dnstr500.RDS"))
names(Q_H2AV_HypB_KD) = paste0(names(Q_H2AV_HypB_KD),".1")
Q_H2AV_GB_PN_R1 = readRDS(paste0(workdir,"DATA/QUANTIF/Q_H2AV_GB_PN_R1_f.RDS"))
Q_H2AV_GB_PN_R2 = readRDS(paste0(workdir,"DATA/QUANTIF/Q_H2AV_GB_PN_R2_f.RDS"))

Q_H2AV_R3_HypB_KD = Q_H2AV_R3_HypB_KD[names(Q_H2AV_R3_HypB_KD) %in% GNref]
Q_H2AV_HypB_KD = Q_H2AV_HypB_KD[names(Q_H2AV_HypB_KD) %in% GNref]

pol2_nelf_N = readRDS(paste0(workdir,"DATA/QUANTIF/Q_pol2_nelf_N_filt_sort_RPGC_profmat_readsCounts_1000_1500.RDS"))
pol2_ctrl_N = readRDS(paste0(workdir,"DATA/QUANTIF/Q_pol2_ctrl_N_filt_sort_RPGC_profmat_readsCounts_1000_1500.RDS"))

pol2_ser2P_nelf_N = readRDS(paste0(workdir,"DATA/QUANTIF/Q_pol2ser2P_nelf_N_filt_sort_RPGC_profmat_readsCounts_1000_1500.RDS"))
pol2_ser2P_ctrl_N = readRDS(paste0(workdir,"DATA/QUANTIF/Q_pol2ser2P_ctrl_N_filt_sort_RPGC_profmat_readsCounts_1000_1500.RDS"))

ratio_ctrl_pol2S2P_pol2 = pol2_ser2P_ctrl_N/(pol2_ctrl_N + pol2_ser2P_ctrl_N)
ratio_nelf_pol2S2P_pol2 = pol2_ser2P_nelf_N/(pol2_nelf_N + pol2_ser2P_nelf_N)

ratio_ctrl_pol2S2P_pol2[is.na(ratio_ctrl_pol2S2P_pol2)] = 0
ratio_nelf_pol2S2P_pol2[is.na(ratio_nelf_pol2S2P_pol2)] = 0

diff_ratio_ctrl_nelf_pol2S2P_pol2 = ratio_ctrl_pol2S2P_pol2 - ratio_nelf_pol2S2P_pol2

#####################################################################################-

##### GENES GROUPES 


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

LIST_QUANTIF=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF.RDS"))
ZSCORE_H2AV_NELF_WT = LIST_QUANTIF$ZSCORE_H2AV_NELF_WT 
ZSCORE_H2AV_HYPB_WT = LIST_QUANTIF$ZSCORE_H2AV_HYPB_WT
ZSCORE_H2AV_R3_HYPB_WT = LIST_QUANTIF$ZSCORE_H2AV_R3_HYPB_WT


LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
Q_H3K36me3_2C4_GB_f=LIST_QUANTIF_K36$Q_H3K36me3_2C4_GB_f

ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f=LIST_QUANTIF_K36$ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f


PAUSE_INDICE_VEC = readRDS(paste0(workdir, "DATA/LIST_FEATURES/PAUSE_INDICE_VEC.RDS"))
PAUSE_IND_pol2_start_TSS_WT_KD = PAUSE_INDICE_VEC$PAUSE_IND_pol2_start_TSS_WT_KD

ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD = PAUSE_INDICE_VEC$ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD
ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD_BY_WT = PAUSE_INDICE_VEC$ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD_BY_WT

len = length(PAUSE_IND_pol2_start_TSS_WT_KD)
QUARTILE_75_100 = names(PAUSE_IND_pol2_start_TSS_WT_KD[(0) : (len*0.25)])
QUARTILE_50_75 = names(PAUSE_IND_pol2_start_TSS_WT_KD[(len*0.25) : (len*0.5)])
QUARTILE_25_50 = names(PAUSE_IND_pol2_start_TSS_WT_KD[(len*0.5) : (len*0.75)])
QUARTILE_0_25 = names(PAUSE_IND_pol2_start_TSS_WT_KD[(len*0.75) : (len)])


LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD =list(
  QUARTILE_75_100 = QUARTILE_75_100,
  QUARTILE_50_75 = QUARTILE_50_75,
  QUARTILE_25_50 = QUARTILE_25_50,
  QUARTILE_0_25 = QUARTILE_0_25)





compare_cluster_Function(feat1 = ZSCORE_H2AV_NELF_WT, feat2 = ZSCORE_H2AV_R3_HYPB_WT, LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                         namefeat1 = "ZSCORE_H2AV_NELF_WT", namefeat2 = "ZSCORE_H2AV_R3_HYPB_WT",
                         logT=F, outpath = outfig, yl=c(-10,10), xl=c(-20,40))



compare_cluster_Function(feat1 = ZSCORE_H2AV_NELF_WT, feat2 = ZSCORE_H2AV_HYPB_WT, LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                         namefeat1 = "ZSCORE_H2AV_NELF_WT", namefeat2 = "ZSCORE_H2AV_R3_HYPB_WT",
                         logT=F, outpath = outfig, yl=c(-26,26), xl=c(-20,30))


compare_cluster_Function(feat1 = ZSCORE_H2AV_R3_HYPB_WT, feat2 = ZSCORE_H2AV_HYPB_WT, LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                         namefeat1 = "ZSCORE_H2AV_R3_HYPB_WT", namefeat2 = "ZSCORE_H2AV_HYPB_WT",
                         logT=F, outpath = outfig, yl=c(-20,30), xl=c(-20,30))



########################

len = length(Q_H2AV_R3_HypB_KD)
Q_H2AV_R3_HypB_KD <- sort(Q_H2AV_R3_HypB_KD, decreasing = TRUE)
Q_H2AV_R3_HypB_KD_rmoutiler = Q_H2AV_R3_HypB_KD[names(Q_H2AV_R3_HypB_KD[(len*0.005) : (len)])]
Q_H2AV_HypB_KD <- sort(Q_H2AV_HypB_KD, decreasing = TRUE)
Q_H2AV_HypB_KD_rmoutiler = Q_H2AV_HypB_KD[names(Q_H2AV_HypB_KD[(len*0.005) : (len)])]
Q_H2AV_GB_PN_R2 <- sort(Q_H2AV_GB_PN_R2, decreasing = TRUE)
Q_H2AV_GB_PN_R2_rmoutiler = Q_H2AV_GB_PN_R2[names(Q_H2AV_GB_PN_R2[(len*0.005) : (len)])]

compare_cluster_Function(feat1 = Q_H2AV_GB_PN_R2, feat2 = Q_H2AV_GB_PN_R1, LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                         namefeat1 = "Q_H2AV_GB_PN_R2", namefeat2 = "Q_H2AV_GB_PN_R1",
                         logT=F, outpath = outfig, yl=c(0,10000), xl=c(0,10000))



compare_cluster_Function(feat1 = Q_H2AV_GB_PN_R2_rmoutiler, feat2 = Q_H2AV_R3_HypB_KD_rmoutiler, LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                         namefeat1 = "Q_H2AV_GB_PN_R2_rmoutiler", namefeat2 = "Q_H2AV_R3_HypB_KD_rmoutlier",
                         logT=F, outpath = outfig, yl=c(0,4000), xl=c(0,4000))


compare_cluster_Function(feat1 = Q_H2AV_GB_PN_R2, feat2 = Q_H2AV_R3_HypB_KD, LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                         namefeat1 = "Q_H2AV_GB_PN_R2", namefeat2 = "Q_H2AV_R3_HypB_KD",
                         logT=F, outpath = outfig, yl=c(0,10000), xl=c(0,10000))



compare_cluster_Function(feat1 = Q_H2AV_GB_PN_R2_rmoutiler, feat2 = Q_H2AV_HypB_KD_rmoutiler, LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                         namefeat1 = "Q_H2AV_GB_PN_R2_rmoutiler", namefeat2 = "Q_H2AV_HypB_KD_rmoutlier",
                         logT=F, outpath = outfig, yl=c(0,3000), xl=c(0,3000))

compare_cluster_Function(feat1 = Q_H2AV_GB_PN_R2, feat2 = Q_H2AV_HypB_KD, LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                         namefeat1 = "Q_H2AV_GB_PN_R2", namefeat2 = "Q_H2AV_HypB_KD",
                         logT=F, outpath = outfig, yl=c(0,10000), xl=c(0,10000))

#############################################################################################"""

compare_cluster_Function(feat1 = ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, feat2 = ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD_BY_WT,LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                         namefeat1 = "ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f", namefeat2 = "ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD_BY_WT",
                         logT=F, outpath = outfig, yl=c(-25,15), xl=c(-40,75))


compare_cluster_Function(feat1 = ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, feat2 = ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD,LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                         namefeat1 = "ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f", namefeat2 = "ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD",
                         logT=F, outpath = outfig, yl=c(-22,22), xl=c(-40,75))



compare_cluster_Function(feat1 = ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, feat2 = diff_ratio_ctrl_nelf_pol2S2P_pol2,LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                         namefeat1 = "ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f", namefeat2 = "diff_ratio_ctrl_nelf_pol2S2P_pol2",
                         logT=F, outpath = outfig, yl=c(-1,1), xl=c(-40,75))










get_means_by_chunk <- function(list1, list2, chunk_size = 50) {
  mean_list1 <- c()
  mean_list2 <- c()
  all_names <- names(list1)
  total <- length(all_names)
  chunk_number <- 1
  
  for (i in seq(1, total, by = chunk_size)) {
    chunk_indices <- i:min(i + chunk_size - 1, total)
    chunk_names <- all_names[chunk_indices]
    
    mean1 <- mean(list1[chunk_names], na.rm = TRUE)
    
    valid_names <- chunk_names[chunk_names %in% names(list2)]

    mean2 <- mean(list2[valid_names], na.rm = TRUE)

    chunk_name <- paste0("chunk_", chunk_number)
    mean_list1[chunk_name] <- mean1
    mean_list2[chunk_name] <- mean2
    chunk_number <- chunk_number + 1
  }
  
  return(list(mean_list1 = mean_list1,mean_list2 = mean_list2))
}


m_list= get_means_by_chunk(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f,ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD_BY_WT)

compare_cluster_Function(feat1 = m_list$mean_list1, feat2 = m_list$mean_list2,LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                         namefeat1 = "ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f", namefeat2 = "ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD_BY_WT",
                         logT=F,info = "mean_by_50" , outpath = outfig, yl=c(-0.5,1.6), xl=c(-30,40))


m_list_10= get_means_by_chunk(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f,ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD_BY_WT,chunk_size = 10)

compare_cluster_Function(feat1 = m_list_10$mean_list1, feat2 = m_list_10$mean_list2,LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                         namefeat1 = "ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f", namefeat2 = "ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD_BY_WT",
                         logT=F,info = "mean_by_10" ,outpath = outfig, yl=c(-2,3), xl=c(-30,55))




m_list2= get_means_by_chunk(ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD_BY_WT,ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f)

compare_cluster_Function(feat1 = m_list2$mean_list2, feat2 = m_list2$mean_list1,LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                         namefeat1 = "ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f", namefeat2 = "ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD_BY_WT",
                         logT=F,info = "mean_by_50_S2P" , outpath = outfig, yl=c(-15,10), xl=c(-3,5))


m_list2_10= get_means_by_chunk(ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD_BY_WT,ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f,chunk_size = 10)

compare_cluster_Function(feat1 = m_list2_10$mean_list2, feat2 = m_list2_10$mean_list1,LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                         namefeat1 = "ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f", namefeat2 = "ZSCORE_pol2ser2P_2C2_2N2_1000_1500_WT_KD_BY_WT",
                         logT=F,info = "mean_by_10_S2P" , outpath = outfig, yl=c(-15,15), xl=c(-7,15))






