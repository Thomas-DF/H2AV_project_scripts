#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

library(ggpubr)
library(rstatix)
library("gplots")
library(textplot)
library("ggplot2")
library(dplyr)
'%ni%' = Negate('%in%')

#####################################################################################-
#         FUNCTIONS  ----
#####################################################################################-

workdir = "~/Bureau/tdefreitas_genobioinfo/PROJET_H2AV_2025/"
source(paste0(workdir,"functionR/Boxplot_wilcoxListFilter_REF.R"))

outfig = paste0(workdir,"FIGURES/BOXPLOT/BOXPLOT_RATIO_SER2P_POL2/")


#####################################################################################-
#         DATA  ----
#####################################################################################-

pol2_nelf_N = readRDS(paste0(workdir,"DATA/QUANTIF/Q_pol2_nelf_N_filt_sort_RPGC_profmat_readsCounts_500_1000.RDS"))
pol2_ctrl_N = readRDS(paste0(workdir,"DATA/QUANTIF/Q_pol2_ctrl_N_filt_sort_RPGC_profmat_readsCounts_500_1000.RDS"))

pol2_ser2P_nelf_N = readRDS(paste0(workdir,"DATA/QUANTIF/Q_pol2ser2P_nelf_N_filt_sort_RPGC_profmat_readsCounts_500_1000.RDS"))
pol2_ser2P_ctrl_N = readRDS(paste0(workdir,"DATA/QUANTIF/Q_pol2ser2P_ctrl_N_filt_sort_RPGC_profmat_readsCounts_500_1000.RDS"))

ratio_ctrl_pol2S2P_pol2 = pol2_ser2P_ctrl_N/(pol2_ctrl_N + pol2_ser2P_ctrl_N)
ratio_nelf_pol2S2P_pol2 = pol2_ser2P_nelf_N/(pol2_nelf_N + pol2_ser2P_nelf_N)


ratio_ctrl_pol2S2P_pol2[is.na(ratio_ctrl_pol2S2P_pol2)] = 0
ratio_nelf_pol2S2P_pol2[is.na(ratio_nelf_pol2S2P_pol2)] = 0

#####################################################################################-

##### GENES GROUPES 

PAUSE_INDICE_VEC = readRDS(paste0(workdir, "DATA/LIST_FEATURES/PAUSE_INDICE_VEC.RDS"))
PAUSE_IND_pol2_start_TSS_WT_KD = PAUSE_INDICE_VEC$PAUSE_IND_pol2_start_TSS_WT_KD

len = length(PAUSE_IND_pol2_start_TSS_WT_KD)
QUARTILE_75_100 = names(PAUSE_IND_pol2_start_TSS_WT_KD[(0) : (len*0.25)])
QUARTILE_50_75 = names(PAUSE_IND_pol2_start_TSS_WT_KD[(len*0.25) : (len*0.5)])
QUARTILE_25_50 = names(PAUSE_IND_pol2_start_TSS_WT_KD[(len*0.5) : (len*0.75)])
QUARTILE_0_25 = names(PAUSE_IND_pol2_start_TSS_WT_KD[(len*0.75) : (len)])

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

UP_5_PAUSE_IND_pol2_start_TSS_WT_KD = getNameList(PAUSE_IND_pol2_start_TSS_WT_KD, topdown = "top", prct = 5)
DN_5_PAUSE_IND_pol2_start_TSS_WT_KD = getNameList(PAUSE_IND_pol2_start_TSS_WT_KD, topdown = "down", prct = 5)
#### LIST GENES GROUPS

LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD =list(
  QUARTILE_75_100 = QUARTILE_75_100,
  QUARTILE_50_75 = QUARTILE_50_75,
  QUARTILE_25_50 = QUARTILE_25_50,
  QUARTILE_0_25 = QUARTILE_0_25)


LIST_PAUSE_IND_pol2_start_TSS_WT_KD =list(
  DN_5_PAUSE_IND_pol2_start_TSS_WT_KD = DN_5_PAUSE_IND_pol2_start_TSS_WT_KD,
  UP_5_PAUSE_IND_pol2_start_TSS_WT_KD = UP_5_PAUSE_IND_pol2_start_TSS_WT_KD)


#####################################################################################-
#         PLOT  ----
#####################################################################################-


Boxplot_wilcoxListFilter_REF(quantifWT = ratio_ctrl_pol2S2P_pol2, quantifKD = ratio_nelf_pol2S2P_pol2, cond1 = "ratio_ctrl_pol2S2P_pol2", cond2="ratio_nelf_pol2S2P_pol2", filterGNList = LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM = c(0,1), bxplt_color = c("#545454", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "ratio_ctrl_pol2S2P_pol2", Cond = "ratio_nelf_pol2S2P_pol2", select = "LIST_QUARTILE_PAUSE_IND_pol2_start_TSS_WT_KD", info = NULL)



Boxplot_wilcoxListFilter_REF(quantifWT = ratio_ctrl_pol2S2P_pol2, quantifKD = ratio_nelf_pol2S2P_pol2, cond1 = "ratio_ctrl_pol2S2P_pol2", cond2="ratio_nelf_pol2S2P_pol2", filterGNList = LIST_PAUSE_IND_pol2_start_TSS_WT_KD,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM = c(0,1), bxplt_color = c("#545454", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "ratio_ctrl_pol2S2P_pol2", Cond = "ratio_nelf_pol2S2P_pol2", select = "LIST_PAUSE_IND_pol2_start_TSS_WT_KD", info = NULL)


