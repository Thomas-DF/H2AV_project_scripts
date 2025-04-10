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

outfig = paste0(workdir,"FIGURES/BOXPLOT/BOXPLOT_POL2/")


#####################################################################################-
#         DATA  ----
#####################################################################################-

pol2_ctrl_N_RPGC = readRDS(paste0(workdir, "DATA/QUANTIF/pol2_ctrl_N_RPGC_readsCounts_GB_SCALED.RDS"))
pol2_nelf_N_RPGC = readRDS(paste0(workdir, "DATA/QUANTIF/pol2_nelf_N_RPGC_readsCounts_GB_SCALED.RDS"))

#####################################################################################-

##### GENES GROUPES 

PAUSE_INDICE_VEC = readRDS(paste0(workdir, "DATA/LIST_FEATURES/PAUSE_INDICE_VEC.RDS"))
PAUSE_IND_pol2_start_TSS_WT_KD = PAUSE_INDICE_VEC$PAUSE_IND_pol2_start_TSS_WT_KD



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


len = length(PAUSE_IND_pol2_start_TSS_WT_KD)
CTRL_60_65_PAUSE_IND_pol2_start_TSS_WT_KD = names(PAUSE_IND_pol2_start_TSS_WT_KD[(len*0.6) : (len*0.65)])
UP_25_PAUSE_IND_pol2_start_TSS_WT_KD = getNameList(PAUSE_IND_pol2_start_TSS_WT_KD, topdown = "top", prct = 25)
UP_5_PAUSE_IND_pol2_start_TSS_WT_KD = getNameList(PAUSE_IND_pol2_start_TSS_WT_KD, topdown = "top", prct = 5)
UP_1_PAUSE_IND_pol2_start_TSS_WT_KD = getNameList(PAUSE_IND_pol2_start_TSS_WT_KD, topdown = "top", prct = 1)
DN_25_PAUSE_IND_pol2_start_TSS_WT_KD = getNameList(PAUSE_IND_pol2_start_TSS_WT_KD, topdown = "down", prct = 25)
DN_5_PAUSE_IND_pol2_start_TSS_WT_KD = getNameList(PAUSE_IND_pol2_start_TSS_WT_KD, topdown = "down", prct = 5)
DN_1_PAUSE_IND_pol2_start_TSS_WT_KD = getNameList(PAUSE_IND_pol2_start_TSS_WT_KD, topdown = "down", prct = 1)

#### LIST GENES GROUPS

LIST_5_PAUSE_IND_pol2_start_TSS_WT_KD = list(
  DN_5_PAUSE_IND_pol2_start_TSS_WT_KD = DN_5_PAUSE_IND_pol2_start_TSS_WT_KD,
  UP_5_PAUSE_IND_pol2_start_TSS_WT_KD = UP_5_PAUSE_IND_pol2_start_TSS_WT_KD
)

LIST_25_UP_DN_PAUSE_IND_pol2_start_TSS_WT_KD = list(
  DN_25_PAUSE_IND_pol2_start_TSS_WT_KD = DN_25_PAUSE_IND_pol2_start_TSS_WT_KD,
  UP_25_PAUSE_IND_pol2_start_TSS_WT_KD = UP_25_PAUSE_IND_pol2_start_TSS_WT_KD
)

LIST_1_UP_DN_PAUSE_IND_pol2_start_TSS_WT_KD = list(
  DN_1_PAUSE_IND_pol2_start_TSS_WT_KD = DN_1_PAUSE_IND_pol2_start_TSS_WT_KD,
  UP_1_PAUSE_IND_pol2_start_TSS_WT_KD = UP_1_PAUSE_IND_pol2_start_TSS_WT_KD
)


LIST_5_UP_CTRL_PAUSE_IND_pol2_start_TSS_WT_KD = list(
  DN_5_PAUSE_IND_pol2_start_TSS_WT_KD = DN_5_PAUSE_IND_pol2_start_TSS_WT_KD,
  CTRL_60_65_PAUSE_IND_pol2_start_TSS_WT_KD = CTRL_60_65_PAUSE_IND_pol2_start_TSS_WT_KD,
  UP_5_PAUSE_IND_pol2_start_TSS_WT_KD = UP_5_PAUSE_IND_pol2_start_TSS_WT_KD
)


#####################################################################################-
#         PLOT  ----
#####################################################################################-


Boxplot_wilcoxListFilter_REF(quantifWT = pol2_ctrl_N_RPGC, quantifKD = pol2_nelf_N_RPGC, cond1 = "pol2_ctrl_N_RPGC", cond2="pol2_nelf_N_RPGC", filterGNList = LIST_5_PAUSE_IND_pol2_start_TSS_WT_KD,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "pol2_ctrl_N_RPGC", Cond = "pol2_nelf_N_RPGC", select = "LIST_5_PAUSE_IND_pol2_start_TSS_WT_KD", info = NULL)

Boxplot_wilcoxListFilter_REF(quantifWT = pol2_ctrl_N_RPGC, quantifKD = pol2_nelf_N_RPGC, cond1 = "pol2_ctrl_N_RPGC", cond2="pol2_nelf_N_RPGC", filterGNList = LIST_25_UP_DN_PAUSE_IND_pol2_start_TSS_WT_KD,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "pol2_ctrl_N_RPGC", Cond = "pol2_nelf_N_RPGC", select = "LIST_25_UP_DN_PAUSE_IND_pol2_start_TSS_WT_KD", info = NULL)

Boxplot_wilcoxListFilter_REF(quantifWT = pol2_ctrl_N_RPGC, quantifKD = pol2_nelf_N_RPGC, cond1 = "pol2_ctrl_N_RPGC", cond2="pol2_nelf_N_RPGC", filterGNList = LIST_1_UP_DN_PAUSE_IND_pol2_start_TSS_WT_KD,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "pol2_ctrl_N_RPGC", Cond = "pol2_nelf_N_RPGC", select = "LIST_1_UP_DN_PAUSE_IND_pol2_start_TSS_WT_KD", info = NULL)

Boxplot_wilcoxListFilter_REF(quantifWT = pol2_ctrl_N_RPGC, quantifKD = pol2_nelf_N_RPGC, cond1 = "pol2_ctrl_N_RPGC", cond2="pol2_nelf_N_RPGC", filterGNList = LIST_5_UP_CTRL_PAUSE_IND_pol2_start_TSS_WT_KD,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "pol2_ctrl_N_RPGC", Cond = "pol2_nelf_N_RPGC", select = "LIST_5_UP_CTRL_PAUSE_IND_pol2_start_TSS_WT_KD", info = NULL)

