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

outfig = paste0(workdir,"FIGURES/BOXPLOT/BOXPLOT_Q_H2AV/")


#####################################################################################-
#         DATA  ----
#####################################################################################-
Q_H2AV_GB_PN_R1 = readRDS(paste0(workdir,"DATA/QUANTIF/Q_H2AV_GB_PN_R1_f.RDS"))
Q_H2AV_GB_PN_R2 = readRDS(paste0(workdir,"DATA/QUANTIF/Q_H2AV_GB_PN_R2_f.RDS"))

Q_H2AV_GB_PW_R1 = readRDS(paste0(workdir,"DATA/QUANTIF/Q_H2AV_GB_PW_R1_f.RDS"))
Q_H2AV_GB_PW_R2 = readRDS(paste0(workdir,"DATA/QUANTIF/Q_H2AV_GB_PW_R2_f.RDS"))



#####################################################################################-

##### GENES GROUPES 

LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))

Q_H3K36me3_2C4_GB_f=LIST_QUANTIF_K36$Q_H3K36me3_2C4_GB_f

ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f = LIST_QUANTIF_K36$ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f


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


DN10_Q_H3K36me3_2C4=getNameList(Q_H3K36me3_2C4_GB_f, topdown = "down", prct = 10)
UP10_Q_H3K36me3_2C4=getNameList(Q_H3K36me3_2C4_GB_f, topdown = "top", prct = 10)

len = length(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f)

CTRL5_ZSCORE_H3K36me3_2C4_vs_2N4 = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0.6) : (len*0.65)])

DN5_ZSCORE_H3K36me3_2C4_vs_2N4 = getNameList(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, topdown = "down", prct = 5)
UP5_ZSCORE_H3K36me3_2C4_vs_2N4 = getNameList(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, topdown = "top", prct = 5)


DN1_ZSCORE_H3K36me3_2C4_vs_2N4 = getNameList(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, topdown = "down", prct = 1)
UP1_ZSCORE_H3K36me3_2C4_vs_2N4 = getNameList(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, topdown = "top", prct = 1)

#### LIST GENES GROUPS

LIST_5_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4=list(
  DN5_ZSCORE_H3K36me3_2C4_vs_2N4=DN5_ZSCORE_H3K36me3_2C4_vs_2N4,
  delta_UP5_ZSCORE_H3K36me3_2C4_vs_2N4=UP5_ZSCORE_H3K36me3_2C4_vs_2N4
)

LIST_1_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4=list(
  DN1_ZSCORE_H3K36me3_2C4_vs_2N4=DN1_ZSCORE_H3K36me3_2C4_vs_2N4,
  delta_UP1_ZSCORE_H3K36me3_2C4_vs_2N4=UP1_ZSCORE_H3K36me3_2C4_vs_2N4
)


LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4=list(
  DN5_H3K36me3_2C4_vs_2N4=DN5_ZSCORE_H3K36me3_2C4_vs_2N4,
  CTRL_60_65_H3K36me3_2C4_vs_2N4 = CTRL5_ZSCORE_H3K36me3_2C4_vs_2N4,
  delta_UP5_H3K36me3_2C4_vs_2N4 = UP5_ZSCORE_H3K36me3_2C4_vs_2N4
)


#####################################################################################-
#         PLOT  ----
#####################################################################################-


Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_GB_PW_R1, quantifKD = Q_H2AV_GB_PN_R1, cond1 = "Q_H2AV_GB_PW_R1", cond2="Q_H2AV_GB_PN_R1", filterGNList = LIST_5_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_H2AV_GB_PW_R1", Cond = "Q_H2AV_GB_PN_R1", select = "LIST_5_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_GB_PW_R1, quantifKD = Q_H2AV_GB_PN_R1, cond1 = "Q_H2AV_GB_PW_R1", cond2="Q_H2AV_GB_PN_R1", filterGNList = LIST_1_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_H2AV_GB_PW_R1", Cond = "Q_H2AV_GB_PN_R1", select = "LIST_1_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_GB_PW_R1, quantifKD = Q_H2AV_GB_PN_R1, cond1 = "Q_H2AV_GB_PW_R1", cond2="Q_H2AV_GB_PN_R1", filterGNList = LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_H2AV_GB_PW_R1", Cond = "Q_H2AV_GB_PN_R1", select = "LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)



Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_GB_PW_R2, quantifKD = Q_H2AV_GB_PN_R2, cond1 = "Q_H2AV_GB_PW_R2", cond2="Q_H2AV_GB_PN_R2", filterGNList = LIST_5_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_H2AV_GB_PW_R2", Cond = "Q_H2AV_GB_PN_R2", select = "LIST_5_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_GB_PW_R2, quantifKD = Q_H2AV_GB_PN_R2, cond1 = "Q_H2AV_GB_PW_R2", cond2="Q_H2AV_GB_PN_R2", filterGNList = LIST_1_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_H2AV_GB_PW_R2", Cond = "Q_H2AV_GB_PN_R2", select = "LIST_1_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_GB_PW_R2, quantifKD = Q_H2AV_GB_PN_R2, cond1 = "Q_H2AV_GB_PW_R2", cond2="Q_H2AV_GB_PN_R2", filterGNList = LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_H2AV_GB_PW_R2", Cond = "Q_H2AV_GB_PN_R2", select = "LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)
