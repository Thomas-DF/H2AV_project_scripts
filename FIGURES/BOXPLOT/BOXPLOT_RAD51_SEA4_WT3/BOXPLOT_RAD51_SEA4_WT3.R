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

outfig = paste0(workdir,"FIGURES/BOXPLOT/BOXPLOT_RAD51_SEA4_WT3/")


#####################################################################################-
#         DATA  ----
#####################################################################################-

Q_WT3_Rad51_N1 = readRDS(paste0(workdir,"DATA/QUANTIF/Q_WT3_Rad51_N1_RPGC_profmat_readsCounts_GB_SCALED.RDS"))
Q_WT3_Rad51_N2 = readRDS(paste0(workdir,"DATA/QUANTIF/Q_WT3_Rad51_N2_RPGC_profmat_readsCounts_GB_SCALED.RDS"))

Q_WT3_Rad51_L1 = readRDS(paste0(workdir,"DATA/QUANTIF/Q_WT3_Rad51_L1_RPGC_profmat_readsCounts_GB_SCALED.RDS"))
Q_WT3_Rad51_L2 = readRDS(paste0(workdir,"DATA/QUANTIF/Q_WT3_Rad51_L2_RPGC_profmat_readsCounts_GB_SCALED.RDS"))


Q_SEA4_Rad51_N1 = readRDS(paste0(workdir,"DATA/QUANTIF/Q_SEA4_Rad51_N1_RPGC_profmat_readsCounts_GB_SCALED.RDS"))
Q_SEA4_Rad51_N2 = readRDS(paste0(workdir,"DATA/QUANTIF/Q_SEA4_Rad51_N2_RPGC_profmat_readsCounts_GB_SCALED.RDS"))

Q_SEA4_Rad51_L1 = readRDS(paste0(workdir,"DATA/QUANTIF/Q_SEA4_Rad51_L1_RPGC_profmat_readsCounts_GB_SCALED.RDS"))
Q_SEA4_Rad51_L2 = readRDS(paste0(workdir,"DATA/QUANTIF/Q_SEA4_Rad51_L2_RPGC_profmat_readsCounts_GB_SCALED.RDS"))


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


len = length(Q_H3K36me3_2C4_GB_f)
CTRL_60_65_Q_H3K36me3_2C4_GB_f = names(Q_H3K36me3_2C4_GB_f[(len*0.6) : (len*0.65)])
UP_5_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "top", prct = 5)
UP_25_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "top", prct = 25)
DN_5_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "down", prct = 5)
DN_25_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "down", prct = 25)

#### LIST GENES GROUPS


LIST_5_UP_DN_Q_H3K36me3_2C4=list(
  DN_5_Q_H3K36me3_2C4=DN_5_Q_H3K36me3_2C4_GB_f,
  delta_UP_5_Q_H3K36me3_2C4=UP_5_Q_H3K36me3_2C4_GB_f
)

LIST_25_UP_DN_Q_H3K36me3_2C4=list(
  DN_25_Q_H3K36me3_2C4=DN_25_Q_H3K36me3_2C4_GB_f,
  delta_UP_25_Q_H3K36me3_2C4=UP_25_Q_H3K36me3_2C4_GB_f
)


LIST_Q_H3K36me3_2C4=list(
  DN_5_H3K36me3_2C4_vs_2N4=DN_5_Q_H3K36me3_2C4_GB_f,
  DN_25_Q_H3K36me3_2C4=DN_25_Q_H3K36me3_2C4_GB_f,
  CTRL_60_65_H3K36me3_2C4_vs_2N4 = CTRL_60_65_Q_H3K36me3_2C4_GB_f,
  delta_UP_25_Q_H3K36me3_2C4=UP_25_Q_H3K36me3_2C4_GB_f,
  delta_U_P5_H3K36me3_2C4_vs_2N4 = UP_5_Q_H3K36me3_2C4_GB_f
)


#####################################################################################-
#         PLOT  ----
#####################################################################################-

### REPLICAT 1
# NELF KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_WT3_Rad51_N1, quantifKD = Q_SEA4_Rad51_N1, cond1 = "Q_WT3_Rad51_N1", cond2="Q_SEA4_Rad51_N1", filterGNList = LIST_5_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_WT3_Rad51_N1", Cond = "Q_SEA4_Rad51_N1", select = "LIST_5_UP_DN_Q_H3K36me3_2C4", info = NULL)

Boxplot_wilcoxListFilter_REF(quantifWT = Q_WT3_Rad51_N1, quantifKD = Q_SEA4_Rad51_N1, cond1 = "Q_WT3_Rad51_N1", cond2="Q_SEA4_Rad51_N1", filterGNList = LIST_25_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_WT3_Rad51_N1", Cond = "Q_SEA4_Rad51_N1", select = "LIST_25_UP_DN_Q_H3K36me3_2C4", info = NULL)

Boxplot_wilcoxListFilter_REF(quantifWT = Q_WT3_Rad51_N1, quantifKD = Q_SEA4_Rad51_N1, cond1 = "Q_WT3_Rad51_N1", cond2="Q_SEA4_Rad51_N1", filterGNList = LIST_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_WT3_Rad51_N1", Cond = "Q_SEA4_Rad51_N1", select = "LIST_Q_H3K36me3_2C4", info = NULL)


# LUC KD
Boxplot_wilcoxListFilter_REF(quantifWT = Q_WT3_Rad51_L1, quantifKD = Q_SEA4_Rad51_L1, cond1 = "Q_WT3_Rad51_L1", cond2="Q_SEA4_Rad51_L1", filterGNList = LIST_5_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_WT3_Rad51_L1", Cond = "Q_SEA4_Rad51_L1", select = "LIST_5_UP_DN_Q_H3K36me3_2C4", info = NULL)

Boxplot_wilcoxListFilter_REF(quantifWT = Q_WT3_Rad51_L1, quantifKD = Q_SEA4_Rad51_L1, cond1 = "Q_WT3_Rad51_L1", cond2="Q_SEA4_Rad51_L1", filterGNList = LIST_25_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_WT3_Rad51_L1", Cond = "Q_SEA4_Rad51_L1", select = "LIST_25_UP_DN_Q_H3K36me3_2C4", info = NULL)

Boxplot_wilcoxListFilter_REF(quantifWT = Q_WT3_Rad51_L1, quantifKD = Q_SEA4_Rad51_L1, cond1 = "Q_WT3_Rad51_L1", cond2="Q_SEA4_Rad51_L1", filterGNList = LIST_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_WT3_Rad51_L1", Cond = "Q_SEA4_Rad51_L1", select = "LIST_Q_H3K36me3_2C4", info = NULL)





### REPLICAT 2
# NELF KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_WT3_Rad51_N2, quantifKD = Q_SEA4_Rad51_N2, cond1 = "Q_WT3_Rad51_N2", cond2="Q_SEA4_Rad51_N2", filterGNList = LIST_5_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_WT3_Rad51_N2", Cond = "Q_SEA4_Rad51_N2", select = "LIST_5_UP_DN_Q_H3K36me3_2C4", info = NULL)

Boxplot_wilcoxListFilter_REF(quantifWT = Q_WT3_Rad51_N2, quantifKD = Q_SEA4_Rad51_N2, cond1 = "Q_WT3_Rad51_N2", cond2="Q_SEA4_Rad51_N2", filterGNList = LIST_25_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_WT3_Rad51_N2", Cond = "Q_SEA4_Rad51_N2", select = "LIST_25_UP_DN_Q_H3K36me3_2C4", info = NULL)

Boxplot_wilcoxListFilter_REF(quantifWT = Q_WT3_Rad51_N2, quantifKD = Q_SEA4_Rad51_N2, cond1 = "Q_WT3_Rad51_N2", cond2="Q_SEA4_Rad51_N2", filterGNList = LIST_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_WT3_Rad51_N2", Cond = "Q_SEA4_Rad51_N2", select = "LIST_Q_H3K36me3_2C4", info = NULL)


# LUC KD
Boxplot_wilcoxListFilter_REF(quantifWT = Q_WT3_Rad51_L2, quantifKD = Q_SEA4_Rad51_L2, cond1 = "Q_WT3_Rad51_L2", cond2="Q_SEA4_Rad51_L2", filterGNList = LIST_5_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_WT3_Rad51_L2", Cond = "Q_SEA4_Rad51_L2", select = "LIST_5_UP_DN_Q_H3K36me3_2C4", info = NULL)

Boxplot_wilcoxListFilter_REF(quantifWT = Q_WT3_Rad51_L2, quantifKD = Q_SEA4_Rad51_L2, cond1 = "Q_WT3_Rad51_L2", cond2="Q_SEA4_Rad51_L2", filterGNList = LIST_25_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_WT3_Rad51_L2", Cond = "Q_SEA4_Rad51_L2", select = "LIST_25_UP_DN_Q_H3K36me3_2C4", info = NULL)

Boxplot_wilcoxListFilter_REF(quantifWT = Q_WT3_Rad51_L2, quantifKD = Q_SEA4_Rad51_L2, cond1 = "Q_WT3_Rad51_L2", cond2="Q_SEA4_Rad51_L2", filterGNList = LIST_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_WT3_Rad51_L2", Cond = "Q_SEA4_Rad51_L2", select = "LIST_Q_H3K36me3_2C4", info = NULL)



