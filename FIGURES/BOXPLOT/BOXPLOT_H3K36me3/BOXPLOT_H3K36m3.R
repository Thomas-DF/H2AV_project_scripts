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
library(plyranges)

#####################################################################################-
#         FUNCTIONS  ----
#####################################################################################-

workdir = "~/Bureau/tdefreitas_genobioinfo/PROJET_H2AV_2025/"
source(paste0(workdir,"functionR/Boxplot_wilcoxListFilter_REF.R"))

outfig = paste0(workdir,"FIGURES/BOXPLOT/BOXPLOT_H3K36me3/")

#####################################################################################-
#         DATA  ----
#####################################################################################-


### CHIPSEQ


Q_ChIP_H3K36me3_FSEA4_Luc = readRDS(paste0(workdir,"DATA/QUANTIF/Q_ChIP_H3K36me3_FSEA4_Luc_readsCounts_GB_SCALED.RDS"))
Q_ChIP_H3K36me3_FSEA4_Nelf = readRDS(paste0(workdir,"DATA/QUANTIF/Q_ChIP_H3K36me3_FSEA4_Nelf_readsCounts_GB_SCALED.RDS"))

Q_ChIP_H3K36me3_FWT3_Luc = readRDS(paste0(workdir,"DATA/QUANTIF/Q_ChIP_H3K36me3_FWT3_Luc_readsCounts_GB_SCALED.RDS"))
Q_ChIP_H3K36me3_FWT3_Nelf = readRDS(paste0(workdir,"DATA/QUANTIF/Q_ChIP_H3K36me3_FWT3_Nelf_readsCounts_GB_SCALED.RDS"))

Q_ChIP_H3K36me3_FSEA4_Luc[is.na(Q_ChIP_H3K36me3_FSEA4_Luc)] = 0
Q_ChIP_H3K36me3_FSEA4_Nelf[is.na(Q_ChIP_H3K36me3_FSEA4_Nelf)] = 0

Q_ChIP_H3K36me3_FWT3_Luc[is.na(Q_ChIP_H3K36me3_FWT3_Luc)] = 0
Q_ChIP_H3K36me3_FWT3_Nelf[is.na(Q_ChIP_H3K36me3_FWT3_Nelf)] = 0


#####################################################################################-

##### GENES GROUPS

LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))

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


len = length(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f)
CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4 = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0.6) : (len*0.65)])

UP_5_ZSCORE_H3K36me3_2C4_vs_2N4 = getNameList(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, topdown = "top", prct = 5)

UP_1_ZSCORE_H3K36me3_2C4_vs_2N4 = getNameList(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, topdown = "top", prct = 1)


#### LIST GENES GROUPS

LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4 = list(
  CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4 = CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4,
  delta_UP_5_ZSCORE_H3K36me3_2C4_vs_2N4 = UP_5_ZSCORE_H3K36me3_2C4_vs_2N4
)

LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4 = list(
  CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4 = CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4,
  delta_UP_1_ZSCORE_H3K36me3_2C4_vs_2N4 = UP_1_ZSCORE_H3K36me3_2C4_vs_2N4
)


#####################################################################################-
#         PLOT  ----
#####################################################################################-

#### CHIPSEQ

### LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4

## FSEA4 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_ChIP_H3K36me3_FSEA4_Luc, quantifKD = Q_ChIP_H3K36me3_FSEA4_Nelf, cond1 = "Q_ChIP_H3K36me3_FSEA4_Luc", cond2="Q_ChIP_H3K36me3_FSEA4_Nelf", filterGNList = LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_ChIP_H3K36me3_FSEA4_Luc", Cond = "Q_ChIP_H3K36me3_FSEA4_Nelf", select = "LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## FWT3 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_ChIP_H3K36me3_FWT3_Luc, quantifKD = Q_ChIP_H3K36me3_FWT3_Nelf, cond1 = "Q_ChIP_H3K36me3_FWT3_Luc", cond2="Q_ChIP_H3K36me3_FWT3_Nelf", filterGNList = LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_ChIP_H3K36me3_FWT3_Luc", Cond = "Q_ChIP_H3K36me3_FWT3_Nelf", select = "LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## Nelf : FWT3 / FSEA4 

Boxplot_wilcoxListFilter_REF(quantifWT = Q_ChIP_H3K36me3_FWT3_Nelf, quantifKD = Q_ChIP_H3K36me3_FSEA4_Nelf, cond1 = "Q_ChIP_H3K36me3_FWT3_Nelf", cond2="Q_ChIP_H3K36me3_FSEA4_Nelf", filterGNList = LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_ChIP_H3K36me3_FWT3_Nelf", Cond = "Q_ChIP_H3K36me3_FSEA4_Nelf", select = "LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## Luc :  FWT3 / FSEA4 

Boxplot_wilcoxListFilter_REF(quantifWT = Q_ChIP_H3K36me3_FWT3_Luc, quantifKD = Q_ChIP_H3K36me3_FSEA4_Luc, cond1 = "Q_ChIP_H3K36me3_FWT3_Luc", cond2="Q_ChIP_H3K36me3_FSEA4_Luc", filterGNList = LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_ChIP_H3K36me3_FWT3_Luc", Cond = "Q_ChIP_H3K36me3_FSEA4_Luc", select = "LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)



#### LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4

## FSEA4 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_ChIP_H3K36me3_FSEA4_Luc, quantifKD = Q_ChIP_H3K36me3_FSEA4_Nelf, cond1 = "Q_ChIP_H3K36me3_FSEA4_Luc", cond2="Q_ChIP_H3K36me3_FSEA4_Nelf", filterGNList = LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_ChIP_H3K36me3_FSEA4_Luc", Cond = "Q_ChIP_H3K36me3_FSEA4_Nelf", select = "LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## FWT3 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_ChIP_H3K36me3_FWT3_Luc, quantifKD = Q_ChIP_H3K36me3_FWT3_Nelf, cond1 = "Q_ChIP_H3K36me3_FWT3_Luc", cond2="Q_ChIP_H3K36me3_FWT3_Nelf", filterGNList = LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_ChIP_H3K36me3_FWT3_Luc", Cond = "Q_ChIP_H3K36me3_FWT3_Nelf", select = "LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)


## Nelf : FWT3 / FSEA4 

Boxplot_wilcoxListFilter_REF(quantifWT = Q_ChIP_H3K36me3_FWT3_Nelf, quantifKD = Q_ChIP_H3K36me3_FSEA4_Nelf, cond1 = "Q_ChIP_H3K36me3_FWT3_Nelf", cond2="Q_ChIP_H3K36me3_FSEA4_Nelf", filterGNList = LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_ChIP_H3K36me3_FWT3_Nelf", Cond = "Q_ChIP_H3K36me3_FSEA4_Nelf", select = "LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## Luc :  FWT3 / FSEA4 

Boxplot_wilcoxListFilter_REF(quantifWT = Q_ChIP_H3K36me3_FWT3_Luc, quantifKD = Q_ChIP_H3K36me3_FSEA4_Luc, cond1 = "Q_ChIP_H3K36me3_FWT3_Luc", cond2="Q_ChIP_H3K36me3_FSEA4_Luc", filterGNList = LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_ChIP_H3K36me3_FWT3_Luc", Cond = "Q_ChIP_H3K36me3_FSEA4_Luc", select = "LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)


