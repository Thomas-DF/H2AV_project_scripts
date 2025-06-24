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

outfig = paste0(workdir,"FIGURES/BOXPLOT/BOXPLOT_DRIP/")

#####################################################################################-
#         DATA  ----
#####################################################################################-

### DRIPSEQ

Q_DRIP_HLuc = readRDS(paste0(workdir,"DATA/QUANTIF/Q_DRIP_HLuc_profmat_readsCounts_GB_SCALED.RDS"))
Q_DRIP_HLuc_RNAse = readRDS(paste0(workdir,"DATA/QUANTIF/Q_DRIP_HLuc_RNAse_profmat_readsCounts_GB_SCALED.RDS"))

Q_DRIP_Nelf = readRDS(paste0(workdir,"DATA/QUANTIF/Q_DRIP_Nelf_profmat_readsCounts_GB_SCALED.RDS"))
Q_DRIP_Nelf_RNAse = readRDS(paste0(workdir,"DATA/QUANTIF/Q_DRIP_Nelf_RNAse_profmat_readsCounts_GB_SCALED.RDS"))

Q_DRIP_HypB = readRDS(paste0(workdir,"DATA/QUANTIF/Q_DRIP_HypB_profmat_readsCounts_GB_SCALED.RDS"))
Q_DRIP_HypB_RNAse = readRDS(paste0(workdir,"DATA/QUANTIF/Q_DRIP_HypB_RNAse_profmat_readsCounts_GB_SCALED.RDS"))

Q_DRIP_HLuc[is.na(Q_DRIP_HLuc)] = 0
Q_DRIP_HLuc_RNAse[is.na(Q_DRIP_HLuc_RNAse)] = 0

Q_DRIP_Nelf[is.na(Q_DRIP_Nelf)] = 0
Q_DRIP_Nelf_RNAse[is.na(Q_DRIP_Nelf_RNAse)] = 0

Q_DRIP_HypB[is.na(Q_DRIP_HypB)] = 0
Q_DRIP_HypB_RNAse[is.na(Q_DRIP_HypB_RNAse)] = 0


ZSCORE_PROFMAT_Q_DRIP_HLuc_vs_HLuc_RNAse = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_Q_DRIP_HLuc_vs_HLuc_RNAse.RDS"))
ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse =  readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse.RDS"))

ZSCORE_PROFMAT_Q_DRIP_HLuc_vs_HLuc_RNAse[is.na(ZSCORE_PROFMAT_Q_DRIP_HLuc_vs_HLuc_RNAse)] = 0
ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse[is.na(ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse)] = 0


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

DOWN_5_ZSCORE_H3K36me3_2C4_vs_2N4 = getNameList(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, topdown = "down", prct = 5)

#### LIST GENES GROUPS

LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4 = list(
  CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4 = CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4,
  delta_UP_5_ZSCORE_H3K36me3_2C4_vs_2N4 = UP_5_ZSCORE_H3K36me3_2C4_vs_2N4
)

LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4 = list(
  CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4 = CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4,
  delta_UP_1_ZSCORE_H3K36me3_2C4_vs_2N4 = UP_1_ZSCORE_H3K36me3_2C4_vs_2N4
)

LIST_5_UP_DOWN_ZSCORE_H3K36me3_2C4_vs_2N4 = list(
  DOWN_5_ZSCORE_H3K36me3_2C4_vs_2N4 = DOWN_5_ZSCORE_H3K36me3_2C4_vs_2N4,
  delta_UP_5_ZSCORE_H3K36me3_2C4_vs_2N4 = UP_5_ZSCORE_H3K36me3_2C4_vs_2N4
)
#####################################################################################-
#         PLOT  ----
#####################################################################################-
#### DRIPSEQ

### LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4

## HLuc /  HLuc RNAse

Boxplot_wilcoxListFilter_REF(quantifWT = Q_DRIP_HLuc, quantifKD = Q_DRIP_HLuc_RNAse, cond1 = "Q_DRIP_HLuc", cond2="Q_DRIP_HLuc_RNAse", filterGNList = LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_DRIP_HLuc", Cond = "Q_DRIP_HLuc_RNAse", select = "LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## Nelf / Nelf RNAse

Boxplot_wilcoxListFilter_REF(quantifWT = Q_DRIP_Nelf, quantifKD = Q_DRIP_Nelf_RNAse, cond1 = "Q_DRIP_Nelf", cond2="Q_DRIP_Nelf_RNAse", filterGNList = LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_DRIP_Nelf", Cond = "Q_DRIP_Nelf_RNAse", select = "LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## HypB / HypB RNAse

Boxplot_wilcoxListFilter_REF(quantifWT = Q_DRIP_HypB, quantifKD = Q_DRIP_HypB_RNAse, cond1 = "Q_DRIP_HypB", cond2="Q_DRIP_HypB_RNAse", filterGNList = LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_DRIP_HypB", Cond = "Q_DRIP_HypB_RNAse", select = "LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## ZSCORE HLuc_vs_HLuc_RNAse / ZSCORE Nelf_vs_Nelf_RNAse

Boxplot_wilcoxListFilter_REF(quantifWT = ZSCORE_PROFMAT_Q_DRIP_HLuc_vs_HLuc_RNAse, quantifKD = ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse, cond1 = "ZSCORE_PROFMAT_Q_DRIP_HLuc_vs_HLuc_RNAse", cond2="ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse", filterGNList = LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=F,
                             outdir = outfig, readQuantif = "ZSCORE_PROFMAT_Q_DRIP_HLuc_vs_HLuc_RNAse", Cond = "ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse", select = "LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)



### LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4

## HLuc /  HLuc RNAse

Boxplot_wilcoxListFilter_REF(quantifWT = Q_DRIP_HLuc, quantifKD = Q_DRIP_HLuc_RNAse, cond1 = "Q_DRIP_HLuc", cond2="Q_DRIP_HLuc_RNAse", filterGNList = LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_DRIP_HLuc", Cond = "Q_DRIP_HLuc_RNAse", select = "LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## Nelf / Nelf RNAse

Boxplot_wilcoxListFilter_REF(quantifWT = Q_DRIP_Nelf, quantifKD = Q_DRIP_Nelf_RNAse, cond1 = "Q_DRIP_Nelf", cond2="Q_DRIP_Nelf_RNAse", filterGNList = LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_DRIP_Nelf", Cond = "Q_DRIP_Nelf_RNAse", select = "LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## HypB / HypB RNAse

Boxplot_wilcoxListFilter_REF(quantifWT = Q_DRIP_HypB, quantifKD = Q_DRIP_HypB_RNAse, cond1 = "Q_DRIP_HypB", cond2="Q_DRIP_HypB_RNAse", filterGNList = LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_DRIP_HypB", Cond = "Q_DRIP_HypB_RNAse", select = "LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)


## ZSCORE HLuc_vs_HLuc_RNAse / ZSCORE Nelf_vs_Nelf_RNAse

Boxplot_wilcoxListFilter_REF(quantifWT = ZSCORE_PROFMAT_Q_DRIP_HLuc_vs_HLuc_RNAse, quantifKD = ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse, cond1 = "ZSCORE_PROFMAT_Q_DRIP_HLuc_vs_HLuc_RNAse", cond2="ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse", filterGNList = LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=F,
                             outdir = outfig, readQuantif = "ZSCORE_PROFMAT_Q_DRIP_HLuc_vs_HLuc_RNAse", Cond = "ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse", select = "LIST_1_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)



#### LIST_5_UP_DOWN_ZSCORE_H3K36me3_2C4_vs_2N4

Boxplot_wilcoxListFilter_REF(quantifWT = ZSCORE_PROFMAT_Q_DRIP_HLuc_vs_HLuc_RNAse, quantifKD = ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse, cond1 = "ZSCORE_PROFMAT_Q_DRIP_HLuc_vs_HLuc_RNAse", cond2="ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse", filterGNList = LIST_5_UP_DOWN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=F,
                             outdir = outfig, readQuantif = "ZSCORE_PROFMAT_Q_DRIP_HLuc_vs_HLuc_RNAse", Cond = "ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse", select = "LIST_5_UP_DOWN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)


## HLuc /  HLuc RNAse

Boxplot_wilcoxListFilter_REF(quantifWT = Q_DRIP_HLuc, quantifKD = Q_DRIP_HLuc_RNAse, cond1 = "Q_DRIP_HLuc", cond2="Q_DRIP_HLuc_RNAse", filterGNList = LIST_5_UP_DOWN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_DRIP_HLuc", Cond = "Q_DRIP_HLuc_RNAse", select = "LIST_5_UP_DOWN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## Nelf / Nelf RNAse

Boxplot_wilcoxListFilter_REF(quantifWT = Q_DRIP_Nelf, quantifKD = Q_DRIP_Nelf_RNAse, cond1 = "Q_DRIP_Nelf", cond2="Q_DRIP_Nelf_RNAse", filterGNList = LIST_5_UP_DOWN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_DRIP_Nelf", Cond = "Q_DRIP_Nelf_RNAse", select = "LIST_5_UP_DOWN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## HypB / HypB RNAse

Boxplot_wilcoxListFilter_REF(quantifWT = Q_DRIP_HypB, quantifKD = Q_DRIP_HypB_RNAse, cond1 = "Q_DRIP_HypB", cond2="Q_DRIP_HypB_RNAse", filterGNList = LIST_5_UP_DOWN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_DRIP_HypB", Cond = "Q_DRIP_HypB_RNAse", select = "LIST_5_UP_DOWN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)
