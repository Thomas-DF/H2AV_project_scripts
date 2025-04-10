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

outfig = paste0(workdir,"FIGURES/BOXPLOT/BOXPLOT_Q_H2AV_SEA4_WT3/BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f/")
outfig2 = paste0(workdir,"FIGURES/BOXPLOT/BOXPLOT_Q_H2AV_SEA4_WT3/BY_Q_H3K36me3_2C4_GB_f/")

#####################################################################################-
#         DATA  ----
#####################################################################################-

### REPLICATS 3

Q_GB_SEA4_NELFKD_R3_f = readRDS(paste0(workdir,"DATA/QUANTIF/Q_GB_SEA4_NELFKD_R3_f.RDS"))
Q_GB_SEA4_WT_R3_f = readRDS(paste0(workdir,"DATA/QUANTIF/Q_GB_SEA4_WT_R3_f.RDS"))

Q_GB_WT3_NELFKD_R3_f = readRDS(paste0(workdir,"DATA/QUANTIF/Q_GB_WT3_NELFKD_R3_f.RDS"))
Q_GB_WT3_WT_R3_f = readRDS(paste0(workdir,"DATA/QUANTIF/Q_GB_WT3_WT_R3_f.RDS"))


### REPLICATS 2

Q_H2AV_SEA4_LucKD=readRDS(paste0(workdir,"DATA/QUANTIF/Q_H2AV-R2_SEA4_Luc-KD_RPKM_upstr500_dnstr500.RDS"))
Q_H2AV_SEA4_NelfKD=readRDS(paste0(workdir,"DATA/QUANTIF/Q_H2AV-R2_SEA4_Nelf-KD_RPKM_upstr500_dnstr500.RDS"))


Q_H2AV_WT3_LucKD=readRDS(paste0(workdir,"DATA/QUANTIF/Q_H2AV-R2_WT3_Luc-KD_RPKM_upstr500_dnstr500.RDS"))
Q_H2AV_WT3_NelfKD=readRDS(paste0(workdir,"DATA/QUANTIF/Q_H2AV-R2_WT3_Nelf-KD_RPKM_upstr500_dnstr500.RDS"))


### FORMAT ADAPTATION

names(Q_H2AV_SEA4_LucKD) <- paste0(names(Q_H2AV_SEA4_LucKD), ".1")
names(Q_H2AV_SEA4_NelfKD) <- paste0(names(Q_H2AV_SEA4_NelfKD), ".1")

names(Q_H2AV_WT3_LucKD) <- paste0(names(Q_H2AV_WT3_LucKD), ".1")
names(Q_H2AV_WT3_NelfKD) <- paste0(names(Q_H2AV_WT3_NelfKD), ".1")



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



len = length(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f)

CTRL5_ZSCORE_H3K36me3_2C4_vs_2N4 = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0.6) : (len*0.65)])
DN5_ZSCORE_H3K36me3_2C4_vs_2N4 = getNameList(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, topdown = "down", prct = 5)
UP5_ZSCORE_H3K36me3_2C4_vs_2N4 = getNameList(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, topdown = "top", prct = 5)
DN1_ZSCORE_H3K36me3_2C4_vs_2N4 = getNameList(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, topdown = "down", prct = 1)
UP1_ZSCORE_H3K36me3_2C4_vs_2N4 = getNameList(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, topdown = "top", prct = 1)

len = length(Q_H3K36me3_2C4_GB_f)
CTRL_60_65_Q_H3K36me3_2C4_GB_f = names(Q_H3K36me3_2C4_GB_f[(len*0.6) : (len*0.65)])
UP_5_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "top", prct = 5)
UP_25_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "top", prct = 25)
DN_5_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "down", prct = 5)
DN_25_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "down", prct = 25)

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

###


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
########## ZSCORE_H3K36me3_2C4_vs_2N4 ########## 

### REPLICATS 2
### TOP / BOT 5%

## SEA4 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_SEA4_LucKD, quantifKD = Q_H2AV_SEA4_NelfKD, cond1 = "Q_H2AV_SEA4_LucKD", cond2="Q_H2AV_SEA4_NelfKD", filterGNList = LIST_5_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_H2AV_SEA4_LucKD", Cond = "Q_H2AV_SEA4_NelfKD", select = "LIST_5_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## WT3 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_WT3_LucKD, quantifKD = Q_H2AV_WT3_NelfKD, cond1 = "Q_H2AV_WT3_LucKD", cond2="Q_H2AV_WT3_NelfKD", filterGNList = LIST_5_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_H2AV_WT3_LucKD", Cond = "Q_H2AV_WT3_NelfKD", select = "LIST_5_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)


### TOP / BOT 1%

## SEA4 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_SEA4_LucKD, quantifKD = Q_H2AV_SEA4_NelfKD, cond1 = "Q_H2AV_SEA4_LucKD", cond2="Q_H2AV_SEA4_NelfKD", filterGNList = LIST_1_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_H2AV_SEA4_LucKD", Cond = "Q_H2AV_SEA4_NelfKD", select = "LIST_1_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## WT3 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_WT3_LucKD, quantifKD = Q_H2AV_WT3_NelfKD, cond1 = "Q_H2AV_WT3_LucKD", cond2="Q_H2AV_WT3_NelfKD", filterGNList = LIST_1_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_H2AV_WT3_LucKD", Cond = "Q_H2AV_WT3_NelfKD", select = "LIST_1_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)



### TOP / BOT / CTRL 60-65%

## SEA4 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_SEA4_LucKD, quantifKD = Q_H2AV_SEA4_NelfKD, cond1 = "Q_H2AV_SEA4_LucKD", cond2="Q_H2AV_SEA4_NelfKD", filterGNList = LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_H2AV_SEA4_LucKD", Cond = "Q_H2AV_SEA4_NelfKD", select = "LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## WT3 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_WT3_LucKD, quantifKD = Q_H2AV_WT3_NelfKD, cond1 = "Q_H2AV_WT3_LucKD", cond2="Q_H2AV_WT3_NelfKD", filterGNList = LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_H2AV_WT3_LucKD", Cond = "Q_H2AV_WT3_NelfKD", select = "LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)






#### REPLICATS 3
### TOP / BOT 5%

## SEA4 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_GB_SEA4_WT_R3_f, quantifKD = Q_GB_SEA4_NELFKD_R3_f, cond1 = "Q_GB_SEA4_WT_R3_f", cond2="Q_GB_SEA4_NELFKD_R3_f", filterGNList = LIST_5_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_GB_SEA4_WT_R3_f", Cond = "Q_GB_SEA4_NELFKD_R3_f", select = "LIST_5_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## WT3 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_GB_WT3_WT_R3_f, quantifKD = Q_GB_WT3_NELFKD_R3_f, cond1 = "Q_GB_WT3_WT_R3_f", cond2="Q_GB_WT3_NELFKD_R3_f", filterGNList = LIST_5_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_GB_WT3_WT_R3_f", Cond = "Q_GB_WT3_NELFKD_R3_f", select = "LIST_5_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

### TOP / BOT 1%

## SEA4 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_GB_SEA4_WT_R3_f, quantifKD = Q_GB_SEA4_NELFKD_R3_f, cond1 = "Q_GB_SEA4_WT_R3_f", cond2="Q_GB_SEA4_NELFKD_R3_f", filterGNList = LIST_1_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_GB_SEA4_WT_R3_f", Cond = "Q_GB_SEA4_NELFKD_R3_f", select = "LIST_1_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## WT3 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_GB_WT3_WT_R3_f, quantifKD = Q_GB_WT3_NELFKD_R3_f, cond1 = "Q_GB_WT3_WT_R3_f", cond2="Q_GB_WT3_NELFKD_R3_f", filterGNList = LIST_1_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_GB_WT3_WT_R3_f", Cond = "Q_GB_WT3_NELFKD_R3_f", select = "LIST_1_UP_DN_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

### TOP / BOT / CTRL 60-65%

## SEA4 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_GB_SEA4_WT_R3_f, quantifKD = Q_GB_SEA4_NELFKD_R3_f, cond1 = "Q_GB_SEA4_WT_R3_f", cond2="Q_GB_SEA4_NELFKD_R3_f", filterGNList = LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_GB_SEA4_WT_R3_f", Cond = "Q_GB_SEA4_NELFKD_R3_f", select = "LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)

## WT3 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_GB_WT3_WT_R3_f, quantifKD = Q_GB_WT3_NELFKD_R3_f, cond1 = "Q_GB_WT3_WT_R3_f", cond2="Q_GB_WT3_NELFKD_R3_f", filterGNList = LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "Q_GB_WT3_WT_R3_f", Cond = "Q_GB_WT3_NELFKD_R3_f", select = "LIST_5_UP_CTRL_ZSCORE_H3K36me3_2C4_vs_2N4", info = NULL)



########## Q_H3K36me3_2C4_GB_f ########## 


### REPLICATS 2
### TOP / BOT 5%

## SEA4 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_SEA4_LucKD, quantifKD = Q_H2AV_SEA4_NelfKD, cond1 = "Q_H2AV_SEA4_LucKD", cond2="Q_H2AV_SEA4_NelfKD", filterGNList = LIST_5_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig2, readQuantif = "Q_H2AV_SEA4_LucKD", Cond = "Q_H2AV_SEA4_NelfKD", select = "LIST_5_UP_DN_Q_H3K36me3_2C4", info = NULL)

## WT3 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_WT3_LucKD, quantifKD = Q_H2AV_WT3_NelfKD, cond1 = "Q_H2AV_WT3_LucKD", cond2="Q_H2AV_WT3_NelfKD", filterGNList = LIST_5_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig2, readQuantif = "Q_H2AV_WT3_LucKD", Cond = "Q_H2AV_WT3_NelfKD", select = "LIST_5_UP_DN_Q_H3K36me3_2C4", info = NULL)


### TOP / BOT 25%

## SEA4 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_SEA4_LucKD, quantifKD = Q_H2AV_SEA4_NelfKD, cond1 = "Q_H2AV_SEA4_LucKD", cond2="Q_H2AV_SEA4_NelfKD", filterGNList = LIST_25_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig2, readQuantif = "Q_H2AV_SEA4_LucKD", Cond = "Q_H2AV_SEA4_NelfKD", select = "LIST_25_UP_DN_Q_H3K36me3_2C4", info = NULL)

## WT3 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_WT3_LucKD, quantifKD = Q_H2AV_WT3_NelfKD, cond1 = "Q_H2AV_WT3_LucKD", cond2="Q_H2AV_WT3_NelfKD", filterGNList = LIST_25_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig2, readQuantif = "Q_H2AV_WT3_LucKD", Cond = "Q_H2AV_WT3_NelfKD", select = "LIST_25_UP_DN_Q_H3K36me3_2C4", info = NULL)



### TOP / BOT / CTRL 60-65%

## SEA4 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_SEA4_LucKD, quantifKD = Q_H2AV_SEA4_NelfKD, cond1 = "Q_H2AV_SEA4_LucKD", cond2="Q_H2AV_SEA4_NelfKD", filterGNList = LIST_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig2, readQuantif = "Q_H2AV_SEA4_LucKD", Cond = "Q_H2AV_SEA4_NelfKD", select = "LIST_Q_H3K36me3_2C4", info = NULL)

## WT3 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_H2AV_WT3_LucKD, quantifKD = Q_H2AV_WT3_NelfKD, cond1 = "Q_H2AV_WT3_LucKD", cond2="Q_H2AV_WT3_NelfKD", filterGNList = LIST_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig2, readQuantif = "Q_H2AV_WT3_LucKD", Cond = "Q_H2AV_WT3_NelfKD", select = "LIST_Q_H3K36me3_2C4", info = NULL)






#### REPLICATS 3
### TOP / BOT 5%

## SEA4 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_GB_SEA4_WT_R3_f, quantifKD = Q_GB_SEA4_NELFKD_R3_f, cond1 = "Q_GB_SEA4_WT_R3_f", cond2="Q_GB_SEA4_NELFKD_R3_f", filterGNList = LIST_5_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig2, readQuantif = "Q_GB_SEA4_WT_R3_f", Cond = "Q_GB_SEA4_NELFKD_R3_f", select = "LIST_5_UP_DN_Q_H3K36me3_2C4", info = NULL)

## WT3 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_GB_WT3_WT_R3_f, quantifKD = Q_GB_WT3_NELFKD_R3_f, cond1 = "Q_GB_WT3_WT_R3_f", cond2="Q_GB_WT3_NELFKD_R3_f", filterGNList = LIST_5_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig2, readQuantif = "Q_GB_WT3_WT_R3_f", Cond = "Q_GB_WT3_NELFKD_R3_f", select = "LIST_5_UP_DN_Q_H3K36me3_2C4", info = NULL)

### TOP / BOT 1%

## SEA4 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_GB_SEA4_WT_R3_f, quantifKD = Q_GB_SEA4_NELFKD_R3_f, cond1 = "Q_GB_SEA4_WT_R3_f", cond2="Q_GB_SEA4_NELFKD_R3_f", filterGNList = LIST_25_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig2, readQuantif = "Q_GB_SEA4_WT_R3_f", Cond = "Q_GB_SEA4_NELFKD_R3_f", select = "LIST_25_UP_DN_Q_H3K36me3_2C4", info = NULL)

## WT3 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_GB_WT3_WT_R3_f, quantifKD = Q_GB_WT3_NELFKD_R3_f, cond1 = "Q_GB_WT3_WT_R3_f", cond2="Q_GB_WT3_NELFKD_R3_f", filterGNList = LIST_25_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig2, readQuantif = "Q_GB_WT3_WT_R3_f", Cond = "Q_GB_WT3_NELFKD_R3_f", select = "LIST_25_UP_DN_Q_H3K36me3_2C4", info = NULL)

### TOP / BOT / CTRL 60-65%

## SEA4 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_GB_SEA4_WT_R3_f, quantifKD = Q_GB_SEA4_NELFKD_R3_f, cond1 = "Q_GB_SEA4_WT_R3_f", cond2="Q_GB_SEA4_NELFKD_R3_f", filterGNList = LIST_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig2, readQuantif = "Q_GB_SEA4_WT_R3_f", Cond = "Q_GB_SEA4_NELFKD_R3_f", select = "LIST_Q_H3K36me3_2C4", info = NULL)

## WT3 : Luc KD / Nelf KD

Boxplot_wilcoxListFilter_REF(quantifWT = Q_GB_WT3_WT_R3_f, quantifKD = Q_GB_WT3_NELFKD_R3_f, cond1 = "Q_GB_WT3_WT_R3_f", cond2="Q_GB_WT3_NELFKD_R3_f", filterGNList = LIST_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig2, readQuantif = "Q_GB_WT3_WT_R3_f", Cond = "Q_GB_WT3_NELFKD_R3_f", select = "LIST_Q_H3K36me3_2C4", info = NULL)



