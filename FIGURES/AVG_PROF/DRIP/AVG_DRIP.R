#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(seqplots)

#####################################################################################-
#         FUNCTIONS  ----
#####################################################################################-

workdir = "~/Bureau/tdefreitas_genobioinfo/PROJET_H2AV_2025/"

source(paste0(workdir, "functionR/AVG_PROFILE.R"))

tmp <- paste0(workdir,"FIGURES/AVG_PROF/TMPgetPlotSetArray")

#####################################################################################-
#         DATA  ----
#####################################################################################-

### GENOME REF 

r6_ref_genes = readRDS(paste0(workdir,"DATA/r6.13/TxDb.GR.dm6.RDS"))



### DRIP

DRIP_HLuc = paste0(workdir,"DATA/DRIPseq/DRIP_HLuc.bigWig")
DRIP_HLuc_RNAse = paste0(workdir,"DATA/DRIPseq/DRIP_HLuc_RNAse.bigWig")

DRIP_Nelf = paste0(workdir,"DATA/DRIPseq/DRIP_Nelf.bigWig")
DRIP_Nelf_RNAse = paste0(workdir,"DATA/DRIPseq/DRIP_Nelf_RNAse.bigWig")

#####################################################################################-

## GENE GROUPES 

LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f=LIST_QUANTIF_K36$ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f

Q_H3K36me3_2C4_GB_f=LIST_QUANTIF_K36$Q_H3K36me3_2C4_GB_f

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

UP_5_ZSCORE_H3K36me3_2C4_vs_2N4 = getNameList(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, topdown = "top", prct = 5)
UP_1_ZSCORE_H3K36me3_2C4_vs_2N4 = getNameList(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, topdown = "top", prct = 1)
len = length(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f)
CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4 = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0.6) : (len*0.65)])
DN_5_ZSCORE_H3K36me3_2C4_vs_2N4 = getNameList(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, topdown = "down", prct = 5)

UP_5_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "top", prct = 5)
UP_1_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "top", prct = 1)
len = length(Q_H3K36me3_2C4_GB_f)
CTRL_60_65_Q_H3K36me3_2C4_GB_f = names(Q_H3K36me3_2C4_GB_f[(len*0.6) : (len*0.65)])
DN_5_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "down", prct = 5)


get_GR_feat = function(refGN, GNlist){
  myovlp=refGN[refGN$name %in% GNlist]
  return(myovlp)
}

GR_UP_5_ZSCORE_H3K36me3_2C4_vs_2N4=get_GR_feat(r6_ref_genes,UP_5_ZSCORE_H3K36me3_2C4_vs_2N4)
GR_UP_1_ZSCORE_H3K36me3_2C4_vs_2N4=get_GR_feat(r6_ref_genes,UP_1_ZSCORE_H3K36me3_2C4_vs_2N4)
GR_CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4=get_GR_feat(r6_ref_genes,CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4)
GR_DN_5_ZSCORE_H3K36me3_2C4_vs_2N4=get_GR_feat(r6_ref_genes,DN_5_ZSCORE_H3K36me3_2C4_vs_2N4)

GR_list_toPlot = c("GR_UP_5_ZSCORE_H3K36me3_2C4_vs_2N4", "GR_UP_1_ZSCORE_H3K36me3_2C4_vs_2N4", "GR_CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4", "GR_DN_5_ZSCORE_H3K36me3_2C4_vs_2N4")

GR_UP_5_Q_H3K36me3_2C4_GB_f=get_GR_feat(r6_ref_genes,UP_5_Q_H3K36me3_2C4_GB_f)
GR_UP_1_Q_H3K36me3_2C4_GB_f=get_GR_feat(r6_ref_genes,UP_1_Q_H3K36me3_2C4_GB_f)
GR_CTRL_60_65_Q_H3K36me3_2C4_GB_f=get_GR_feat(r6_ref_genes,CTRL_60_65_Q_H3K36me3_2C4_GB_f)
GR_DN_5_Q_H3K36me3_2C4_GB_f=get_GR_feat(r6_ref_genes,DN_5_Q_H3K36me3_2C4_GB_f)

GR_list_toPlot_H3K36me3 = c("GR_UP_5_Q_H3K36me3_2C4_GB_f", "GR_UP_1_Q_H3K36me3_2C4_GB_f", "GR_CTRL_60_65_Q_H3K36me3_2C4_GB_f", "GR_DN_5_Q_H3K36me3_2C4_GB_f")

#####################################################################################-
#         PLOT  ----
#####################################################################################-

### ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f

## DRIP HLuc - AVERAGE PLOT

pdf(paste0(workdir,"FIGURES/AVG_PROF/DRIP/AVG_PROF_DRIP_HLuc_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
for(GR in GR_list_toPlot){
  seqPlotSDoutliers_scaleFact(c(DRIP_HLuc,DRIP_HLuc_RNAse),tmp,GR,c(0,500),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
  }
dev.off()


## DRIP Nelf - AVERAGE PLOT

pdf(paste0(workdir,"FIGURES/AVG_PROF/DRIP/AVG_PROF_DRIP_Nelf_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
for(GR in GR_list_toPlot){
  seqPlotSDoutliers_scaleFact(c(DRIP_Nelf,DRIP_Nelf_RNAse),tmp,GR,c(0,500),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()



### Q_H3K36me3_2C4_GB_f

## DRIP HLuc - AVERAGE PLOT

pdf(paste0(workdir,"FIGURES/AVG_PROF/DRIP/AVG_PROF_DRIP_HLuc_Q_H3K36me3_2C4_GB_f.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(DRIP_HLuc,DRIP_HLuc_RNAse),tmp,GR,c(0,700),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()


## DRIP Nelf - AVERAGE PLOT

pdf(paste0(workdir,"FIGURES/AVG_PROF/DRIP/AVG_PROF_DRIP_Nelf_Q_H3K36me3_2C4_GB_f.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(DRIP_Nelf,DRIP_Nelf_RNAse),tmp,GR,c(0,700),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()
