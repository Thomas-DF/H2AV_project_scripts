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
seqlevels(r6_ref_genes) <- gsub("chr", "", seqlevels(r6_ref_genes))


### BIGWIG FILES
H2AV_PN_2 = paste0(workdir,"DATA/CHIPSEQ/H2AV_2022/PN_trimmed_filt_sort_RPGC.bw")
H2AV_PN_1 = paste0(workdir,"DATA/CHIPSEQ/H2AV_2022/H2AV_PN_J_L1_RPGC.bw")


H2AV_PW_1 = paste0(workdir,"DATA/CHIPSEQ/H2AV_2022/H2AV_PW_1_L1_RPGC.bw")
H2AV_PW_2 = paste0(workdir,"DATA/CHIPSEQ/H2AV_2022/H2AV_PW_2_L1_RPGC.bw")

H2AV_PHYPB_R1=paste0(workdir,"DATA/CHIPSEQ/H2AV_2022/PHYP_B_L1_trimmed_filt_sort_RPGC.bw")
H2AV_PHYPB_R2=paste0(workdir,"DATA/CHIPSEQ/H2AV_2022/H2AV_PHYPB_A_L1_RPGC.bw")

IHYP_B_L1=paste0(workdir,"DATA/CHIPSEQ/H2AV_2019/IHYP_B_L1_trimmed_filt_sort_RPGC.bw")
IHYP_B_L2=paste0(workdir,"DATA/CHIPSEQ/H2AV_2019/IHYP_B_L2_trimmed_filt_sort_RPGC.bw")

IHYPH_B_L1=paste0(workdir,"DATA/CHIPSEQ/H2AV_2019/IHYPH_B_L1_trimmed_filt_sort_RPGC.bw")
IHYPH_B_L2=paste0(workdir,"DATA/CHIPSEQ/H2AV_2019/IHYPH_B_L2_trimmed_filt_sort_RPGC.bw")

PHYP_B_L1=paste0(workdir,"DATA/CHIPSEQ/H2AV_2019/PHYP_B_L1_trimmed_filt_sort_RPGC.bw")
PHYP_B_L2=paste0(workdir,"DATA/CHIPSEQ/H2AV_2019/PHYP_B_L2_trimmed_filt_sort_RPGC.bw")

PHYP_B_L1_inputnorm=paste0(workdir,"DATA/CHIPSEQ/H2AV_2019/PHYP_B_L1_trimmed_filt_sort_inputnormRPGC.bw")
PHYP_B_L2_inputnorm=paste0(workdir,"DATA/CHIPSEQ/H2AV_2019/PHYP_B_L2_trimmed_filt_sort_inputnormRPGC.bw")

H2AV_R3_Luc_KD = paste0(workdir,"DATA/CHIPSEQ/H2AV-R3_Luc-KD_CPM.bw")
H2AV_R3_Luc_HU_KD = paste0(workdir,"DATA/CHIPSEQ/H2AV-R3_Luc-KD-HU_CPM.bw")

H2AV_R3_HypB_KD = paste0(workdir,"DATA/CHIPSEQ/H2AV-R3_HypB-KD_CPM.bw")
H2AV_R3_HypB_HU_KD = paste0(workdir,"DATA/CHIPSEQ/H2AV-R3_HypB-KD-HU_CPM.bw")


#####################################################################################-

## GENE GROUPES 

LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
Q_H3K36me3_2C4_GB_f = LIST_QUANTIF_K36$Q_H3K36me3_2C4_GB_f

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


UP_5_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "top", prct = 5)
UP_1_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "top", prct = 1)
len = length(Q_H3K36me3_2C4_GB_f)
CTRL_60_65_Q_H3K36me3_2C4_GB_f = names(Q_H3K36me3_2C4_GB_f[(len*0.6) : (len*0.65)])
DN_5_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "down", prct = 5)

get_GR_feat = function(refGN, GNlist){
  myovlp=refGN[refGN$name %in% GNlist]
  return(myovlp)
}


GR_UP_5_Q_H3K36me3_2C4_GB_f=get_GR_feat(r6_ref_genes,UP_5_Q_H3K36me3_2C4_GB_f)
GR_UP_1_Q_H3K36me3_2C4_GB_f=get_GR_feat(r6_ref_genes,UP_1_Q_H3K36me3_2C4_GB_f)
GR_CTRL_60_65_Q_H3K36me3_2C4_GB_f=get_GR_feat(r6_ref_genes,CTRL_60_65_Q_H3K36me3_2C4_GB_f)
GR_DN_5_Q_H3K36me3_2C4_GB_f=get_GR_feat(r6_ref_genes,DN_5_Q_H3K36me3_2C4_GB_f)

GR_list_toPlot_H3K36me3 = c("GR_UP_5_Q_H3K36me3_2C4_GB_f", "GR_UP_1_Q_H3K36me3_2C4_GB_f", "GR_CTRL_60_65_Q_H3K36me3_2C4_GB_f", "GR_DN_5_Q_H3K36me3_2C4_GB_f")

#####################################################################################-
#         PLOT  ----
#####################################################################################-

#### Q_H3K36me3_2C4_GB_f
## H2AV REPLICAT 1 - AVERAGE PLOT

pdf(paste0(workdir,"FIGURES/AVG_PROF/H2AV/AVG_PROF_H2AV_WT_NELF_KD_R1.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H2AV_PW_1,H2AV_PN_1),tmp,GR,c(0,6),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()



## H2AV REPLICAT 2 - AVERAGE PLOT

pdf(paste0(workdir,"FIGURES/AVG_PROF/H2AV/AVG_PROF_H2AV_WT_NELF_KD_R2.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H2AV_PW_2,H2AV_PN_2),tmp,GR,c(0,6),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()



### HYPB KD

#### Q_H3K36me3_2C4_GB_f
## H2AV REPLICAT 1 - AVERAGE PLOT

pdf(paste0(workdir,"FIGURES/AVG_PROF/H2AV/AVG_PROF_H2AV_WT_HYPB_KD_R1.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H2AV_PW_1,H2AV_PHYPB_R1),tmp,GR,c(0,6),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()



## H2AV REPLICAT 2 - AVERAGE PLOT

pdf(paste0(workdir,"FIGURES/AVG_PROF/H2AV/AVG_PROF_H2AV_WT_HYPB_KD_R2.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H2AV_PW_2,H2AV_PHYPB_R2),tmp,GR,c(0,6),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()






#### Q_H3K36me3_2C4_GB_f
## H2AV REPLICAT 1 - AVERAGE PLOT

pdf(paste0(workdir,"FIGURES/AVG_PROF/H2AV/AVG_PROF_H2AV_WT_IHYP_B_L1_2019.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H2AV_PW_1,IHYP_B_L1),tmp,GR,c(0,6),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()

## H2AV REPLICAT 2 - AVERAGE PLOT

pdf(paste0(workdir,"FIGURES/AVG_PROF/H2AV/AVG_PROF_H2AV_WT_IHYP_B_L2_2019.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H2AV_PW_2,IHYP_B_L2),tmp,GR,c(0,6),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()



#######################"

pdf(paste0(workdir,"FIGURES/AVG_PROF/H2AV/AVG_PROF_H2AV_WT_IHYPH_B_L1_2019.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H2AV_PW_1,IHYPH_B_L1),tmp,GR,c(0,6),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()


pdf(paste0(workdir,"FIGURES/AVG_PROF/H2AV/AVG_PROF_H2AV_WT_IHYPH_B_L2_2019.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H2AV_PW_2,IHYPH_B_L2),tmp,GR,c(0,6),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()


##############################


pdf(paste0(workdir,"FIGURES/AVG_PROF/H2AV/AVG_PROF_H2AV_WT_PHYP_B_L1_2019.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H2AV_PW_1,PHYP_B_L1),tmp,GR,c(0,6),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()


pdf(paste0(workdir,"FIGURES/AVG_PROF/H2AV/AVG_PROF_H2AV_WT_PHYP_B_L2_2019.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H2AV_PW_2,PHYP_B_L2),tmp,GR,c(0,6),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()


##############################


pdf(paste0(workdir,"FIGURES/AVG_PROF/H2AV/AVG_PROF_H2AV_WT_PHYP_B_L1_inputnorm_2019.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H2AV_PW_1,PHYP_B_L1_inputnorm),tmp,GR,c(0,6),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()


pdf(paste0(workdir,"FIGURES/AVG_PROF/H2AV/AVG_PROF_H2AV_WT_PHYP_B_L2_inputnorm_2019.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H2AV_PW_2,PHYP_B_L2_inputnorm),tmp,GR,c(0,6),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()


##############################


pdf(paste0(workdir,"FIGURES/AVG_PROF/H2AV/AVG_PROF_H2AV_WT_H2AV_R3_HypB_KD.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H2AV_R3_Luc_KD,H2AV_R3_HypB_KD),tmp,GR,c(0,5),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()


pdf(paste0(workdir,"FIGURES/AVG_PROF/H2AV/AVG_PROF_H2AV_WT_H2AV_R3_HypB_HU_KD.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H2AV_R3_Luc_HU_KD,H2AV_R3_HypB_HU_KD),tmp,GR,c(0,5),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()



