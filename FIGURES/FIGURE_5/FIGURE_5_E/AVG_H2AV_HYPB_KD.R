#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-


library(dplyr) # used for recover names(coverages)
library(Rsamtools)
library(GenomicRanges)
library(GenomeInfoDb)
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library(seqplots.tweaked)
library("BSgenome.Dmelanogaster.UCSC.dm6")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(rtracklayer)
library(devtools)
'%ni%' = Negate('%in%')
library("BSgenome")


#####################################################################################-
#         FUNCTIONS  ----
#####################################################################################-
# Load FUNCTION

workdir="/work/user/cperrois/"
fonctiondir="/work/user/cperrois/functionR/"
source(paste0(fonctiondir, "AVG_PROFILE.R"))



#####################################################################################-
#         DATA  ----
#####################################################################################-

r6_ref_genes=readRDS('/work/user/cperrois/PROJET_H2AV_2023/DATA/r6.13/TxDb.GR.dm6.RDS')

tmp <- paste0(workdir,"PROJET_H2AV_2023/FIGURE/AVG_PROF/TMPgetPlotSetArray")
bwdir = paste0(workdir, "PROJET_H2AV_2023/DATA/CHIPSEQ/BIGWIG/")

### LOAD BIGWIG

BW_PHYPBH_R3_CPM=paste0(bwdir,"HypB_KD_HU_None_norm.bw")
BW_PHYPB_R3_CPM=paste0(bwdir,"HypB_KD_None_norm.bw")
BW_PWH_R3_CPM=paste0(bwdir,"Luc_KD_HU_None_norm.bw")
BW_PW_R3_CPM=paste0(bwdir,"LucKD_None_norm.bw")


## GENENS GROUPES 

LIST_QUANTIF_K36=readRDS("/work/user/cperrois/PROJET_H2AV_2023/DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS")
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
UP10_Q_H3K36me3_2C4=getNameList(Q_H3K36me3_2C4_GB_f, topdown = "top", prct = 10)
DN10_Q_H3K36me3_2C4=getNameList(Q_H3K36me3_2C4_GB_f, topdown = "down", prct = 10)

## GET GRANGES COORDS OF FEATURES TO PLOT AROUND
get_GR_feat = function(refGN, GNlist){
myovlp=refGN[refGN$name %in% GNlist]
  return(myovlp)
}

GR_UP10_Q_H3K36me3_2C4=get_GR_feat(r6_ref_genes,UP10_Q_H3K36me3_2C4)
GR_DN10_Q_H3K36me3_2C4=get_GR_feat(r6_ref_genes,DN10_Q_H3K36me3_2C4)
 
seqlevelsStyle(GR_UP10_Q_H3K36me3_2C4) = "ensembl"
seqlevelsStyle(GR_DN10_Q_H3K36me3_2C4) = "ensembl"


GR_list_toPlot=c(
"GR_UP10_Q_H3K36me3_2C4",
"GR_DN10_Q_H3K36me3_2C4"
)


#####################################################################################-
#         PLOT  ----
#####################################################################################-

for(GR in GR_list_toPlot){
  pdf(paste0(workdir,"PROJET_H2AV_2023/FIGURE/AVG_PROF/HYPB/AVG_PLOT_H2AV_HYPB_KD_R3_CPM_" ,GR,".pdf"))
  par(lwd=2)
  seqPlotSDoutliers_scaleFact(c(BW_PW_R3_CPM, BW_PHYPB_R3_CPM),tmp,GR,c(-1,10),c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
  seqPlotSDoutliers_scaleFact(c(BW_PWH_R3_CPM, BW_PHYPBH_R3_CPM),tmp,GR,c(-1,20),c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
  dev.off()
}
