library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)

workdir = "~/Bureau/tdefreitas_genobioinfo/PROJET_H2AV_2025/"

r6_ref_genes = readRDS(paste0(workdir,"DATA/r6.13/TxDb.GR.dm6.RDS"))
seqlevels(r6_ref_genes) = gsub("chr", "", seqlevels(r6_ref_genes))


H2AV_SEA4_Luc_KD = import(paste0(workdir,"DATA/CHIPSEQ/H2AV-R2_SEA4_Luc-KD_RPGC.bw"), which = r6_ref_genes, as = "NumericList")
H2AV_SEA4_Nelf_KD = import(paste0(workdir,"DATA/CHIPSEQ/H2AV-R2_SEA4_Nelf-KD_RPGC.bw"), which = r6_ref_genes, as = "NumericList")

H2AV_WT3_Luc_KD = import(paste0(workdir,"DATA/CHIPSEQ/H2AV-R2_WT3_Luc-KD_RPGC.bw"), which = r6_ref_genes, as = "NumericList")
H2AV_WT3_Nelf_KD = import(paste0(workdir,"DATA/CHIPSEQ/H2AV-R2_WT3_Nelf-KD_RPGC.bw"), which = r6_ref_genes, as = "NumericList")



ggplot(H2AV_SEA4_Luc_KD) + geom_line()



LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f=LIST_QUANTIF_K36$ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f

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



## GET GRANGES COORDS OF FEATURES TO PLOT AROUND
get_GR_feat = function(refGN, GNlist){
  myovlp=refGN[refGN$name %in% GNlist]
  return(myovlp)
}

GR_UP_5_ZSCORE_H3K36me3_2C4_vs_2N4=get_GR_feat(r6_ref_genes,UP_5_ZSCORE_H3K36me3_2C4_vs_2N4)
GR_UP_1_ZSCORE_H3K36me3_2C4_vs_2N4=get_GR_feat(r6_ref_genes,UP_1_ZSCORE_H3K36me3_2C4_vs_2N4)
GR_CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4=get_GR_feat(r6_ref_genes,CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4)

seqlevels(GR_UP_5_ZSCORE_H3K36me3_2C4_vs_2N4) = gsub("chr", "", seqlevels(GR_UP_5_ZSCORE_H3K36me3_2C4_vs_2N4))
seqlevels(GR_UP_1_ZSCORE_H3K36me3_2C4_vs_2N4) = gsub("chr", "", seqlevels(GR_UP_1_ZSCORE_H3K36me3_2C4_vs_2N4))
seqlevels(GR_CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4) = gsub("chr", "", seqlevels(GR_CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4))


