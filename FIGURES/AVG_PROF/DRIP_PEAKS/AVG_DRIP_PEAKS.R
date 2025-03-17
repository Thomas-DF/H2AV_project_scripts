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


### BED FILES

pol2_nelf = import(paste0(workdir,"DATA/macs2/pol2_nelf_summits.bed"))
pol2_ctrl = import(paste0(workdir,"DATA/macs2/pol2_ctrl_summits.bed"))

### RESIZE PEAKS

pol2_nelf = resize(pol2_nelf, fix="start", 50)
pol2_nelf= resize(pol2_nelf, fix="end", 100)

pol2_ctrl = resize(pol2_ctrl, fix="start", 50)
pol2_ctrl= resize(pol2_ctrl, fix="end", 100)


### PEAKS FILTER

overlaps = findOverlaps(pol2_ctrl, pol2_nelf)

common_peaks = pol2_ctrl[queryHits(overlaps)]
luc_only = pol2_ctrl[!pol2_ctrl %in% common_peaks]
nelf_only = pol2_nelf[!pol2_nelf %in% common_peaks]

seqlevels(common_peaks) = paste0("chr", seqlevels(common_peaks))  
seqlevels(luc_only) = paste0("chr", seqlevels(luc_only))
seqlevels(nelf_only) = paste0("chr", seqlevels(nelf_only))


### SAVING 

chromosomes_to_keep <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")
common_peaks = common_peaks[seqnames(common_peaks) %in% chromosomes_to_keep]
luc_only = luc_only[seqnames(luc_only) %in% chromosomes_to_keep]
nelf_only = nelf_only[seqnames(nelf_only) %in% chromosomes_to_keep]

seqlevels(common_peaks) <- chromosomes_to_keep
seqlevels(luc_only) <- chromosomes_to_keep
seqlevels(nelf_only) <- chromosomes_to_keep

seqlengths(common_peaks) <- c(chr2L=23513712, chr2R = 25286936, chr3L = 28110227, chr3R =32079331, chr4 = 1348131, chrX=23542271, chrY=3667352)
export(common_peaks, paste0(workdir,"DATA/macs2/common_peaks_summit.bw"), format = "BigWig")

seqlengths(luc_only) <- c(chr2L=23513712, chr2R = 25286936, chr3L = 28110227, chr3R =32079331, chr4 = 1348131, chrX=23542271, chrY=3667352)
export(luc_only, paste0(workdir,"DATA/macs2/common_peaks_summit.bw"), format = "BigWig")

seqlengths(nelf_only) <- c(chr2L=23513712, chr2R = 25286936, chr3L = 28110227, chr3R =32079331,chr4 = 1348131, chrX=23542271, chrY=3667352)
export(nelf_only, paste0(workdir,"DATA/macs2/common_peaks_summit.bw"), format = "BigWig")


### LOADING DATA

common_peaks = paste0(workdir,"DATA/macs2/common_peaks_summit.bw")
luc_only = paste0(workdir,"DATA/macs2/luc_only_summit.bw")
nelf_only = paste0(workdir,"DATA/macs2/nelf_only_summit.bw")


#####################################################################################-

### GENE GROUPES 

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

DN_5_ZSCORE_H3K36me3_2C4_vs_2N4 = getNameList(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, topdown = "down", prct = 5)



## GET GRANGES COORDS OF FEATURES TO PLOT AROUND
get_GR_feat = function(refGN, GNlist){
  myovlp=refGN[refGN$name %in% GNlist]
  return(myovlp)
}

GR_UP_5_ZSCORE_H3K36me3_2C4_vs_2N4=get_GR_feat(r6_ref_genes,UP_5_ZSCORE_H3K36me3_2C4_vs_2N4)
GR_UP_1_ZSCORE_H3K36me3_2C4_vs_2N4=get_GR_feat(r6_ref_genes,UP_1_ZSCORE_H3K36me3_2C4_vs_2N4)
GR_CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4=get_GR_feat(r6_ref_genes,CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4)
GR_DN_5_ZSCORE_H3K36me3_2C4_vs_2N4=get_GR_feat(r6_ref_genes,DN_5_ZSCORE_H3K36me3_2C4_vs_2N4)

GR_list_toPlot = c("GR_UP_5_ZSCORE_H3K36me3_2C4_vs_2N4", "GR_UP_1_ZSCORE_H3K36me3_2C4_vs_2N4", "GR_CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4", "GR_DN_5_ZSCORE_H3K36me3_2C4_vs_2N4")


#####################################################################################-
#         PLOT  ----
#####################################################################################-


pdf(paste0(workdir,"FIGURES/AVG_PROF/DRIP_PEAKS/AVG_PROF_PEAKS_common_nelf_luc.pdf"))
for(GR in GR_list_toPlot){
  seqPlotSDoutliers_scaleFact(c(common_peaks, luc_only, nelf_only),tmp,GR,c(0,100),c(5000,5000),type="af",bin=10,
                              smooth=FALSE,spar=0.3, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434","#228B22")) 
}
dev.off()


