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
seqlevels(r6_ref_genes) = gsub("chr", "", seqlevels(r6_ref_genes))

### BIGWIG FILES

pol2_nelf = import(paste0(workdir,"DATA/macs2/pol2_nelf_summits.bed"))
pol2_ctrl = import(paste0(workdir,"DATA/macs2/pol2_ctrl_summits.bed"))


pol2_nelf = resize(pol2_nelf, fix="center", width=100)
pol2_ctrl = resize(pol2_ctrl, fix="center", width=100)


H2AV_pol2_nelf_file = paste0(workdir,"DATA/CHIPSEQ/H2AV_2022/H2AV_PN_J_L1_RPGC.bw")
H2AV_pol2_wt_file = paste0(workdir,"DATA/CHIPSEQ/H2AV_2022/H2AV_PW_1_L1_RPGC.bw")

# chromosomes_to_keep <- c("2L", "2R", "3L", "3R", "4", "X", "Y")
# seqlevels(pol2_nelf) = chromosomes_to_keep
# seqlengths(pol2_nelf) <- c("2L"=23513712, "2R" = 25286936, "3L" = 28110227, "3R" =32079331,"4" = 1348131, "X"=23542271, "Y"=3667352)
# export(pol2_nelf, paste0(workdir,"DATA/macs2/pol2_nelf_summit.bw"), format = "BigWig")
# 
# seqlevels(pol2_ctrl, pruning.mode="coarse") = chromosomes_to_keep
# seqlengths(pol2_ctrl) <- c("2L"=23513712, "2R" = 25286936, "3L" = 28110227, "3R" =32079331,"4" = 1348131, "X"=23542271, "Y"=3667352)
# export(pol2_ctrl, paste0(workdir,"DATA/macs2/pol2_ctrl_summit.bw"), format = "BigWig")


pol2_nelf_file = paste0(workdir,"DATA/macs2/pol2_nelf_summits.bw")
pol2_ctrl_file = paste0(workdir,"DATA/macs2/pol2_ctrl_summits.bw")
#####################################################################################-
#         PEAK FILTER  ----
#####################################################################################-

overlaps = findOverlaps(pol2_ctrl, pol2_nelf)

common_peaks = pol2_nelf[unique(overlaps@to)]
luc_only = pol2_ctrl[-unique(overlaps@from)]
nelf_only = pol2_nelf[-unique(overlaps@to)]


#####################################################################################-
#         GROUP SELECTION  ----
#####################################################################################-

r6_ref_genes_GB = r6_ref_genes [width(r6_ref_genes)>500]
starrseq_ENHANCER = readRDS(paste0(workdir,"DATA/r6.13/starrseq_ENHANCER_zabidi.RDS"))
r6_ref_genes_PROM = promoters(r6_ref_genes, upstream=250, downstream=0 )
r6_ref_genes_TES = terminators(r6_ref_genes, upstream=0, downstream=250 )
r6_ref_genes_ALL = c(r6_ref_genes_GB,starrseq_ENHANCER,r6_ref_genes_PROM,r6_ref_genes_TES)

luc_only_GB = luc_only[queryHits(findOverlaps(luc_only, r6_ref_genes_GB))]
luc_only_Other = luc_only[-queryHits(findOverlaps(luc_only, r6_ref_genes_ALL))]

nelf_only_GB = nelf_only[queryHits(findOverlaps(nelf_only, r6_ref_genes_GB))]
nelf_only_Other = nelf_only[-queryHits(findOverlaps(nelf_only, r6_ref_genes_ALL))]


#####################################################################################-

### GENE GROUPES 
common_peaks = resize(common_peaks, fix = "center", width = 1)

nelf_only = resize(nelf_only, fix = "center", width = 1)
nelf_only_GB = resize(nelf_only_GB, fix = "center", width = 1)
nelf_only_Other = resize(nelf_only_Other, fix = "center", width = 1)

luc_only = resize(luc_only, fix = "center", width = 1)
luc_only_GB = resize(luc_only_GB, fix = "center", width = 1)
luc_only_Other = resize(luc_only_Other, fix = "center", width = 1)

GR_list_toPlot_peaks = c("common_peaks","nelf_only","nelf_only_GB", "nelf_only_Other","luc_only", "luc_only_GB", "luc_only_Other")

GR_list_toPlot_nelf = c("nelf_only","nelf_only_GB", "nelf_only_Other")

GR_list_toPlot_luc = c("luc_only", "luc_only_GB", "luc_only_Other")


#####################################################################################-
#         PLOT  ----
#####################################################################################-


pdf(paste0(workdir,"FIGURES/AVG_PROF/POL2_PEAKS/AVG_H2AV_POL2_PEAKS.pdf"))
for(GR in GR_list_toPlot_peaks){
  seqPlotSDoutliers_scaleFact(c(H2AV_pol2_nelf_file,H2AV_pol2_wt_file),tmp,GR,c(0,3),c(2000,2000),type="mf",bin=10,
                              smooth=TRUE,spar=0.3, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#594d32", "#e09d00")) 
}
dev.off()



pdf(paste0(workdir,"FIGURES/AVG_PROF/POL2_PEAKS/AVG_H2AV_POL2_PEAKS_NELF.pdf"))
for(GR in GR_list_toPlot_nelf){
  seqPlotSDoutliers_scaleFact(c(H2AV_pol2_nelf_file,H2AV_pol2_wt_file),tmp,GR,c(0,3),c(2000,2000),type="mf",bin=10,
                              smooth=TRUE,spar=0.3, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#594d32", "#e09d00")) 
}
dev.off()


pdf(paste0(workdir,"FIGURES/AVG_PROF/POL2_PEAKS/AVG_H2AV_POL2_PEAKS_LUC.pdf"))
for(GR in GR_list_toPlot_luc){
  seqPlotSDoutliers_scaleFact(c(H2AV_pol2_nelf_file,H2AV_pol2_wt_file),tmp,GR,c(0,3),c(2000,2000),type="mf",bin=10,
                              smooth=TRUE,spar=0.3, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#594d32", "#e09d00")) 
}
dev.off()







pdf(paste0(workdir,"FIGURES/AVG_PROF/POL2_PEAKS/AVG_POL2_PEAKS.pdf"))
for(GR in GR_list_toPlot_peaks){
  seqPlotSDoutliers_scaleFact(c(pol2_ctrl_file,pol2_nelf_file),tmp,GR,c(0,3),c(2000,2000),type="mf",bin=10,
                              smooth=TRUE,spar=0.3, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#594d32", "#e09d00")) 
}
dev.off()

