#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

library(GenomicRanges)
library(ggplot2)
library(GenomicFeatures)
library(rtracklayer)
library(ggrepel)
library(dplyr)
'%ni%' = Negate('%in%')
library(gridExtra)
library(grid)

#####################################################################################-
#         DATA  ----
#####################################################################################-

workdir = "~/Bureau/tdefreitas_genobioinfo/PROJET_H2AV_2025/"


DRIP_Nelf_summits = import(paste0(workdir,"DATA/Peaks/DRIP_Nelf_summits.bed"))
DRIP_Nelf_RNAse_summits = import(paste0(workdir,"DATA/Peaks/DRIP_Nelf_RNAse_summits.bed"))

DRIP_Luc_summits = import(paste0(workdir,"DATA/Peaks/DRIP_Luc_summits.bed"))
DRIP_Luc_RNAse_summits = import(paste0(workdir,"DATA/Peaks/DRIP_Luc_RNAse_summits.bed"))


#####################################################################################-
#         PEAK FILTER  ----
#####################################################################################-

overlaps = findOverlaps(DRIP_Luc_summits, DRIP_Luc_RNAse_summits)
filtered_Luc = DRIP_Luc_summits[-queryHits(overlaps)]
filtered_Luc = resize(filtered_Luc, fix="start", 50)
filtered_Luc= resize(filtered_Luc, fix="end", 100)

overlaps = findOverlaps(DRIP_Nelf_summits, DRIP_Nelf_RNAse_summits)
filtered_Nelf = DRIP_Nelf_summits[-queryHits(overlaps)]
filtered_Nelf = resize(filtered_Nelf, fix="start", 50)
filtered_Nelf= resize(filtered_Nelf, fix="end", 100)


#####################################################################################-
#         GROUP SELECTION  ----
#####################################################################################-

GN_actifs = readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))
r6_ref_genes = readRDS(paste0(workdir,"DATA/r6.13/TxDb.GR.dm6.RDS"))
r6_ref_genes_GB = r6_ref_genes [width(r6_ref_genes)>500]

starrseq_ENHANCER = readRDS(paste0(workdir,"DATA/r6.13/starrseq_ENHANCER_zabidi.RDS"))
seqlevels(starrseq_ENHANCER) = paste0("chr", seqlevels(starrseq_ENHANCER))

r6_ref_genes_PROM = promoters(r6_ref_genes, upstream=250, downstream=0 )

r6_ref_genes_TES = terminators(r6_ref_genes, upstream=0, downstream=250 )

r6_ref_genes_GB_actif = unique(r6_ref_genes_GB [r6_ref_genes_GB$name %in% GN_actifs])
r6_ref_genes_PROM_actif = r6_ref_genes_PROM [r6_ref_genes_PROM$name %in% GN_actifs]
r6_ref_genes_TES_actif = r6_ref_genes_TES [r6_ref_genes_TES$name %in% GN_actifs]


#####################################################################################-
#         FUNCTIONS  ----
#####################################################################################-

count_peak_categories <- function(peaks, promoters, gene_bodies, terminators, enhancers) {
  
  peak_categories <- rep("Other", length(peaks))
  
  peak_categories[queryHits(findOverlaps(peaks, gene_bodies))] <- "Gene Body"
  peak_categories[queryHits(findOverlaps(peaks, enhancers))] <- "Enhancer"
  peak_categories[queryHits(findOverlaps(peaks, promoters))] <- "Promoter"
  peak_categories[queryHits(findOverlaps(peaks, terminators))] <- "Terminator"
  
  df_summary <- as.data.frame(table(peak_categories))
  colnames(df_summary) <- c("Category", "Count")
  df_summary$Percentage <- round(df_summary$Count / sum(df_summary$Count) * 100, 1)
  
  return(df_summary)
}


Plot_piechart_peaks = function(peaks, promoters, gene_bodies, terminators, enhancers, title = NULL) {
  
  df = count_peak_categories(peaks, promoters, gene_bodies, terminators, enhancers)
  print(df)
  df = df %>% arrange(desc(Category)) %>% mutate(ypos = cumsum(Count) - Count / 2)
  
  ggplot(df, aes(x = "", y = Count, fill = Category)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
    coord_polar("y", start = 0) +
    theme_void() +
    ggtitle(title) +
    scale_fill_manual(values = c("#FF6F61", "#92A8D1", "#88B04B", "#F7CAE9", "#6B5B95")) +
    
    geom_label_repel(aes(y = ypos, label = paste0(Percentage, "%")),
                     nudge_x = 0.7,
                     size = 5, fontface = "bold",
                     color = "black",
                     fill = "white",
                     segment.color = "black",
                     segment.size = 0.5) +  
    
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    )
}



#####################################################################################-
#         PLOT  ----
#####################################################################################-

### NELF

pdf(paste0(workdir, "FIGURES/PIE_CHART/PIE_CHART_DRIP/PIECHART_DRIP_Nelf.pdf"))
Plot_piechart_peaks(filtered_Nelf, r6_ref_genes_PROM, r6_ref_genes_GB, r6_ref_genes_TES, starrseq_ENHANCER, title = "Peaks repartition pie chart - DRIP_Nelf_summits")
grid.newpage()
grid.table(count_peak_categories(filtered_Nelf, r6_ref_genes_PROM, r6_ref_genes_GB, r6_ref_genes_TES, starrseq_ENHANCER))
dev.off()


### LUC

pdf(paste0(workdir, "FIGURES/PIE_CHART/PIE_CHART_DRIP/PIECHART_DRIP_Luc.pdf"))
Plot_piechart_peaks(filtered_Luc, r6_ref_genes_PROM, r6_ref_genes_GB, r6_ref_genes_TES, starrseq_ENHANCER, title = "Peaks repartition pie chart - DRIP_Luc_summits")
grid.newpage()
grid.table(count_peak_categories(filtered_Luc, r6_ref_genes_PROM, r6_ref_genes_GB, r6_ref_genes_TES, starrseq_ENHANCER))
dev.off()






#####################################################################################-
#
#         COMPLEMENTARY DISTRIBUTION ANALYSE  ----
#
#####################################################################################-

######### GENE GROUPES ######### 

LIST_QUANTIF=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF.RDS"))

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
notop5 = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0) : (len*0.95)])
noctrl = c(names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0) : (len*0.6)]), names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0.65) : (len)]))

ZSCORE_00_10 = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0.0) : (len*0.1)])
ZSCORE_10_20 = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0.1) : (len*0.2)])
ZSCORE_20_30 = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0.2) : (len*0.3)])
ZSCORE_30_40 = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0.3) : (len*0.4)])
ZSCORE_40_50 = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0.4) : (len*0.5)])
ZSCORE_50_60 = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0.5) : (len*0.6)])
ZSCORE_60_70 = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0.6) : (len*0.7)])
ZSCORE_70_80 = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0.7) : (len*0.8)])
ZSCORE_80_90 = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0.8) : (len*0.9)])
ZSCORE_90_100 = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[(len*0.9) : (len*1)])

ZSCORE_ALL = names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f)


ol = findOverlaps(r6_ref_genes_GB_actif, filtered_Nelf)
peaks_GB = unique(r6_ref_genes_GB_actif[queryHits(ol)])
unique(peaks_GB[(peaks_GB$name %in% UP_5_ZSCORE_H3K36me3_2C4_vs_2N4)])

unique(peaks_GB[(peaks_GB$name %in% CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4)])

ol = findOverlaps(r6_ref_genes_GB_actif, filtered_Nelf)
peaks_GB = r6_ref_genes_GB_actif[queryHits(ol)]
unique(peaks_GB[(peaks_GB$name %in% notop5)])


# m = matrix(c(69,280,4845,1773), nrow = 2)
# fisher.test(m,alternative = "less")


######### FUNCTIONS ######### 

peaks_repartition <- function(peaks, group, selection){
  ol = findOverlaps(group, peaks)
  peaks_group = unique(group[queryHits(ol)])
  selected_peaks = unique(peaks_group[peaks_group$name %in% selection])
  peaks_percentage = length(selected_peaks)/length(peaks_group)*100
  peaks_percentage = round(peaks_percentage, 1)
  return(length(selected_peaks))
}


peaks_repartition_table <- function(peaks, promoters, gene_bodies, terminators, enhancers, ZSCORE){
  result = list()
  result["PROMOTER"] = peaks_repartition(peaks, promoters,ZSCORE)
  result["GENE BODY"] = peaks_repartition(peaks, gene_bodies,ZSCORE)
  result["TERMINATOR"] = peaks_repartition(peaks, terminators,ZSCORE)
  result["ENHANCER"] = peaks_repartition(peaks, enhancers,ZSCORE)
  df = as.data.frame(result)
  rownames(df) = deparse(substitute(ZSCORE))
  return(df)
}

p1 = peaks_repartition_table(filtered_Nelf, r6_ref_genes_PROM_actif, r6_ref_genes_GB_actif, r6_ref_genes_TES_actif, starrseq_ENHANCER,UP_5_ZSCORE_H3K36me3_2C4_vs_2N4)
p2 = peaks_repartition_table(filtered_Nelf, r6_ref_genes_PROM_actif, r6_ref_genes_GB_actif, r6_ref_genes_TES_actif, starrseq_ENHANCER,notop5)
p3 = peaks_repartition_table(filtered_Nelf, r6_ref_genes_PROM_actif, r6_ref_genes_GB_actif, r6_ref_genes_TES_actif, starrseq_ENHANCER,CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4)
p4 = peaks_repartition_table(filtered_Nelf, r6_ref_genes_PROM_actif, r6_ref_genes_GB_actif, r6_ref_genes_TES_actif, starrseq_ENHANCER,noctrl)
df = rbind(p1,p2,p3,p4)
df

p0 = peaks_repartition_table(filtered_Nelf, r6_ref_genes_PROM_actif, r6_ref_genes_GB_actif, r6_ref_genes_TES_actif, starrseq_ENHANCER,ZSCORE_00_10)
p1 = peaks_repartition_table(filtered_Nelf, r6_ref_genes_PROM_actif, r6_ref_genes_GB_actif, r6_ref_genes_TES_actif, starrseq_ENHANCER,ZSCORE_10_20)
p2 = peaks_repartition_table(filtered_Nelf, r6_ref_genes_PROM_actif, r6_ref_genes_GB_actif, r6_ref_genes_TES_actif, starrseq_ENHANCER,ZSCORE_20_30)
p3 = peaks_repartition_table(filtered_Nelf, r6_ref_genes_PROM_actif, r6_ref_genes_GB_actif, r6_ref_genes_TES_actif, starrseq_ENHANCER,ZSCORE_30_40)
p4 = peaks_repartition_table(filtered_Nelf, r6_ref_genes_PROM_actif, r6_ref_genes_GB_actif, r6_ref_genes_TES_actif, starrseq_ENHANCER,ZSCORE_40_50)
p5 = peaks_repartition_table(filtered_Nelf, r6_ref_genes_PROM_actif, r6_ref_genes_GB_actif, r6_ref_genes_TES_actif, starrseq_ENHANCER,ZSCORE_50_60)
p6 = peaks_repartition_table(filtered_Nelf, r6_ref_genes_PROM_actif, r6_ref_genes_GB_actif, r6_ref_genes_TES_actif, starrseq_ENHANCER,ZSCORE_60_70)
p7 = peaks_repartition_table(filtered_Nelf, r6_ref_genes_PROM_actif, r6_ref_genes_GB_actif, r6_ref_genes_TES_actif, starrseq_ENHANCER,ZSCORE_70_80)
p8 = peaks_repartition_table(filtered_Nelf, r6_ref_genes_PROM_actif, r6_ref_genes_GB_actif, r6_ref_genes_TES_actif, starrseq_ENHANCER,ZSCORE_80_90)
p9 = peaks_repartition_table(filtered_Nelf, r6_ref_genes_PROM_actif, r6_ref_genes_GB_actif, r6_ref_genes_TES_actif, starrseq_ENHANCER,ZSCORE_90_100)
df_total = rbind(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9)
df_total


######### SAVE TABLE #########

pdf(paste0(workdir, "FIGURES/PIE_CHART/PIE_CHART_DRIP/PIECHART_Nelf_Zscore_activ_gene_count_repartition.pdf"),width = 10, height = 10)
grid.table(df)
grid.newpage()
grid.table(df_total)
dev.off()


