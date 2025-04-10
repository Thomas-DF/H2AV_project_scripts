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

pol2_ctrl_N_RPGC = import(paste0(workdir, "DATA/CHIPSEQ/pol2_ctrl_N_filt_sort_RPGC.bw"))
pol2_nelf_N_RPGC = import(paste0(workdir, "DATA/CHIPSEQ/pol2_nelf_N_filt_sort_RPGC.bw"))


seqlevels(pol2_ctrl_N_RPGC) = paste0("chr", seqlevels(pol2_ctrl_N_RPGC))
seqlevels(pol2_nelf_N_RPGC) = paste0("chr", seqlevels(pol2_nelf_N_RPGC))

#####################################################################################-

### GENE GROUPES 

PAUSE_INDICE_VEC = readRDS(paste0(workdir, "DATA/LIST_FEATURES/PAUSE_INDICE_VEC.RDS"))
DELTA_PAUSE_IND_pol2_nelf_N = PAUSE_INDICE_VEC$DELTA_PAUSE_IND_pol2_nelf_N

PAUSE_IND_pol2_ctrl_N = PAUSE_INDICE_VEC$PAUSE_IND_pol2_ctrl_N


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


len = length(DELTA_PAUSE_IND_pol2_nelf_N)
CTRL_60_65_DELTA_PAUSE_IND_pol2_nelf_N = names(DELTA_PAUSE_IND_pol2_nelf_N[(len*0.6) : (len*0.65)])
UP_5_DELTA_PAUSE_IND_pol2_nelf_N = getNameList(DELTA_PAUSE_IND_pol2_nelf_N, topdown = "top", prct = 5)
UP_1_DELTA_PAUSE_IND_pol2_nelf_N = getNameList(DELTA_PAUSE_IND_pol2_nelf_N, topdown = "top", prct = 1)
DN_5_DELTA_PAUSE_IND_pol2_nelf_N = getNameList(DELTA_PAUSE_IND_pol2_nelf_N, topdown = "down", prct = 5)
DN_1_DELTA_PAUSE_IND_pol2_nelf_N = getNameList(DELTA_PAUSE_IND_pol2_nelf_N, topdown = "down", prct = 1)

len = length(PAUSE_IND_pol2_ctrl_N)
CTRL_60_65_PAUSE_IND_pol2_ctrl_N = names(PAUSE_IND_pol2_ctrl_N[(len*0.6) : (len*0.65)])
UP_5_PAUSE_IND_pol2_ctrl_N = getNameList(PAUSE_IND_pol2_ctrl_N, topdown = "top", prct = 5)
UP_1_PAUSE_IND_pol2_ctrl_N = getNameList(PAUSE_IND_pol2_ctrl_N, topdown = "top", prct = 1)
DN_5_PAUSE_IND_pol2_ctrl_N = getNameList(PAUSE_IND_pol2_ctrl_N, topdown = "down", prct = 5)




## GET GRANGES COORDS OF FEATURES TO PLOT AROUND
get_GR_feat = function(refGN, GNlist){
  myovlp=refGN[refGN$name %in% GNlist]
  return(myovlp)
}

GR_UP_5_DELTA_PAUSE_IND_pol2_nelf_N = get_GR_feat(r6_ref_genes,UP_5_DELTA_PAUSE_IND_pol2_nelf_N)
GR_UP_1_DELTA_PAUSE_IND_pol2_nelf_N = get_GR_feat(r6_ref_genes,UP_1_DELTA_PAUSE_IND_pol2_nelf_N)
GR_CTRL_60_65_DELTA_PAUSE_IND_pol2_nelf_N = get_GR_feat(r6_ref_genes,CTRL_60_65_DELTA_PAUSE_IND_pol2_nelf_N)
GR_DN_5_DELTA_PAUSE_IND_pol2_nelf_N = get_GR_feat(r6_ref_genes,DN_5_DELTA_PAUSE_IND_pol2_nelf_N)
GR_DN_1_DELTA_PAUSE_IND_pol2_nelf_N = get_GR_feat(r6_ref_genes,DN_1_DELTA_PAUSE_IND_pol2_nelf_N)


GR_UP_5_PAUSE_IND_pol2_ctrl_N = get_GR_feat(r6_ref_genes,UP_5_PAUSE_IND_pol2_ctrl_N)
GR_UP_1_PAUSE_IND_pol2_ctrl_N = get_GR_feat(r6_ref_genes,UP_1_PAUSE_IND_pol2_ctrl_N)
GR_CTRL_60_65_PAUSE_IND_pol2_ctrl_N = get_GR_feat(r6_ref_genes,CTRL_60_65_PAUSE_IND_pol2_ctrl_N)
GR_DN_5_PAUSE_IND_pol2_ctrl_N = get_GR_feat(r6_ref_genes,DN_5_PAUSE_IND_pol2_ctrl_N)



#####################################################################################-
#         GROUP SELECTION  ----
#####################################################################################-

r6_ref_genes = readRDS(paste0(workdir,"DATA/r6.13/TxDb.GR.dm6.RDS"))
r6_ref_genes_GB = r6_ref_genes [width(r6_ref_genes)>500]

starrseq_ENHANCER = readRDS(paste0(workdir,"DATA/r6.13/starrseq_ENHANCER_zabidi.RDS"))
seqlevels(starrseq_ENHANCER) = paste0("chr", seqlevels(starrseq_ENHANCER))

r6_ref_genes_PROM = promoters(r6_ref_genes, upstream=250, downstream=0 )

r6_ref_genes_TES = terminators(r6_ref_genes, upstream=0, downstream=250 )


#####################################################################################-
#         FUNCTIONS  ----
#####################################################################################-

count_peak_categories <- function(peaks, promoters, gene_bodies, terminators, enhancers) {
  
  peak_categories <- rep("Other", length(peaks))
  
  peak_categories[queryHits(findOverlaps(peaks, enhancers))] <- "Enhancer"
  peak_categories[queryHits(findOverlaps(peaks, promoters))] <- "Promoter"
  peak_categories[queryHits(findOverlaps(peaks, terminators))] <- "Terminator"
  peak_categories[queryHits(findOverlaps(peaks, gene_bodies))] <- "Gene Body"
  
  df_summary <- as.data.frame(table(peak_categories))
  colnames(df_summary) <- c("Category", "Count")
  df_summary$Percentage <- round(df_summary$Count / sum(df_summary$Count) * 100, 1)
  
  return(df_summary)
}


Plot_piechart_peaks = function(peaks, promoters, gene_bodies, terminators, enhancers, title = NULL) {
  
  df = count_peak_categories(peaks, promoters, gene_bodies, terminators, enhancers)
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

### NELF / LUC COMMON PEAKS

pdf(paste0(workdir, "FIGURES/PIE_CHART/PIE_CHART_POL2/"))
Plot_piechart_peaks(pol2_ctrl_N_RPGC, r6_ref_genes_PROM, r6_ref_genes_GB, r6_ref_genes_TES, starrseq_ENHANCER,title = "Peaks repartition pie chart -")
grid.newpage()
grid.table(count_peak_categories(common_peaks, r6_ref_genes_PROM, r6_ref_genes_GB, r6_ref_genes_TES, starrseq_ENHANCER))
dev.off()


### LUC ONLY PEAKS

pdf(paste0(workdir, "FIGURES/PIE_CHART/PIE_CHART_POL2/"))
Plot_piechart_peaks(luc_only, r6_ref_genes_PROM, r6_ref_genes_GB, r6_ref_genes_TES, starrseq_ENHANCER,title = "Peaks repartition pie chart - luc_only")
grid.newpage()
grid.table(count_peak_categories(luc_only, r6_ref_genes_PROM, r6_ref_genes_GB, r6_ref_genes_TES, starrseq_ENHANCER))
dev.off()
















#############################################################################
df = data.frame(Category = c("Luc only", "Nelf only", "Common peaks"), Count = c(length(luc_only),length(nelf_only),length(common_peaks)))
df$Percentage = round(df$Count / sum(df$Count) * 100, 1) 
df


df = df %>% arrange(desc(Category)) %>% mutate(ypos = cumsum(Count) - Count / 2)


pdf(paste0(workdir, "FIGURES/PIE_CHART/PIE_CHART_POL2/PIECHART_Nelf_Luc.pdf"))

ggplot(df, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
  coord_polar("y", start = 0) +
  theme_void() +
  ggtitle("Peaks repartition pie chart - Luc / Nelf") +
  scale_fill_manual(values = c("#FF6F61", "#92A8D1", "#88B04B", "#F7CAE9", "#6B5B95")) +
  
  geom_label_repel(aes(y = ypos, label = paste0(Percentage, "%")),
                   nudge_x = 0.,
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


grid.newpage()
grid.table(select(df,-ypos))
dev.off()
