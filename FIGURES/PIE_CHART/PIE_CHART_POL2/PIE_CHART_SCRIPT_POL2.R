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

pol2_nelf = import(paste0(workdir,"DATA/macs2/pol2_nelf_summits.bed"))
pol2_ctrl = import(paste0(workdir,"DATA/macs2/pol2_ctrl_summits.bed"))


pol2_nelf = resize(pol2_nelf, fix="start", 10)
pol2_nelf= resize(pol2_nelf, fix="end", 20)

pol2_ctrl = resize(pol2_ctrl, fix="start", 10)
pol2_ctrl= resize(pol2_ctrl, fix="end", 20)


#####################################################################################-
#         PEAK FILTER  ----
#####################################################################################-

overlaps = findOverlaps(pol2_ctrl, pol2_nelf)

common_peaks = pol2_ctrl[queryHits(overlaps)]
luc_only = pol2_ctrl[!pol2_ctrl %in% common_peaks]
nelf_only = pol2_nelf[!pol2_nelf %in% common_peaks]

seqlevels(common_peaks) = paste0("chr", seqlevels(common_peaks))  
seqlevels(luc_only) = paste0("chr", seqlevels(luc_only))
seqlevels(nelf_only) = paste0("chr", seqlevels(nelf_only))


#####################################################################################-
#         GROUP SELECTION  ----
#####################################################################################-

r6_ref_genes = readRDS(paste0(workdir,"DATA/r6.13/TxDb.GR.dm6.RDS"))
r6_ref_genes_GB = r6_ref_genes [width(r6_ref_genes)>500]

starrseq_ENHANCER = readRDS(paste0(workdir,"DATA/r6.13/starrseq_ENHANCER_zabidi.RDS"))
seqlevels(starrseq_ENHANCER) = paste0("chr", seqlevels(starrseq_ENHANCER))

r6_ref_genes_PROM = promoters(r6_ref_genes, upstream=250, downstream=0 )

r6_ref_genes_TES = promoters(r6_ref_genes, upstream=0, downstream=250 )


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

pdf(paste0(workdir, "FIGURES/PIE_CHART/PIE_CHART_POL2/PIECHART_common_peaks.pdf"))
Plot_piechart_peaks(common_peaks, r6_ref_genes_PROM, r6_ref_genes_GB, r6_ref_genes_TES, starrseq_ENHANCER,title = "Peaks repartition pie chart - common_peaks")
grid.newpage()
grid.table(count_peak_categories(common_peaks, r6_ref_genes_PROM, r6_ref_genes_GB, r6_ref_genes_TES, starrseq_ENHANCER))
dev.off()


### LUC ONLY PEAKS

pdf(paste0(workdir, "FIGURES/PIE_CHART/PIE_CHART_POL2/PIECHART_luc_only_peaks.pdf"))
Plot_piechart_peaks(luc_only, r6_ref_genes_PROM, r6_ref_genes_GB, r6_ref_genes_TES, starrseq_ENHANCER,title = "Peaks repartition pie chart - luc_only")
grid.newpage()
grid.table(count_peak_categories(luc_only, r6_ref_genes_PROM, r6_ref_genes_GB, r6_ref_genes_TES, starrseq_ENHANCER))
dev.off()


### NELF ONLY PEAKS

pdf(paste0(workdir, "FIGURES/PIE_CHART/PIE_CHART_POL2/PIECHART_nelf_only_peaks.pdf"))
Plot_piechart_peaks(nelf_only, r6_ref_genes_PROM, r6_ref_genes_GB, r6_ref_genes_TES, starrseq_ENHANCER,title = "Peaks repartition pie chart - nelf_only")
grid.newpage()
grid.table(count_peak_categories(nelf_only, r6_ref_genes_PROM, r6_ref_genes_GB, r6_ref_genes_TES, starrseq_ENHANCER))
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
