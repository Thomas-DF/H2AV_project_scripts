library(GenomicRanges)
library(rtracklayer)
library(EnrichedHeatmap)
library(ggplot2)
library(dplyr)

workdir = "~/Bureau/tdefreitas_genobioinfo/PROJET_H2AV_2025/"

r6_ref_genes = readRDS(paste0(workdir,"DATA/r6.13/TxDb.GR.dm6.RDS"))


H2AV_SEA4_Luc_KD = paste0(workdir,"DATA/CHIPSEQ/H2AV-R2_SEA4_Luc-KD_RPGC.bw")
H2AV_SEA4_Nelf_KD = paste0(workdir,"DATA/CHIPSEQ/H2AV-R2_SEA4_Nelf-KD_RPGC.bw")

H2AV_WT3_Luc_KD = paste0(workdir,"DATA/CHIPSEQ/H2AV-R2_WT3_Luc-KD_RPGC.bw")
H2AV_WT3_Nelf_KD = paste0(workdir,"DATA/CHIPSEQ/H2AV-R2_WT3_Nelf-KD_RPGC.bw")


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

plot_avg_profile <- function(bw_list, gr_list, bin = 50, anchor = 10000, smooth = TRUE, spar = 0.2, colvec = c("black", "firebrick2"), ylim = c(0, 7)) {
  
  mat_list <- list()
  bw_names <- sapply(substitute(bw_list)[-1], deparse)
  
  for (i in seq_along(bw_list)) {
    bw <- bw_list[i]
    mat <- normalizeToMatrix(import(bw), gr_list, extend = c(anchor, anchor), w = bin, mean_mode = "w0", value_column = "score")
    mat[is.na(mat)] <- 0
    sd_threshold <- 3
    mat_clean <- mat
    mat_clean[abs(scale(mat)) > sd_threshold] <- NA
    avg_signal <- colMeans(mat_clean, na.rm = TRUE)
    mat_list[[bw_names[i]]] <- avg_signal
  }
  
  df <- do.call(cbind, mat_list)
  df <- as.data.frame(df)
  df$position <- seq(-anchor, anchor, length.out = nrow(df))
  
  df_long <- reshape2::melt(df, id.vars = "position", variable.name = "Condition", value.name = "Signal")
  
  p <- ggplot(df_long, aes(x = position, y = Signal, color = Condition, fill = Condition)) +
    scale_x_continuous(name = "Relative position [bp]", limits = c(-anchor, anchor)) +
    scale_y_continuous(name = "Signal", limits = ylim) +
    scale_color_manual(values = colvec) +
    scale_fill_manual(values = colvec) + 
    theme_minimal() +
    ggtitle(paste0("Average Profile - ", deparse(substitute(gr_list))))+
    theme(legend.position = c(0.95, 0.95), legend.justification = c(1, 1)) +  
    guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))
  
  if (smooth) {
    p <- p + geom_smooth(method = "loess", span = spar, alpha = 0.2, level = 0.99)
  } else {
    p <- p + geom_line(size = 1)
  }
  
  print(p)
}





pdf(paste0(workdir,"FIGURES/AVG_PROF/AVG_PROF_H2AV_SEA4_Luc_KD_Nelf_KD.pdf"))
plot_avg_profile(
  bw_list = c(H2AV_SEA4_Luc_KD, H2AV_SEA4_Nelf_KD), gr_list = GR_UP_5_ZSCORE_H3K36me3_2C4_vs_2N4, bin = 100, anchor = 10000, 
  smooth = TRUE, spar = 0.20, colvec = c("#285bad", "#eb3434"), ylim = c(0, 9))
plot_avg_profile(
  bw_list = c(H2AV_SEA4_Luc_KD, H2AV_SEA4_Nelf_KD), gr_list = GR_UP_1_ZSCORE_H3K36me3_2C4_vs_2N4, bin = 100, anchor = 10000, 
  smooth = TRUE, spar = 0.20, colvec = c("#285bad", "#eb3434"), ylim = c(0, 9))
plot_avg_profile(
  bw_list = c(H2AV_SEA4_Luc_KD, H2AV_SEA4_Nelf_KD), gr_list = GR_CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4, bin = 100, anchor = 10000, 
  smooth = TRUE, spar = 0.20, colvec = c("#285bad", "#eb3434"), ylim = c(0, 9)
)
dev.off()





pdf(paste0(workdir,"FIGURES/AVG_PROF/AVG_PROF_H2AV_WT3_Luc_KD_Nelf_KD.pdf"))
plot_avg_profile(
  bw_list = c(H2AV_WT3_Luc_KD, H2AV_WT3_Nelf_KD), gr_list = GR_UP_5_ZSCORE_H3K36me3_2C4_vs_2N4, bin = 100, anchor = 10000, 
  smooth = TRUE, spar = 0.20, colvec = c("#285bad", "#eb3434"), ylim = c(0, 9))
plot_avg_profile(
  bw_list = c(H2AV_WT3_Luc_KD, H2AV_WT3_Nelf_KD), gr_list = GR_UP_1_ZSCORE_H3K36me3_2C4_vs_2N4, bin = 100, anchor = 10000, 
  smooth = TRUE, spar = 0.20, colvec = c("#285bad", "#eb3434"), ylim = c(0, 9))
plot_avg_profile(
  bw_list = c(H2AV_WT3_Luc_KD, H2AV_WT3_Nelf_KD), gr_list = GR_CTRL_60_65_ZSCORE_H3K36me3_2C4_vs_2N4, bin = 100, anchor = 10000, 
  smooth = TRUE, spar = 0.20, colvec = c("#285bad", "#eb3434"), ylim = c(0, 9))
dev.off()

