
# Heterochromatin spread / Single Cell RNA seq project
# Cuvier''s Team
# Depierre*
# 2018

# This script aims to plot boxplot for comparison of HTRCHR reads quanitf


######################################################
# Load LIBRARY
library(ggplot2)
library(ggpubr)
library(gplots)
'%ni%' = Negate('%in%')
require(GenomicRanges)
require(BiocGenerics)
require(parallel)


boxplot_wilcox_dom = function(quantifWT, quantifKD, filterGN, outdir, readQuantif = "readQuantif", Cond = "KD"){
  filterGN_IN = filterGN$IN
  filterGN_OUT = filterGN$OUT
  filterGN_BORDER = c(filterGN$UPSTREAM, filterGN$DOWNSTREAM)
  filterGN_DOWNSTREAM = filterGN$DOWNSTREAM
  filterGN_UPSTREAM = filterGN$UPSTREAM
  ## GGplot method -> DATA FRAME
  df_toPlot = as.data.frame(cbind(as.numeric(unname(quantifWT)), rep("all_Genes", length(quantifWT)), rep("WT", length(quantifWT))))
  df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD)), rep("all_Genes", length(quantifKD)), rep("KD", length(quantifKD)))))
  df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifWT[names(quantifWT) %in% filterGN_IN])), rep("IN", length(quantifWT[names(quantifWT) %in% filterGN_IN])), rep("WT", length(quantifWT[names(quantifWT) %in% filterGN_IN])))))
  df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD[names(quantifKD) %in% filterGN_IN])), rep("IN", length(quantifKD[names(quantifKD) %in% filterGN_IN])), rep("KD", length(quantifKD[names(quantifKD) %in% filterGN_IN])))))
  df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifWT[names(quantifWT) %in% filterGN_OUT])), rep("OUT", length(quantifWT[names(quantifWT) %in% filterGN_OUT])), rep("WT", length(quantifWT[names(quantifWT) %in% filterGN_OUT])))))
  df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD[names(quantifKD) %in% filterGN_OUT])), rep("OUT", length(quantifKD[names(quantifKD) %in% filterGN_OUT])), rep("KD", length(quantifKD[names(quantifKD) %in% filterGN_OUT])))))
  df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifWT[names(quantifWT) %in% filterGN_BORDER])), rep("BORDER", length(quantifWT[names(quantifWT) %in% filterGN_BORDER])), rep("WT", length(quantifWT[names(quantifWT) %in% filterGN_BORDER])))))
  df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD[names(quantifKD) %in% filterGN_BORDER])), rep("BORDER", length(quantifKD[names(quantifKD) %in% filterGN_BORDER])), rep("KD", length(quantifKD[names(quantifKD) %in% filterGN_BORDER])))))
  df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifWT[names(quantifWT) %in% filterGN_DOWNSTREAM])), rep("DOWNSTREAM", length(quantifWT[names(quantifWT) %in% filterGN_DOWNSTREAM])), rep("WT", length(quantifWT[names(quantifWT) %in% filterGN_DOWNSTREAM])))))
  df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD[names(quantifKD) %in% filterGN_DOWNSTREAM])), rep("DOWNSTREAM", length(quantifKD[names(quantifKD) %in% filterGN_DOWNSTREAM])), rep("KD", length(quantifKD[names(quantifKD) %in% filterGN_DOWNSTREAM])))))
  df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifWT[names(quantifWT) %in% filterGN_UPSTREAM])), rep("UPSTREAM", length(quantifWT[names(quantifWT) %in% filterGN_UPSTREAM])), rep("WT", length(quantifWT[names(quantifWT) %in% filterGN_UPSTREAM])))))
  df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD[names(quantifKD) %in% filterGN_UPSTREAM])), rep("UPSTREAM", length(quantifKD[names(quantifKD) %in% filterGN_UPSTREAM])), rep("KD", length(quantifKD[names(quantifKD) %in% filterGN_UPSTREAM])))))
	colnames(df_toPlot) = c("readsCount", "select", "condition")
  df_toPlot$readsCount = log(as.numeric(as.character(df_toPlot$readsCount))+1)
  ## BOXPLOT
  p = ggplot(df_toPlot, aes(x=select, y=readsCount, fill=condition))+
      geom_boxplot() + stat_compare_means(method = "wilcox.test")
      # + scale_fill_manual(breaks=c("WT", "KD"),values=c("grey", "red"))
  pdf(paste0(outdir, "BOXPLOT_",readQuantif ,"_reads_WTvs",Cond,".pdf"), width = 16, height = 10 )
  print(p)
  dev.off()

}






boxplot_wilcox_splitted = function(quantifWT, quantifKD, sort_value, outdir, readQuantif = "readQuantif", Cond = "KD", Nsplit = 5, info = NULL){
  	sort_value_splitted = split(rownames(sort_value), ceiling(seq_along(rownames(sort_value))/ceiling(length(rownames(sort_value))/Nsplit)))
  ## GGplot method -> DATA FRAME
  df_toPlot = as.data.frame(cbind(as.numeric(unname(quantifWT)), rep("all_Genes", length(quantifWT)), rep("WT", length(quantifWT))))
  df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD)), rep("all_Genes", length(quantifKD)), rep("KD", length(quantifKD)))))
  	for(i in 1:length(sort_value_splitted)){
		df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifWT[rownames(quantifWT) %in% sort_value_splitted[[i]],1])), rep(paste0("Q",i), length(quantifWT[rownames(quantifWT) %in% sort_value_splitted[[i]],1])), rep("WT", length(quantifWT[rownames(quantifWT) %in% sort_value_splitted[[i]],1])))))
		df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD[rownames(quantifKD) %in% sort_value_splitted[[i]],1])), rep(paste0("Q",i), length(quantifKD[rownames(quantifKD) %in% sort_value_splitted[[i]],1])), rep("KD", length(quantifKD[rownames(quantifKD) %in% sort_value_splitted[[i]],1])))))
	}

	colnames(df_toPlot) = c("readsCount", "select", "condition")
  df_toPlot$readsCount = log(as.numeric(as.character(df_toPlot$readsCount))+1)
  # df_toPlot$readsCount = as.numeric(as.character(df_toPlot$readsCount))
  ## BOXPLOT
  p = ggplot(df_toPlot, aes(x=select, y=readsCount, fill=condition))+
      geom_boxplot() + stat_compare_means(method = "wilcox.test", label = "p.signif") +
      geom_boxplot(outlier.shape = NA)
			# coord_cartesian(ylim = quantile(df_toPlot$readsCount, c(0.01, 0.99)))
      # + scale_fill_manual(breaks=c("WT", "KD"),values=c("grey", "red"))
  pdf(paste0(outdir, "BOXPLOT_",readQuantif ,"_reads_WTvs",Cond,"_",info,".pdf"), width = 16, height = 10 )
  print(p)
  dev.off()

}


boxplot_split = function(quantif_1, quantif_2 ,Nquantif_1, Nquantif_2, outdir, Nsplit = 5, info = NULL){
	# Split sorting vector
	quantif_2S = split(rownames(quantif_2), ceiling(seq_along(rownames(quantif_2))/ceiling(length(rownames(quantif_2))/Nsplit)))
	## GGplot method -> DATA FRAME
	df_toPlot = data.frame(matrix(ncol = 2, nrow = 0))
	for(i in 1:length(quantif_2S)){
		df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantif_1[rownames(quantif_1) %in% quantif_2S[[i]]])), rep(i, length(quantif_1[rownames(quantif_1) %in% quantif_2S[[i]]])))))
	}
	colnames(df_toPlot) = c("readsCount", "select")
	df_toPlot$select = as.factor(df_toPlot$select)
	df_toPlot$readsCount = log(as.numeric(as.character(df_toPlot$readsCount))+1)
  ## BOXPLOT
  p = ggplot(df_toPlot, aes(x=select, y=readsCount, color=select)) +
			geom_boxplot(outlier.shape = NA) +
			coord_cartesian(ylim = quantile(df_toPlot$readsCount, c(0.05, 0.95))) +
      labs(x=Nquantif_2, y = Nquantif_1)
      # scale_color_grey() + theme_classic() +
	pdf(paste0(outdir, "BOXPLOT_",Nquantif_1 ,"_splittedby_",Nquantif_2,"_",Nsplit,"tile_",info,".pdf"), width = 16, height = 10 )
	print(p)
	dev.off()
}





    # ns: p > 0.05
    #
    # *: p <= 0.05
    #
    # **: p <= 0.01
    #
    # ***: p <= 0.001
    #
    # ****: p <= 0.0001
