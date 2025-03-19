
################################################################################
## Comparison by boxplot of reads quantif upon 2 conditions / 2 features for a list of genes subset with wilcoxon test
##

print('USAGE : boxplot_wilcoxListFilter(quantifWT, quantifKD, cond1 = "WT", cond2="KD", filterGNList, effMin = NULL, YLIM = c(0, 20), bxplt_color=NULL,outlierTH = 0.01, logTrans =T, outdir, readQuantif = "readQuantif", Cond = "KD", select = "", info = "")')
## > str(quantifWT)
## Named num [1:24156] 1143 4523 1525 2299 1700 ...
## - attr(*, "names")= chr [1:24156] "100038246" "10006" "100126326" "100126327" ...
#
# > str(filterGNList)
## List of 6
 # $ ALL_TSS      : chr [1:24156] "100038246" "10006" "100126326" "100126327" ...
 # $ NO_Notch     : chr [1:16580] "100038246" "100126326" "100126327" "100127889" ...
 # $ NOTCH        : chr [1:7576] "79854" "26155" "339451" "57801" ...
 # $ BMI1         : chr [1:640] "5293" "9249" "3399" "57822" ...
 # $ BMI1xNOTCH   : chr [1:220] "9372" "257194" "164045" "26191" ...
 # $ NOTCH_no_BMI1: chr [1:7356] "79854" "26155" "339451" "57801" ...
 # outlierTH = 0.01 means that min AND max 1% outliers are removed (so 2% in total)
#
################################################################################
# quantifWT = LIST_QUANTIF_K27$Q_K27C_RPGC_GB_f
# quantifKD = LIST_QUANTIF_K27$Q_K27H_RPGC_GB_f
# cond1 = "Q_K27C_RPGC_GB"
# cond2="Q_K27H_RPGC_GB"
# filterGNList = List_genes_DOM_withISLAND_withOrientationGene
# effMin = 150
# SampleNorm = c(T, "noUP_noDN")
# YLIM = c(1,5)
# bxplt_color = c("#3d3d3d", "#820002")
# outlierTH = 0.01
# logTrans =T
# outdir = paste0(workdir, "PROJET_K27K9K36/FIGURE/BOXPLOT/H3K27me3/DELTA_H3K27me3/GB_by_Clust7G_DEA/")
# readQuantif = "Q_K27C_RPGC_GB"
# Cond = "Q_K27H_RPGC_GB"
# select = "DEA_Cluster_7G"
# info = NULL
################################################################################
boxplot_wilcoxListFilter = function(quantifWT, quantifKD, cond1 = "WT", cond2="KD", filterGNList, effMin = NULL, test.side = "two.sided", SampleNorm = c(F, "NULL"), YLIM = c(0, 20), bxplt_color=NULL, outlierTH = 0.01, logTrans =T, outdir, readQuantif = "readQuantif", Cond = "KD", select = "filterValue", info = ""){
  # LOG transform
  if(logTrans %in% T){
        quantifWT = log2(quantifWT+1)
        quantifWT=quantifWT[!is.na(quantifWT)]
        quantifKD = log2(quantifKD+1)
        quantifKD=quantifKD[!is.na(quantifKD)]
  }
  #d <- d[!is.na(d)]
  ## remove OUTLIER
  if(outlierTH %ni% 0){
    for(i in 1:length(filterGNList)){
      QUANT_WT_fortrim = quantifWT[names(quantifWT) %in% filterGNList[[i]]]
      QUANT_KD_fortrim = quantifKD[names(quantifKD) %in% filterGNList[[i]]]
      GN_to_remove = unique(c(names(which(QUANT_WT_fortrim > quantile(QUANT_WT_fortrim, 1-outlierTH))), names(which(QUANT_WT_fortrim < quantile(QUANT_WT_fortrim, outlierTH))),
      names(which(QUANT_KD_fortrim > quantile(QUANT_KD_fortrim, 1-outlierTH))), names(which(QUANT_KD_fortrim < quantile(QUANT_KD_fortrim, outlierTH)))))
      filterGNList[[i]] = filterGNList[[i]][filterGNList[[i]] %ni% GN_to_remove]
    }
  }
  ## Norm inter sample
  if(SampleNorm[1] %in% T){
    GN_for_norm = filterGNList[[SampleNorm[2]]]
    SampleNorm_fact = median(quantifKD[GN_for_norm])/median(quantifWT[GN_for_norm])
    quantifWT = quantifWT*SampleNorm_fact
  }

  ### Normalisation des effectifs
  if(is.null(effMin)){
    effMin =  min(unlist(lapply(filterGNList, length)))
  }
  filterGNList_effMin = lapply(filterGNList, function(filterGN){
    if(length(filterGN) >= effMin){
      filterGN = sample(filterGN, effMin)
    }else{
      filterGN = filterGN
    }
  })
  ## GGplot method -> DATA FRAME
  df_toPlot = data.frame()
  if(class(quantifWT) %in% "numeric"){
    # df_toPlot = as.data.frame(cbind(as.numeric(unname(quantifWT)), rep("all_Genes", length(quantifWT)), rep(cond1, length(quantifWT))))
    # df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD)), rep("all_Genes", length(quantifKD)), rep(cond2, length(quantifKD)))))
    for(i in 1:length(filterGNList_effMin)){
      df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifWT[names(quantifWT) %in% filterGNList_effMin[[i]]])), rep(names(filterGNList_effMin)[i], length(quantifWT[names(quantifWT) %in% filterGNList_effMin[[i]]])), rep(cond1, length(quantifWT[names(quantifWT) %in% filterGNList_effMin[[i]]])))))
      df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD[names(quantifKD) %in% filterGNList_effMin[[i]]])), rep(names(filterGNList_effMin)[i], length(quantifKD[names(quantifKD) %in% filterGNList_effMin[[i]]])), rep(cond2, length(quantifKD[names(quantifKD) %in% filterGNList_effMin[[i]]])))))
    }
  }else{
    # df_toPlot = as.data.frame(cbind(as.numeric(unname(quantifWT)), rep("all_Genes", length(quantifWT)), rep(cond1, length(quantifWT))))
    # df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD)), rep("all_Genes", length(quantifKD)), rep(cond2, length(quantifKD)))))
    for(i in 1:length(filterGNList_effMin)){
      df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifWT[rownames(quantifWT) %in% filterGNList_effMin[[i]],1,drop=F])), rep(names(filterGNList_effMin)[i], length(quantifWT[rownames(quantifWT) %in% filterGNList_effMin[[i]],1,drop=F])), rep(cond1, length(quantifWT[rownames(quantifWT) %in% filterGNList_effMin[[i]],1,drop=F])))))
      df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD[rownames(quantifKD) %in% filterGNList_effMin[[i]],1,drop=F])), rep(names(filterGNList_effMin)[i], length(quantifKD[rownames(quantifKD) %in% filterGNList_effMin[[i]],1,drop=F])), rep(cond2, length(quantifKD[rownames(quantifKD) %in% filterGNList_effMin[[i]],1,drop=F])))))
    }
  }
  colnames(df_toPlot) = c("readsCount", "select", "condition")
  df_toPlot$readsCount = as.numeric(as.character(df_toPlot$readsCount))
  ## ratio
  df_ratio = data.frame()
  for(i in 1:length(levels(df_toPlot[,2]))){
    df_ratioTEMP = df_toPlot[df_toPlot[,2] %in% levels(df_toPlot[,2])[i],,drop=F]
    df_ratioTEMP2 = cbind(as.numeric(as.character(df_ratioTEMP[df_ratioTEMP[,3] %in% levels(df_ratioTEMP[,3])[2],1]))-as.numeric(as.character(df_ratioTEMP[df_ratioTEMP[,3] %in% levels(df_ratioTEMP[,3])[1],1])),
    rep(levels(df_toPlot[,2])[i] ,nrow(df_ratioTEMP)/2), rep("RATIO" ,nrow(df_ratioTEMP)/2))
    df_ratio = rbind(df_ratio, df_ratioTEMP2)
  }
  df_ratio[,1] = as.numeric(as.character(df_ratio[,1]))

  colnames(df_ratio) = c("readsCount", "select", "condition")
  means <- aggregate(readsCount ~  select, df_ratio, median)


  # bxplt_color
  if(is.null(bxplt_color)){
    bxplt_color =  c("#aaaaaa", "#5e5e5e")
  }
  # print(df_toPlot)
  ## BOXPLOT
  p1_wilcox_paired = ggplot(df_toPlot, aes(x=select, y=readsCount, fill=condition))+theme_classic(base_line_size = 1, base_rect_size=1)+
      geom_boxplot(outlier.shape = 4, lwd=1.5, colour="#000000") + stat_compare_means(method = "wilcox.test", paired=T, method.args = list(alternative = test.side)) + ggtitle("auto ylim, outliers plotted and ketp for wilcox paired") + scale_fill_manual(values=bxplt_color)
  p1_wilcox = ggplot(df_toPlot, aes(x=select, y=readsCount, fill=condition))+theme_classic(base_line_size = 1, base_rect_size=1)+
      geom_boxplot(outlier.shape = 4, lwd=1.5, colour="#000000") + stat_compare_means(method = "wilcox.test", paired=F, method.args = list(alternative = test.side)) + ggtitle("auto ylim, outliers plotted and ketp for wilcox NO paired") + scale_fill_manual(values=bxplt_color)
  # p1_ttest_paired = ggplot(df_toPlot, aes(x=select, y=readsCount, fill=condition))+theme_classic(base_line_size = 1, base_rect_size=1)+
  #     geom_boxplot(outlier.shape = 4, lwd=1.5, colour="#000000") + stat_compare_means(method = "t.test", paired=T, method.args = list(alternative = test.side)) + ggtitle("auto ylim, outliers plotted and ketp for t.test paired ") + scale_fill_manual(values=bxplt_color)
  # p1_ttest = ggplot(df_toPlot, aes(x=select, y=readsCount, fill=condition))+theme_classic(base_line_size = 1, base_rect_size=1)+
  #     geom_boxplot(outlier.shape = 4, lwd=1.5, colour="#000000") + stat_compare_means(method = "t.test", paired=F, method.args = list(alternative = test.side)) + ggtitle("auto ylim, outliers plotted and ketp for wilcox NO paired") + scale_fill_manual(values=bxplt_color)
  ## BOXPLOT NO OUTLIER
  p1outlier_wilcox_paired = ggplot(df_toPlot, aes(x=select, y=readsCount, fill=condition))+theme_classic(base_line_size = 1, base_rect_size=1)+
      geom_boxplot(outlier.shape = 4, lwd=1.5, colour="#000000") + scale_y_continuous(limits = quantile(df_toPlot$readsCount, c(outlierTH, 1-outlierTH))) +
      stat_compare_means(method = "wilcox.test", paired=T, method.args = list(alternative = test.side)) + ggtitle(paste0("auto ylim, outlier removed : ",outlierTH*2*100, "% for wilcox paired")) + scale_fill_manual(values=bxplt_color)
  p1outlier_wilcox = ggplot(df_toPlot, aes(x=select, y=readsCount, fill=condition))+theme_classic(base_line_size = 1, base_rect_size=1)+
      geom_boxplot(outlier.shape = 4, lwd=1.5, colour="#000000") + scale_y_continuous(limits = quantile(df_toPlot$readsCount, c(outlierTH, 1-outlierTH))) +
      stat_compare_means(method = "wilcox.test", paired=F, method.args = list(alternative = test.side)) + ggtitle(paste0("auto ylim, outlier removed : ",outlierTH*2*100, "% for wilcox NO paired")) + scale_fill_manual(values=bxplt_color)
  # p1outlier_ttest_paired = ggplot(df_toPlot, aes(x=select, y=readsCount, fill=condition))+theme_classic(base_line_size = 1, base_rect_size=1)+
  #     geom_boxplot(outlier.shape = 4, lwd=1.5, colour="#000000") + scale_y_continuous(limits = quantile(df_toPlot$readsCount, c(outlierTH, 1-outlierTH))) +
  #     stat_compare_means(method = "t.test", paired=T, method.args = list(alternative = test.side)) + ggtitle(paste0("auto ylim, outlier removed : ",outlierTH*2*100, "% for t.test paired ")) + scale_fill_manual(values=bxplt_color)
  # p1outlier_ttest = ggplot(df_toPlot, aes(x=select, y=readsCount, fill=condition))+theme_classic(base_line_size = 1, base_rect_size=1)+
  #     geom_boxplot(outlier.shape = 4, lwd=1.5, colour="#000000") + scale_y_continuous(limits = quantile(df_toPlot$readsCount, c(outlierTH, 1-outlierTH))) +
  #     stat_compare_means(method = "t.test", paired=F, method.args = list(alternative = test.side)) + ggtitle(paste0("auto ylim, outlier removed : ",outlierTH*2*100, "% for wilcox NO paired")) + scale_fill_manual(values=bxplt_color)


  # p1_NOoutlier = ggplot(df_toPlot, aes(x=select, y=readsCount, fill=condition))+
  #     geom_boxplot(outlier.shape = 4, lwd=1.5, colour="#000000") + scale_y_continuous(limits = quantile(df_toPlot$readsCount, c(outlierTH, 1-outlierTH))) +
  #      stat_compare_means(method = "wilcox.test") + ggtitle(paste0("auto ylim, outlier removed : ",outlierTH*2*100, "%")) + scale_fill_manual(values=bxplt_color) + theme_classic(base_line_size = 1, base_rect_size=1)
  ## BOXPLOT of RATIO
  p2 = ggplot(df_ratio, aes(x=select, y=readsCount)) + geom_boxplot(outlier.shape = NA, lwd=1.5, colour="#000000") + stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3,show.legend = FALSE) +
  geom_text(data = means, aes(label = round(readsCount,4), y = readsCount + 0.5)) +  ylim(-2,2) + scale_fill_manual(values=bxplt_color)
  ## BOXPLOT FOR FIGURE
  p1_fig = ggplot(df_toPlot, aes(x=select, y=readsCount, fill=condition))+
      geom_boxplot(outlier.shape = NA, lwd=1.5, colour="#000000") +  ylim(YLIM) + ggtitle(paste0("figure style plot, ylim = ",YLIM)) + scale_fill_manual(values=bxplt_color)  + theme_classic(base_line_size = 1, base_rect_size=1)
  p2_fig = ggplot(df_toPlot, aes(x=select, y=readsCount, fill=condition))+
      geom_violin(scale = "width", trim = F, adjust = .5, lwd=1.5, colour="#000000") +  ylim(c(YLIM[1],YLIM[2]+2)) + ggtitle(paste0("figure style plot, ylim = ",YLIM)) + scale_fill_manual(values=bxplt_color)  + theme_classic(base_line_size = 1, base_rect_size=1)
  pdf(paste0(outdir, "BOXPLOT_",readQuantif ,"_reads_WTvs",Cond, "_by",select ,"_",info,".pdf"), height=10, width=20)
    print(p1_wilcox_paired)
    print(p1_wilcox)
    # print(p1_ttest_paired)
    # print(p1_ttest)
    print(p1outlier_wilcox_paired)
    print(p1outlier_wilcox)
    # print(p1outlier_ttest_paired)
    # print(p1outlier_ttest)
    textplot(lapply(filterGNList, length))
    textplot(lapply(filterGNList_effMin, length))
    print(p2)
    print(p1_fig)
    print(p2_fig)
  dev.off()
}
