

################################################################################
## Comparison by boxplot of 1 feature for a  list of genes subset with wilcoxon test
##

print('USAGE : boxplot_1COND_wilcoxListFilter(quantifWT, cond1 = "WT",filterGNList, effMin = NULL, YLIM = c(0, 20), outlierTH = 0.01, logTrans =T,comp_list=NULL, outdir, readQuantif = "readQuantif", select = "filterValue", info = "")')
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
boxplot_1COND_wilcoxListFilter = function(quantifWT, cond1 = "WT",filterGNList, effMin = NULL, YLIM = c(0, 20), outlierTH = 0.01, logTrans =T,comp_list=NULL, outdir, readQuantif = "readQuantif", select = "filterValue", info = ""){
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
    }
  }else{
    # df_toPlot = as.data.frame(cbind(as.numeric(unname(quantifWT)), rep("all_Genes", length(quantifWT)), rep(cond1, length(quantifWT))))
    # df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD)), rep("all_Genes", length(quantifKD)), rep(cond2, length(quantifKD)))))
    for(i in 1:length(filterGNList_effMin)){
      df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifWT[rownames(quantifWT) %in% filterGNList_effMin[[i]],1,drop=F])), rep(names(filterGNList_effMin)[i], length(quantifWT[rownames(quantifWT) %in% filterGNList_effMin[[i]],1,drop=F])), rep(cond1, length(quantifWT[rownames(quantifWT) %in% filterGNList_effMin[[i]],1,drop=F])))))
    }
  }
  colnames(df_toPlot) = c("readsCount", "select")
  if(logTrans %in% T){
    df_toPlot$readsCount = log(as.numeric(as.character(df_toPlot$readsCount))+1)
  }else{
    df_toPlot$readsCount = as.numeric(as.character(df_toPlot$readsCount))
  }
  ## ratio
  # df_ratio = data.frame()
  # for(i in 1:length(levels(df_toPlot[,2]))){
  #   df_ratioTEMP = df_toPlot[df_toPlot[,2] %in% levels(df_toPlot[,2])[i],,drop=F]
  #   df_ratioTEMP2 = cbind(as.numeric(as.character(df_ratioTEMP[df_ratioTEMP[,3] %in% levels(df_ratioTEMP[,3])[2],1]))-as.numeric(as.character(df_ratioTEMP[df_ratioTEMP[,3] %in% levels(df_ratioTEMP[,3])[1],1])),
  #   rep(levels(df_toPlot[,2])[i] ,nrow(df_ratioTEMP)/2), rep("RATIO" ,nrow(df_ratioTEMP)/2))
  #   df_ratio = rbind(df_ratio, df_ratioTEMP2)
  # }
  # df_ratio[,1] = as.numeric(as.character(df_ratio[,1]))
  #
  # colnames(df_ratio) = c("readsCount", "select", "condition")
  # means <- aggregate(readsCount ~  select, df_ratio, median)
  # print(df_toPlot)
  ## BOXPLOT
  # my_comparisons <- list( c("1", "2"), c("2", "3"), c("3", "4"))
  p1_outlier = ggplot(df_toPlot, aes(x=select, y=readsCount))+
      geom_boxplot(outlier.shape = 4, lwd=2) + stat_compare_means(method = "wilcox.test",comparisons = comp_list) + ggtitle("auto ylim, outliers plotted and ketp for wilcox")
  p1_NOoutlier = ggplot(df_toPlot, aes(x=select, y=readsCount))+
      geom_boxplot(outlier.shape = 4, lwd=2) + scale_y_continuous(limits = quantile(df_toPlot$readsCount, c(outlierTH, 1-outlierTH))) +
       stat_compare_means(method = "wilcox.test",comparisons = comp_list) + ggtitle(paste0("auto ylim, outlier removed : ",outlierTH*2*100, "%"))
  p1_fig = ggplot(df_toPlot, aes(x=select, y=readsCount))+
      geom_boxplot(outlier.shape = NA, lwd=2) +  ylim(YLIM) + ggtitle(paste0("figure style plot, ylim = ",YLIM)) + theme_classic(base_line_size = 1, base_rect_size=1)
  # + scale_fill_manual(breaks=c("WT", "KD"),values=c("grey", "red"))
  pdf(paste0(outdir, "BOXPLOT_",readQuantif ,"_by",select ,"_",info,".pdf"), height=10, width=20)
    print(p1_outlier)
    print(p1_NOoutlier)
    textplot(lapply(filterGNList, length))
    textplot(lapply(filterGNList_effMin, length))
    print(p1_fig)
  dev.off()
}



library("data.table")














#end
