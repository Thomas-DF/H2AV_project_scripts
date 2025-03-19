

################################################################################
## Comparison by boxplot of 1 feature for a  list of genes subset with wilcoxon test
##

print('USAGE : boxplot_1COND_wilcoxListFilter(quantifWT, cond1 = "WT",filterGNList, effMin = NULL,bxplt_color=NULL,YLIM = c(0, 20), outlierTH = 0.01, logTrans =T,comp_list=NULL, outdir, readQuantif = "readQuantif", select = "filterValue", info = "")')
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
boxplot_1COND_wilcoxListFilter_REF = function(quantifWT, cond1 = "WT",filterGNList, effMin = NULL,bxplt_color=NULL, YLIM = c(0, 20), outlierTH = 0.01, logTrans =T,comp_list=NULL, outdir, readQuantif = "readQuantif", select = "filterValue", info = "")
{
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
  colnames(df_toPlot) = c("readsCount", "select","condition")
  if(logTrans %in% T){
    df_toPlot$readsCount = log(as.numeric(as.character(df_toPlot$readsCount))+1)
  }else{
    df_toPlot$readsCount = as.numeric(as.character(df_toPlot$readsCount))
  }

if(is.null(bxplt_color)){
    bxplt_color =  c("#aaaaaa", "#5e5e5e")
  }


y_position_val= max(df_toPlot$readsCount)+1
y_position_val= trunc(y_position_val)
#y_position_vect=rep(y_position_val,length.out= length(filterGNList))

readsCount=df_toPlot$readsCount
condition=df_toPlot$condition
select= df_toPlot$select


### Paired wilcox_test
stat.test_1 <- df_toPlot %>%
group_by(condition) %>%
wilcox_test(readsCount ~ select,paired = TRUE) %>% # paire test 
adjust_pvalue(method = "bonferroni") %>%
add_significance("p")
stat.test_1



##### NO_ paired wilcox_test
stat.test_2 <- df_toPlot %>%
group_by(condition) %>%
wilcox_test(readsCount ~ select,paired = FALSE) %>% # NO paire test 
adjust_pvalue(method = "bonferroni") %>%
add_significance("p")
stat.test_2


pdf(paste0(outdir, "BOXPLOT_",readQuantif ,"_reads",".pdf"), height=7, width=10)
# Créer  boxplot 1 
bxp_1 <- ggboxplot(df_toPlot, x = "select", y = "readsCount", ylim=NULL, color = "black", fill="condition", palette = bxplt_color,title="wilcox paire test")
# Ajoutez des p-values sur les graphiques en box plot
stat.test_1 <- stat.test_1 %>% add_xy_position(x = "select", dodge = 0.8)
# ajouter les crochets
bxp_1 + stat_pvalue_manual(stat.test_1,  label = "{p}{p.signif}", tip.length = 0,remove.bracket = FALSE, y.position=y_position_val)

# Créer box plot 2 
bxp_2 <- ggboxplot(df_toPlot, x = "select", y = "readsCount", color = "black", fill="condition",  palette = bxplt_color,
title= " wilcox NO paire test")
# Ajoutez des p-values sur les graphiques en box plot
stat.test_2 <- stat.test_2 %>% add_xy_position(x = "select", dodge = 0.8)
# ajouter les crochets
bxp_2 + stat_pvalue_manual(stat.test_2,  label = "{p}{p.signif}", tip.length = 0,remove.bracket = FALSE, y.position=y_position_val)
textplot(lapply(filterGNList, length))
textplot(lapply(filterGNList_effMin, length))

dev.off()
}

