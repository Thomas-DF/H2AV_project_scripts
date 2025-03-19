# Nuc Positionning prof
# Cuvier''s Team
# Depierre*
# 2018 16 avril

# plot HISTO with colored tail corresponding to a given percent

print('USAGE : plotHISTO_coltail(HIST_data = vec(), THqtl = 0.05, VarName = "name", info = "info", outDir = paste0(workdir, "David_M2_K9K27DE/integration_DE_DHMG/HISTOGRAM/"))')

plotHISTO_coltail = function(HIST_data, THqtl = 0.05, VarName = NULL, info = NULL, outDir = paste0(workdir, "David_M2_K9K27DE/integration_DE_DHMG/HISTOGRAM/")){
  histgm <- hist(HIST_data, breaks=100, plot=FALSE)
  mycol = ifelse(histgm$breaks < quantile(HIST_data, THqtl), "red", ifelse (histgm$breaks > quantile(HIST_data, 1-THqtl), "blue", "gray50"))
  v1 = round(quantile(HIST_data, THqtl), 2)
  v2 = round(quantile(HIST_data, 1-THqtl), 2)
  pdf(paste0(outDir, "HISTOGRAM_",VarName ,"_", THqtl,"tail_",info,".pdf"), width = 16, height = 10 )
  plot(histgm, col = mycol, border=T, main=paste0("Histogram of ", VarName, " / red > ",v1," / blue < ",v2," / tail=",THqtl))
  dev.off()
}




##################################################################################
## FUNCTION TO COMPARE DISTRIBUTION of a same feature on 2
##################################################################################
# vec = vec
# > str(vec)
#  Named int [1:6967] 84243 65991 65346 62357 61081 60984 60953 57504 56606 55403 ...
#  - attr(*, "names")= chr [1:6967] "FBgn0029123.1" "FBgn0003866.1" "FBgn0030309.1" "FBgn0030310.1" ...
# nameVec = "nameVec"
# vec_GN = vec_GN
# > str(vec_GN) {liste d'interet}
#  chr [1:697] "FBgn0000566.1" "FBgn0027655.1" "FBgn0038067.1" ...
# nameVec_GN = "nameVec_GN"
# xlim = c(0, 30000)
# fillingcol = c("#3d3d3d", "#820002"),
# outDir = paste0(workdir, "path/to/out/")
# info = NULL
##################################################################################
library("data.table")
library(ggplot2)

print('USAGE : plotHisto_listGN(vec = vec, nameVec = "nameVec", vec_GN = vec_GN, nameVec_GN = "nameVec_GN", xlim = c(0, 30000), fillingcol = c("#3d3d3d", "#820002"),
outDir = paste0(workdir, "path/to/out/"), info = NULL)')

plotHisto_listGN = function(vec = vec, nameVec = "nameVec", vec_GN = vec_GN, nameVec_GN = "nameVec_GN", xlim = c(0, 30000), binW = 0.5, alphaTransp = 0.6, fillingcol = c("#3d3d3d", "#820002"),
outDir = paste0(workdir, "path/to/out/"), info = NULL){
  data2GG = data.table(vec)
  rownames(data2GG) = names(vec)
  data2GG$feat1 = "A"
  data2GG$feat1[rownames(data2GG) %in% vec_GN] = "B"
  data2GG$feat1 = as.factor(data2GG$feat1)
  pdf(paste0(outDir, "HISTOGRAM_",nameVec ,"_Spllited_", nameVec_GN,"_",info,".pdf"), width = 16, height = 10 )
  p_hist = ggplot(data=data2GG, aes(data2GG$vec, fill = data2GG$feat1, stat(density))) + xlim(xlim) + xlab(nameVec) +
  geom_histogram(binwidth=binW,alpha=alphaTransp, position="identity")  +
  scale_fill_manual(values = fillingcol,
    name=nameVec_GN,
    labels=c("NO", "YES")) +
    theme_classic(base_line_size = 1, base_rect_size=1) +
    geom_vline(data=data2GG, aes(xintercept=0), linetype="dashed", size=1, colour="#3d3d3d")
  p_hist2 = ggplot(data=data2GG, aes(data2GG$vec, fill = data2GG$feat1, stat(density))) + xlim(xlim) + xlab(nameVec) +
  geom_histogram(binwidth=binW,alpha=1, position="dodge")  +
  scale_fill_manual(values = fillingcol,
    name=nameVec_GN,
    labels=c("NO", "YES")) +
    theme_classic(base_line_size = 1, base_rect_size=1)
    print(p_hist)
    print(p_hist2)

    dev.off()
}

##################################################################################
## FUNCTION TO COMPARE DISTRIBUTION of a same feature on 2
##################################################################################
# vec = LIST_QUANTIF_K27$ZSCORE_K27M_K27C_f
# str(LIST_QUANTIF_K27$ZSCORE_K27M_K27C_f)
#  # Named int [1:6967] 84243 65991 65346 62357 61081 60984 60953 57504 56606 55403 ...
#  # - attr(*, "names")= chr [1:6967] "FBgn0029123.1" "FBgn0003866.1" "FBgn0030309.1" "FBgn0030310.1" ...
# nameVec = "ZSCORE_K27M_K27C_f"
# list_GN = List_genes_DOM_Vreduce2
# str(list_GN) #{liste d'interet}
# # List of 5
# #  $ OUTDOM   : chr [1:5267] "FBgn0000017.1" "FBgn0000018.1" "FBgn0000042.1" "FBgn0000043.1" ...
# #  $ OUTBORDER: chr [1:1495] "FBgn0000032.1" "FBgn0000053.1" "FBgn0000064.1" "FBgn0000092.1" ...
# #  $ INBORDER : chr [1:1577] "FBgn0000024.1" "FBgn0000028.1" "FBgn0000046.1" "FBgn0000075.1" ...
# #  $ INDOM    : chr [1:2806] "FBgn0000014.1" "FBgn0000015.1" "FBgn0000037.1" "FBgn0000038.1" ...
# #  $ ISLAND   : chr [1:702] "FBgn0000153.1" "FBgn0000183.1" "FBgn0000241.1" "FBgn0000250.1" ...
# nameVec_GN = "List_genes_DOM_Vreduce2"
# xlim = c(-20,20)
# fillingcol = c("#3d3d3d","#00660a","#54d426","#d46e26","#820002","#8c1aa1")
# outDir = paste0(workdir, "path/to/out/")
# info = NULL
##################################################################################

print('USAGE : plotHisto_MULTIlistGN(vec = vec, nameVec = "nameVec", list_GN = list_GN, nameVec_GN = "nameVec_GN", xlim = c(-20,20),binW = 0.5, alphaTransp = 0.6, fillingcol = c("#3d3d3d","#00660a","#54d426","#d46e26","#820002","#8c1aa1"),
outDir = paste0(workdir, "path/to/out/"), info = NULL)')


plotHisto_MULTIlistGN = function(vec = vec, nameVec = "nameVec", list_GN = list_GN, nameVec_GN = "nameVec_GN", xlim = c(-20,20), binW = 0.5, alphaTransp = 0.6, fillingcol = c("#3d3d3d","#00660a","#54d426","#d46e26","#820002","#8c1aa1"),
outDir = paste0(workdir, "path/to/out/"), info = NULL){
  data2GG = data.table(vec)
  rownames(data2GG) = names(vec)
  data2GG$feat1 = "feat1"
  for(i in 1:length(list_GN)){
    data2GG$feat1[rownames(data2GG) %in% list_GN[[i]]] = names(list_GN)[i]
  }
  data2GG$feat1 = as.factor(data2GG$feat1)
  p_hist = ggplot(data=data2GG, aes(data2GG$vec, fill = data2GG$feat1, stat(density))) + xlim(xlim) + xlab(nameVec) +
  geom_histogram(binwidth=binW,alpha=alphaTransp, position="identity")  +
  scale_fill_manual(values = fillingcol,
    name=nameVec_GN,
    labels=c("OTHERS", names(list_GN))) +
    theme_classic(base_line_size = 1, base_rect_size=1)
  p_hist2 = ggplot(data=data2GG, aes(data2GG$vec, fill = data2GG$feat1, stat(density))) + xlim(xlim) + xlab(nameVec) +
  geom_histogram(binwidth=binW,alpha=1, position="dodge")  +
  scale_fill_manual(values = fillingcol,
    name=nameVec_GN,
    labels=c("OTHERS", names(list_GN))) +
    theme_classic(base_line_size = 1, base_rect_size=1)
  pdf(paste0(outDir, "HISTOGRAM_",nameVec ,"_Spllited_", nameVec_GN,"_",info,".pdf"), width = 16, height = 10 )
    print(p_hist)
    print(p_hist2)
    textplot(table(data2GG$feat1))
    dev.off()
}

















#end
