#############################
#
# do PCA and plot dimX vs dimY then plot clustering of PCA on dim chosen
#
# return : plot
#############################
invisible(require(dendsort))
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
#
# dat = list(
# 	PAUSE_IND_pol2_ctrl = PAUSE_INDICE_VEC$PAUSE_IND_pol2_ctrl,
# 	LOG2RATIO_K36me3_me2_f = LIST_QUANTIF_K36$LOG2RATIO_K36me3_me2_f,
# 	SUM_DELTA_NUC_C1_L3_f = LIST_QUANTIF_MNASE$SUM_DELTA_NUC_C1_L3_f,
# 	Q_H3K36me3_2C4_GB_f = LIST_QUANTIF_K36$Q_H3K36me3_2C4_GB_f,
# 	Q_H3K36me2_2C3_GB_f = LIST_QUANTIF_K36$Q_H3K36me2_2C3_GB_f,
# 	ZSCORE_K27H_K27C_f = LIST_QUANTIF_K27$ZSCORE_K27H_K27C_f,
# 	ZSCORE_K27M_K27C_f = LIST_QUANTIF_K27$ZSCORE_K27M_K27C_f,
# 	ZSCORE_K27H_K27C_MNASE_f = LIST_QUANTIF_K27$ZSCORE_K27H_K27C_MNASE_f,
# 	ZSCORE_K27M_K27C_MNASE_f = LIST_QUANTIF_K27$ZSCORE_K27M_K27C_MNASE_f,
# 	ZSCORE_nuc_H1_nuc_C1_GB_f = LIST_QUANTIF_MNASE$ZSCORE_nuc_H1_nuc_C1_GB_f*-1,
# 	ZSCORE_nuc_M1_nuc_C1_GB_f = LIST_QUANTIF_MNASE$ZSCORE_nuc_M1_nuc_C1_GB_f*-1,
# 	ZSCORE_SUM_DELTA_H1_L3_C1_L3 = LIST_QUANTIF_MNASE$ZSCORE_SUM_DELTA_H1_L3_C1_L3*-1,
# 	ZSCORE_SUM_DELTA_M1_L3_C1_L3 = LIST_QUANTIF_MNASE$ZSCORE_SUM_DELTA_M1_L3_C1_L3*-1,
# 	LOGFC_HYPBKD_RPGC_pval005_f = LIST_QUANTIF_RNASEQ$LOGFC_HYPBKD_RPGC_pval005_f*-1,
# 	LOGFC_MES4KD_RPGC_pval005_fmoins1 = LIST_QUANTIF_RNASEQ$LOGFC_MES4KD_RPGC_pval005_f*-1,
# 	LOGFC_MES4KD_RPGC_pval005_f = LIST_QUANTIF_RNASEQ$LOGFC_MES4KD_RPGC_pval005_f,
# 	dist_gene2border_f = dist_gene2border_f*-1
# )
# list.clust=list(c(1,2),c(1,3),c(2,3))
# list.cpa=list(c(1,2),c(1,3),c(2,3))
# cpa.title="ACP"
# dist.method="euclidean"
# choix="var"
# col.ind = "1"

do_acp_and_clustering <- function(
                            dat,
                            ncp=length(dat),
                            list.clust=list(c(1,2),c(1,3),c(2,3)),
                            list.cpa=list(c(1,2),c(1,3),c(2,3)),
                            cpa.title="ACP", dist.method="pearson",
                            choix="ind",
                            col.ind = "1"){

  invisible(require(FactoMineR))
  invisible(require(dendextend))


  old.par <- par()
  if(class(dat) %in% "list"){
    comRownames = Reduce(intersect,lapply(dat,names))
    dat = lapply(dat, function(feat1){feat1 = feat1[names(feat1) %in% comRownames]})
    data_df = matrix(NA,ncol = length(dat), nrow = length(dat[[1]]))
    rownames(data_df) = names(dat[[1]])
    colnames(data_df) = names(dat)
    for(i in 1:length(dat)){
      data_df[,i] = dat[[i]][rownames(data_df)]
    }
    dat = data_df
  }

  stopifnot(choix %in% c("var", "ind"))



  ncp=ncp
  dat.PCA=PCA(X=dat, graph=F, ncp=ncp)
  eig_title=paste(c("eigenvalue for",cpa.title),collapse=" ")
  #   plot(dat.PCA$eig[,2], type='h', lwd=5, main=eig_title)
  eig <- dat.PCA$eig[1:ncp,"percentage of variance"]
  barplot(eig, names.arg=paste(round(eig, 2),"%"), las=2, main=eig_title)

  ncp <- max(c(unlist(list.clust), unlist(list.cpa)))
  ncp <- min(ncp, dim(dat)[2])
  for (i in list.cpa){
    cpa_title=paste(cpa.title,": dim",i[1],"vs",i[2])
    plot.PCA(dat.PCA, choix=choix, axe=i, title=cpa_title, col.ind = col.ind)
    # second plot without label
    plot.PCA(dat.PCA, choix=choix, axe=i, title=cpa_title, col.ind = col.ind, label="none", lwd =4)
  }

  ## vérifier que pour chaque ligne de dat.PCA[[2]][[1]][,i]
  ## les valeurs sont différentes (sinon la fonction plante)
  for (i in list.clust){
    #dd => matrice de distance
    # autre methode de distance qu pred en compte la projection su rl'acp
    if(dist.method == "euclidean"){
      dd <-  dist(dat.PCA[[choix]][["coord"]][,i], method = dist.method)
    }else{
      dd <- as.dist((1 - cor(t(dat.PCA[[choix]][["coord"]][,i]), method=dist.method))/2)
    }

    hc <- hclust(dd, method = "ward.D")
    clust_title <- paste("classification sur les donnees projetees par l'ACP :\naxes",paste(i, collapse=", "), paste0("  dist : ",dist.method))

    mar.l <- max(nchar(hc$labels))/2


    hc <- as.dendrogram(hc)
    labels_colors(hc) <- col.ind[order.dendrogram(hc)]
    par(mar= c(mar.l,3,3,1))

    hc <- sort_hclust(hc)
    plot(hc, main = clust_title)
    # Second plot without label

    plot(hc ,labels=FALSE, main = clust_title, lwd = 4)
  }
  par(mar=old.par$mar)
}
