################################
# getProfmatFromBW  function
################################
# source("/home/robel/mount/NAS_CUVIER/grpCuvier/PROJETS/Cuvier/Projet_helicase/Figures/Mtr4/Figure2/Heatmap_profiles/Script_HEATMAP_profile.R")

# examples:
# profmat_Other_IBPs_under_all_peaks_Beaf_CP190_SMC3_smooth_csawNorm.lst = getProfmatFromBW(
#   bw.lst = list(Beaf_KD_M_L = Beaf_KD_Log2_M_L.bw,
#                 Beaf_KD_L = Beaf_KD_Luc.bw,
#                 Beaf_KD_M = Beaf_KD_Mtr4.bw,
#                 CP190_KD_M_L = CP190_KD_Log2_M_L.bw,
#                 CP190_KD_L = CP190_KD_Luc.bw,
#                 CP190_KD_M = CP190_KD_Mtr4.bw,
#                 SMC3_KD_M_L = SMC3_KD_Log2_M_L.bw,
#                 SMC3_KD_L = SMC3_KD_Luc.bw,
#                 SMC3_KD_M = SMC3_KD_Mtr4.bw,
#                 dCTCF_Rep1_WT_RPGC = dCTCF_Rep1_WT_RPGC.bw,
#                 M1BP_WT_Rep1 =  M1BP_WT_Rep1.bw,
#                 M1BP_WT_Rep2 =  M1BP_WT_Rep2.bw,
#                 GAF_WT_Rep1 =  GAF_WT_Rep1.bw,
#                 GAF_WT_Rep2 =  GAF_WT_Rep2.bw,
#                 Pita_WT_Rep1 =  Pita_WT_Rep1.bw,
#                 ZIPIC_WT_Rep1 =  ZIPIC_WT_Rep1.bw,
#                 SuHw_WT_Rep1 =  SuHw_WT_Rep1.bw,
#                 CG8436_WT_Rep1 =  CG8436_WT_Rep1.bw,
#                 CG9740_WT_Rep1 =  CG9740_WT_Rep1.bw,
#                 Mod_mdg4_Rep1 = Mod_mdg4_Rep1.bw,
#                 Chromator_WT_Rep1 = Chromator_WT_Rep1.bw,  
#                 Chromator_WT_Rep2 = Chromator_WT_Rep2.bw,  
#                 Z4_WT_Rep1 = Z4_WT_Rep1.bw),
#   bed.path = "/home/robel/mount/NAS_CUVIER/grpCuvier/PROJETS/Cuvier/Projet_helicase/data/merged_ChIP-Seq/allPeaks_BCS_concat_ML_500_red.bed",
#   bin = 10, xlim = c(5000, 5000), type = 'mf', anchor=NULL, smooth_spar=0.4)

# To Perform kmeans clustering or other clustering
# pdf("~/mount/NAS_CUVIER/grpCuvier/PROJETS/Cuvier/Projet_helicase/data/merged_ChIP-Seq/Ibps_underAllPeaks_reduced_8k_reordered_customScaled_all.pdf",width = 40,height = 12)
# set.seed(12345)
# clusters_8_customScale_seeded3 = multiHeatMatrix_local2(profmat_Other_IBPs_under_all_peaks_Beaf_CP190_SMC3_smooth_csawNorm.lst[c(1:11,13,15:21,23)],xcoords = c(-5000,5000),clustfun = function(x)kmeans(x,centers = 8),winsorize = c(1,99),RangeToCluster = c(475,525),col = list(c("blue","white","red"),c("white","#37A7DB"),c("white","#37A7DB"),c("blue","white","red"),c("white","#004CCC"),c("white","#004CCC"),c("blue","white","red"),c("white","#00D1CA"),c("white","#00D1CA"),c("white","#6656B6"),c("white","#086016"),c("white","#086016"),c("white","#086016"),c("white","#086016"),c("white","#086016"),c("white","#086016"),c("white","#086016"),c("white","#086016"),c("white","#086016"),c("white","#086016")),rangesList = list(c(-1,2),c(0,30),c(0,30),c(-1,2),c(0,30),c(0,30),c(-1.5,1),c(0,15),c(0,15),c(0,15),c(0,30),c(0,50),c(0,9),c(0,6),c(0,30),c(0,13),c(0,13),c(0,30),c(0,8),c(0,20)),KmeansUsed=TRUE,order=TRUE)
# dev.off()

# without clustering, because you've already performed your kmeans and know the order of peaks

# groupcol = rainbow(8)
# groupcol[4]="#0D613A"
# pdf("~/mount/NAS_CUVIER/grpCuvier/PROJETS/Cuvier/Projet_helicase/data/merged_ChIP-Seq/Ibps_underAllPeaks_reduced_8k_reordered_customScaled_CBS_only_allChIPs.pdf",width = 40,height = 12)
# # set.seed(12345)
# clusters_8_customScale_seeded3 = multiHeatMatrix_local2(profmat_Other_IBPs_under_all_peaks_Beaf_CP190_SMC3_smooth_csawNorm.lst[c(1:11,13,15:21,23)],xcoords = c(-5000,5000),winsorize = c(1,99),RangeToCluster = c(475,525),col = list(c("blue","white","red"),c("white","#37A7DB"),c("white","#37A7DB"),c("blue","white","red"),c("white","#004CCC"),c("white","#004CCC"),c("blue","white","red"),c("white","#00D1CA"),c("white","#00D1CA"),c("white","#6656B6"),c("white","#086016"),c("white","#086016"),c("white","#086016"),c("white","#086016"),c("white","#086016"),c("white","#086016"),c("white","#086016"),c("white","#086016"),c("white","#086016"),c("white","#086016")),rangesList = list(c(-1.5,1.5),c(0,30),c(0,30),c(-1.5,1.5),c(0,30),c(0,30),c(-1.5,1.5),c(0,15),c(0,15),c(0,15),c(0,30),c(0,50),c(0,9),c(0,6),c(0,30),c(0,13),c(0,13),c(0,30),c(0,8),c(0,20)),group = factor(clusters_8_customScale_seeded$cluster),order = TRUE,group.col = groupcol)
# dev.off()

getProfmatFromBW = function(
    bw.lst,
    bed.path = NULL,
    bed.gr = NULL,
    gr.name=NULL,
    bin = 10,
    xlim = c(5000, 5000),
    type = 'pf',
    anchor=NULL,
    smooth_spar=NULL, # value between 0 and 1, 0.4 should be enough
    seqplots.tmp = NULL)
    {
        if(is.null(bed.gr) & is.null(bed.path)){
            stop("bed.gr (GRanges object) or bed.path required")
        }


        if(!is.null(bed.gr)){
            if(is.null(gr.name)){
                gr.name = 'bed.gr.seqplot'
            }
            if(is.null(seqplots.tmp)){seqplots.tmp=getwd()}
            bed.path <- toString(rtracklayer::export.bed(gr.name,paste0(seqplots.tmp,"/",gr.name,"_#",length(bed.gr),".bed")))

        }else if(!is.null(bed.path)){
            if(is.null(gr.name)){
                gr.name = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(bed.path))
            }
            bed.gr = rtracklayer::import.bed(bed.path)
        }

        if(is.null(names(bed.gr))){
            if(!is.null(bed.gr$name)){
                names(bed.gr) = bed.gr$name
            }else{
                names(bed.gr) = paste0('elem_', seq_along(bed.gr))
            }
        }

        gpsa.lst <- seqplots::getPlotSetArray(bw.lst,bed.path,refgenome=NA,bin = bin,ignore_strand = T,xmin = xlim[1],xmax=xlim[2],rm0 = F,type=type, xanchored=anchor)

        profmat.lst = list()
        for(i in seq_along(bw.lst)){
            profmat.lst[[i]] = gpsa.lst$data[[sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(bed.path))]][[names(bw.lst)[i]]]$heatmap
            rownames(profmat.lst[[i]]) = names(bed.gr)
        }
        names(profmat.lst) = names(bw.lst)

        # replace NA
        lapply(seq_along(profmat.lst), function(ndx){
            profmat.lst[[ndx]][is.na(profmat.lst[[ndx]])] <<- 0
            # profmat.lst[[ndx]][is.na(profmat.lst[[ndx]])] <- 0
        })

        # smooth profmat
        if(!is.null(smooth_spar)){
            profmat.lst = lapply(profmat.lst, function(profmat){
                t(apply(profmat, 1, function(X){smooth.spline(X, spar = smooth_spar)$y}))
            })
        }
        return(profmat.lst)

}


## récupéré de https://github.com/al2na/genomation/blob/master/R/plotMatrix.R

library(grid)

.heatLegendY<-function(min,max,cols,legend.name,main=TRUE,cex.legend=1,
                       cex.lab=1){

  # get value range as a vector of 100
  vals=seq(min,max,length.out=100)
  rng <- range(vals, na.rm = TRUE) # get min/mqx
  m <- (vals - min)/(max-min) # get normalized range
  rasta= rgb(colorRamp(cols)(m), maxColorValue= 255) # get color for each element of range

  grid.raster( rev(rasta), interpolate=FALSE,height = unit(1, "npc"),
               width=unit(1, "npc")) # make the legend
  # make legend ticks
  at = seq(0,1,length.out=5); label = seq(min,max,length.out=5)

  #make the axis of the legend
  grid.yaxis(at=at,label=formatC(label,digits=2,format="g"),main=main,
             edits=gEdit("labels", rot=90,hjust=0.5),
             gp=gpar(cex=cex.legend)) # make axis for legend
  my.x = -2
  if(main==FALSE)
    my.x=3.4
  grid.text(legend.name,rot=90,x=unit(my.x, "npc"),gp=gpar(cex=cex.lab))

}

.convertToColors <- function(mat,cols, RangeValue = NULL) {
  # Produce 'normalized' version of matrix, with values ranging from 0 to 1
  # rng <- range(mat, na.rm = TRUE)
  if(is.null(RangeValue)){
    rng = range(mat, na.rm = TRUE)
    }else{
    rng = c(RangeValue[1], RangeValue[2])
  }
  m <- (mat - rng[1])/(diff(rng))
  # Convert to a matrix of sRGB color strings
  #m2 <- m; class(m2) <- "character"
  m2<-matrix("transparent",ncol=ncol(m),nrow=nrow(m))
  m2[!is.na(m)] <- rgb(colorRamp(cols)(m[!is.na(m)]), maxColorValue = 255)
  #m2[is.na(m)] <- "transparent"
  return(m2)
}

.gridHeat<-function(mat,col,rng,xcoords,xlab,cex.lab,cex.axis,angle=0,
                    hjust=0,vjust=0){

  mat2=.convertToColors(mat,col,rng)
  ras=grid.raster(mat2,interpolate = FALSE, width= unit(1, "npc"),
                  height=unit(1, "npc"))

  # make legend ticks
  at = seq(0,1,length.out=5); label = paste0(seq(min(xcoords),max(xcoords),length.out=5)/1000,"Kb")

  ax=grid.xaxis(at=at,label=formatC(label,digits=4,format="g"),
                edits=gEdit("labels", rot=angle,hjust=hjust,vjust=vjust),
                gp=gpar(cex=cex.axis)) # make axis for legend

  grid.text(xlab,y=unit(-2.5, "lines"),gp=gpar(cex=cex.lab))
  #grid.draw(ax)
}

.rowSideCol<-function(group.vector,group.names=NULL,group.col=NULL,cex.lab=1){

  if( is.null(group.col) ){
    cols=rainbow(length(unique(group.vector)))
    img=cols[factor(group.vector,levels=unique(group.vector))]
  }else{
    img=group.col[factor(group.vector,levels=unique(group.vector))]
  }
  grid.raster(img,interpolate = FALSE, width= unit(1, "npc"),
              height=unit(1, "npc"))

  # segment heights calculated from group.vector
  # will be used to put group names in the middle of the segment
  segh=as.vector(table(factor(group.vector,levels=unique(group.vector))))
  name.coord=1-((cumsum(segh)-(segh/2))/sum(segh)) # NPC coord

  if( is.null(group.names)){
    grid.text(unique(group.vector), y=unit(name.coord,"npc"),
              x = unit(-0.5, "lines"),
              gp=gpar(cex=cex.lab),just="right")
  }else{
    grid.text(group.names, y=unit(name.coord,"npc"),
              x = unit(-0.5, "lines"),
              gp=gpar(cex=cex.lab),just="right")
  }

}

.heatLegendX<-function(min,max,cols,legend.name,main=TRUE,cex.legend=1,
                       cex.lab=1,hjust=0,vjust=0){

  # get value range as a vector of 100
  vals=seq(min,max,length.out=100)
  rng <- range(vals, na.rm = TRUE) # get min/mqx
  m <- (vals - min)/(max-min) # get normalized range
  rasta= rgb(colorRamp(cols)(m), maxColorValue = 255) # get color for each element of range

  grid.raster( matrix((rasta),nrow=1), interpolate=FALSE,height = unit(1, "npc"),
               width=unit(1, "npc")) # make the legend
  # This is making breaks at sites out of the range of min and max. This is modified because our scales are too low. look eg. down
  # pretty(c(-0.75,0.75),n=3)
  # [1] -1.0 -0.5  0.0  0.5  1.0
  # # make legend ticks
  label = pretty(c(min,max),n=5);at = seq(0,1,length.out=length(label));
  
#   # make legend ticks
#   at = seq(0,1,length.out=5); label = paste0(seq(min(xcoords),max(xcoords),length.out=5)/1000,"Kb")
  #make the axis of the legend
  grid.xaxis(at=at,label=label,main=main,
             edits=gEdit("labels",hjust=hjust,vjust=vjust),
             gp=gpar(cex=cex.legend)) # make axis for legend
  my.y = -3
  grid.text(legend.name,y=unit(my.y, "npc"),gp=gpar(cex=cex.lab))

}



.winsorize<-function(mat,rng){
  hi.th=quantile(mat,rng[2]/100,na.rm=TRUE)
  lo.th=quantile(mat,rng[1]/100,na.rm=TRUE)
  mat[mat>hi.th]=hi.th
  mat[mat<lo.th]=lo.th
  mat
}



heatMatrixMat=function (mat, grid = FALSE, col = NULL, xcoords = NULL, group = NULL,
    group.col = NULL, order = FALSE, RangeValue = NULL, winsorize = c(0, 100), kmeans = FALSE,
    k = 3, main = "", legend.name = NULL, cex.legend = 1, xlab = NULL,
    cex.main = 1, cex.lab = 1, cex.axis = 1, newpage = TRUE, colorSet = NULL)
{
require(gridBase)
if(is.null(colorSet)){
  .jets<-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  # .jets<-colorRampPalette(c("#7F0000" ,   "red", "#FF7F00","yellow","#7FFF7F","cyan","#007FFF","blue", "#00007F"))
}else{
  .jets<-colorRampPalette(colorSet)
}
    if (class(mat)[1] != "matrix") {
        stop("'mat' is not a matrix\n")
    }
    mat2 = mat
    group.vector = NULL
    group.names = NULL
    if(is.null(RangeValue)){
      if (winsorize[2] < 100 | winsorize[1] > 0) {
        hi.th = quantile(mat2, winsorize[2]/100, na.rm = TRUE)
        lo.th = quantile(mat2, winsorize[1]/100, na.rm = TRUE)
        mat2[mat2 > hi.th] = hi.th
        mat2[mat2 < lo.th] = lo.th
      }
    }else{

        mat2[mat2 > RangeValue[2]] = RangeValue[2]
        mat2[mat2 < RangeValue[1]] = RangeValue[1]
    }
    if (kmeans) {
        if (any(is.na(mat2))) {
            mat3 = impute.knn(mat2, k = 10, rowmax = 0.5, colmax = 0.8,
                maxp = 1500)$data
            clu = kmeans(mat3, c = k)
        }
        else {
            clu = kmeans(mat2, c = k)
        }
        group.vector = clu$cluster
        kcenters = clu$centers
        mat2 = mat2[order(group.vector), ]
        group.vector = group.vector[order(group.vector)]
        if (order) {
            g.factor = factor(group.vector, levels = unique(group.vector))
            cent.val = rowSums(kcenters, na.rm = TRUE)[g.factor]
            my.order = order(-cent.val, group.vector, -rowSums(mat2,
                na.rm = TRUE))
            mat2 = mat2[my.order, ]
            group.vector = group.vector[my.order]
        }
        group.names = unique(group.vector)
    }
    if (!is.null(group) & !kmeans) {
        if (is.list(group)) {
            win.numbs = (lapply(group, function(x) unique(x)))
            win.vec = unlist(win.numbs)
            if (any(table(unlist(win.numbs)) > 1)) {
                stop("'group' containing a list must not have duplicated numbers\n")
            }
            row.ids = rownames(mat2)
            group.vector = rep(0, length(row.ids))
            for (i in 1:length(win.numbs)) {
                group.vector[row.ids %in% win.numbs[[i]]] = i
            }
            if (!is.null(names(group))) {
                group.names = names(group)
            }
            if (all(group.vector == 0)) {
                stop("None of the elements in 'group' are a part of rownames(mat) \n")
            }
            if (any(group.vector == 0)) {
                warning("Number of elements in 'group' argument is less then nrow(mat) \n",
                  "Dropping rows from 'mat' that are not contained in 'group'\n")
                mat2 = mat2[group.vector > 0, ]
                group.vector = group.vector[group.vector > 0]
            }
        }
        else if (is.factor(group)) {
            if (length(group) != nrow(mat)) {
                stop("'group' is a factor, and its length should be equal to nrow(mat)\n")
            }
            group = factor(as.character(group), levels = as.character(unique(group)))
            group.names = levels(group)
            levels(group) = 1:length(levels(group))
            group.vector = as.numeric(group)
        }
        else {
            stop("'group' must be a factor or a list\n")
        }
        mat2 = mat2[order(group.vector), ]
        group.vector = group.vector[order(group.vector)]
        if (order) {
            my.order = order(group.vector, -rowSums(mat2, na.rm = TRUE))
            mat2 = mat2[my.order, ]
            group.vector = group.vector[my.order]
            names(group.vector) = rownames(mat2)
        }
    }
    else if (order & !kmeans) {
        order.vector = rep(1, nrow(mat2))
        names(order.vector) = rownames(mat2)
        order.vector = order.vector[order(-rowSums(mat2, na.rm = TRUE))]
        mat2 = mat2[order(-rowSums(mat2, na.rm = TRUE)), ]
    }
    if (!grid) {
        plot.new()
        vps <- baseViewports()
        pushViewport(vps$figure)
    }
    else {
        if (newpage)
            grid.newpage()
    }
    # X axis coordinates legend
    if (!is.null(xcoords) & is.vector(xcoords)) {
        if (length(xcoords) == 2 & xcoords[1] < xcoords[2]) {
            xcoords = seq(xcoords[1], xcoords[2], length.out = ncol(mat2))
        }
        if (length(xcoords) != ncol(mat2))
            stop("xcoords has wrong length: ", length(xcoords),
                " \n", " it should be equal to the number of columns of ScoreMatrix\n",
                " which is", ncol(mat2), "\n")
    }
    else {
        xcoords = 1:ncol(mat2)
    }
    if (is.null(col)) {
        col = .jets(100)
    }
    legendVp <- viewport(width = unit(0.7, "lines"), height = unit(0.4,
        "npc"), x = unit(0.71, "npc"), y = unit(0.5, "npc"),
        just = "left")
    pushViewport(legendVp)
    if(is.null(RangeValue)){
      rng = range(mat2, na.rm = TRUE)
      }else{
      rng = c(RangeValue[1], RangeValue[2])
    }
    .heatLegendY(min = rng[1], max = rng[2], col, legend.name,
        main = FALSE, cex.legend, cex.lab)
    popViewport()
    heatHeightNPC = 0.7
    heatVp <- viewport(width = unit(0.5, "npc"), height = unit(heatHeightNPC,
        "npc"), x = unit(0.2, "npc"), y = unit(0.5, "npc"), just = "left")
    pushViewport(heatVp)
    .gridHeat(mat2, col, RangeValue, xcoords, xlab, cex.lab, cex.axis, hjust = 0.5)
    popViewport()
    if (!is.null(group.vector)) {
        sideVp <- viewport(width = unit(0.05, "npc"), height = unit(heatHeightNPC,
            "npc"), x = unit(0.145, "npc"), y = unit(0.5, "npc"),
            just = "left")
        pushViewport(sideVp)
        grid.rect()
        .rowSideCol(group.vector, group.names = group.names,
            group.col = group.col, cex.lab = cex.lab)
        popViewport()
    }
    title.y = unit(0.9, "npc")
    grid.text(main, y = title.y, x = unit(0.45, "npc"), gp = gpar(cex = cex.main))
    if (!grid)
        popViewport()
    if (kmeans | !is.null(group.vector)) {
        return(invisible(group.vector))
    }
    else if (order & is.null(group.vector)) {
        return(invisible(order.vector))
    }
}
