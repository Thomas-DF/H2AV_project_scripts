# REFKA 2021

#Ce package contient les principalaux fonctions utilisées habituellement 
#ensemble pour l'analyse des données CHIP-seq et qui doit étre analysé dans la meme figure 
# à savoir :HEATMAPS+AVEREGE+BOXPLOT associé 
#l'objectif est de ressembler dans le même PDF une visualisation + test statistique des données

#les fonctions ajoutées dans ce document sont ;
workdir="/home/cperrois/work/"
#source(paste0(workdir, "functionR/Script_HEATMAP_profile.R"))
#source(paste0(workdir, "functionR/Boxplot_wilcoxListFilter_REF.R"))
#source(paste0(workdir, "functionR/AVG_PROFILE.R"))

# Ces PACKAGES doivenent installés dans la session 
.libPaths(c(.libPaths(),"/home/cperrois/work/Rpackages/R-3.4.3/"))
library(GenomicRanges)
library(ggplot2)
library(ggpubr) # magrittr
library(gplots)
'%ni%' = Negate('%in%')
require(BiocGenerics)
require(parallel)
library(Rsamtools)
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library(rtracklayer)
library(nucleR)
library(gridBase)

library(dplyr) # used for recover names(coverages)
library(gsubfn) # used in bam2coverage
library(GenomeInfoDb)
library(normr)
library("htmlwidgets",lib="/home/cperrois/work/Rpackages/R-3.4.3/")
library(seqplots)
library("BSgenome.Dmelanogaster.UCSC.dm6")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(devtools)
library("BSgenome")
ibrary(ggpubr)  # https://github.com/kassambara/ggpubr
library(rstatix)  # https://github.com/kassambara/rstatix
library("gplots")
#install.packages('textplot')
library(textplot)
library("ggplot2")
#######""
workdir="/home/cperrois/work/"
#source(paste0(workdir, "functionR/BOXPLOT_1condWilcox_genesList.R"))
#source(paste0(workdir, "functionR/BOXPLOT_2condWilcox_genesList.R"))




## Comment ça marche ???? 
## DATA for heat maps 
rangeheatmap = c(1:1000)
PRFOMAT_MATRIX=
qauntif_vector=
# DATA preparation pour AVEG plot 
tmp <- create(paste0(workdir,"PROJET_H2AV_2021/FIGURE/AVG_PLOT/AVG_H2AV_HYPB_KD/TMPgetPlotSetArray"))
tmp <- paste0(workdir,"PROJET_H2AV_2021/FIGURE/AVG_PLOT/AVG_H2AV_HYPB_KD/TMPgetPlotSetArray")
bwdir = paste0(workdir, "PROJET_H2AV_2021/DATA/BIGWIG/")
BW_1=paste0(bwdir,"BW_1.bw")
BW_2=paste0(bwdir,"BW_2.bw")


  GR1 #=> Grange of a gene liste of intrest 
  GR2 #=> Grange of a gene liste of intrest 
)
GR_list_toPlot=c(
"GR1",
"GR2"
)

# DATA for boxplot 
List_Pause_Elong = list(
  Paused = PAUSE_IND_pol2_ctrl_N_GN_splitted[[2]],
  Elong = PAUSE_IND_pol2_ctrl_N_GN_splitted[[18]]
)


pdf(paste0("XXXXXXXXXXXXXXXXX.pdf"), height=7, width=10))

heatMatrixMat(PRFOMAT_MATRIX[names(qauntif_vector),rangeheatmap],winsorize=c(5,95),RangeValue=c(-1,1),main="PRFOMAT_MATRIX",legend.name="qauntif_vector")        

for(GR in GR_list_toPlot){
  seqPlotSDoutliers_scaleFact(c(BW_1, BW_1),tmp,GR,c(-1,4),c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
}

Boxplot_wilcoxListFilter_REF
(
quantifWT = Q1, 
quantifKD = Q_2,
cond1 = "Q_POLII_TSS_ctrl_f",
cond2="Q_POLII_TSS_hypbKD_f",
filterGNList = List_Pause_Elong,
SampleNorm = c(F, "NULL"),
effMin = 200, YLIM = c(3,6), bxplt_color = c("#5A5E6B", "#1E7FCB"), outlierTH = 0.01, logTrans=T,
outdir = outfig, 
readQuantif = "Q_POLII_TSS_ctrl_f", 
Cond = "Q_POLII_TSS_hypbKD_f", 
select = "List_Pause_Elong",
 info = NULL)



dev.off()



############################################################################################################
#     																									   #
#										HEAMTMAPS														   #
#
############################################################################################################
## récupéré de https://github.com/al2na/genomation/blob/master/R/plotMatrix.R

library(grid)
.heatLegendY<-function(min,max,cols,legend.name,main=TRUE,cex.legend=10,
                       cex.lab=10){

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

.gridHeat<-function(mat,col, RangeValue,xcoords,xlab,cex.lab,cex.axis,angle=0,
                    hjust=0,vjust=0){

  mat2=.convertToColors(mat,col, RangeValue)
  ras=grid.raster(mat2,interpolate = FALSE, width= unit(1, "npc"),
                  height=unit(1, "npc"))

  # make legend ticks
  at = seq(0,1,length.out=5); label = seq(min(xcoords),max(xcoords)
                                          ,length.out=5)

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
  # make legend ticks
  label = pretty(c(min,max),n=5);at = seq(0,1,length.out=length(label));

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
    if (class(mat) != "matrix") {
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


############################################################################################################
#     																									   #
#										AVERAGE_PLOT														   #
#
############################################################################################################

#####################################################################################-
#          INSIDE FUNCTIONS  ----
#####################################################################################-
 myplotfun <- function(bw.l,tmp,gr.v,ylim,xlim=xlim,bin=bin){
   for(mygr in gr.v){
     sze <- length(get(mygr))
     # print(mygr)
     o.tmp <- toString(rtracklayer::export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed")))
     # o.tmp <- toString(export.bed(mygr,paste0(tmp,"/",mygr,"_#",sze,"peaks.bed")))
     #plotTMP <- getPlotSetArray(bw.l,o.tmp,"hg38",bin = bin,ignore_strand = F,xmin = xlim[1],xmax=xlim[2],rm0 = F, type = "af")
      plotTMP <- getPlotSetArray(bw.l,o.tmp,"dm6",bin = bin,ignore_strand = F,xmin = xlim[1],xmax=xlim[2],rm0 = F, type = "af")
file.remove(paste0(o.tmp))
     print(plotAverage(plotTMP,xlab='Relative position [bp]', ylim=ylim, ylab='number of reads',keepratio = F,error.estimates = T,colvec =c("black","firebrick2")))
     print(plotAverage(plotTMP,xlab='Relative position [bp]', ylim=ylim, ylab='number of reads',keepratio = F,error.estimates = T,colvec =c("black","firebrick2"), legend=F))
   }
 }

#####################################################################################-

seqPlotSDoutliers_scaleFact <- function(bw.l,tmp,gr.v,ylim,xlim=xlim,bin=bin,anchor=10000,sd=c(F,10),err=T,type="pf",smooth=T,spar=0.7, KR = F,gnme=NA,ignore.strand=F, scalingF = c(1,1), colvec = c("black","firebrick2"))
{
  bw.n <- NULL
  o.tmp <- NULL
  for(n in 1:length(bw.l)){
    bw.c <- bw.l[n]
    bw.n[n] <- gsub("(.*).bw","\\1",basename(bw.c))
  }
  for(mygr in gr.v){
    sze <- length(get(mygr))
    print(mygr)
    o.tmp <- c(o.tmp,toString(rtracklayer::export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed"))))
  }

  gpsa <- getPlotSetArray(bw.l,o.tmp,gnme,bin = bin,ignore_strand = ignore.strand,xmin = xlim[1],xmax=xlim[2],rm0 = F,type=type, xanchored=anchor)
  gpsa.data <- gpsa$data
  for(mygr in gr.v){
    for(my.bw in bw.n){
      sze <- length(get(mygr))
      gpsa.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["heatmap"]])
      Value_per_bin = c(gpsa.mtx)
      means <- colMeans(gpsa.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
      means[which(is.nan(means))] <- 0 # change NA in 0
      if(sd[1] %in% T){
        gpsa.scl.mtx <- gpsa.mtx %>% mutate_all(scale) # scale the data (center reduce)
        gpsa.scl.mtx[abs(gpsa.scl.mtx) > sd[2]] <- NA # Remove value X SD away (sd = 3 by default ~ 98% of the data)
        means <- colMeans(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
        stderror <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,2,function(n){
          sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
        })
        conint <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx, 2, function(n) {
          qt(0.95, sum(!is.na(n))) * sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
        })
        stderror[is.na(stderror)] <- 0
        conint[is.na(conint)] <- 0
        if(smooth){
          means = smooth.spline(1:(length(means)), means, spar=spar)$y
       
        }
        gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["stderror"]] <- stderror # change the means vector from getPlotSetArray object
        gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["conint"]] <- conint  # change the means vector from getPlotSetArray object
      }
      if(smooth){
        means = smooth.spline(1:(length(means)), means, spar=spar)$y
      }
      if(my.bw == bw.n[[scalingF[1]]]){
        gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["means"]] <- means*scalingF[2]      #16.33333*1.761243 # change the means vector from getPlotSetArray object + scale fact for comparison

      }else{
        gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["means"]] <- means # change the means vector from getPlotSetArray object
      }
    }

  }
  file.remove(paste0(o.tmp))
  par(mfrow=c(1,1))
  plotAverage(gpsa,xlab='Relative position [bp]', ylim=ylim, ylab='Signal', main = paste0("Plot profile \n SD Removed_",sd[1]," ",sd[2]), keepratio = KR,error.estimates = err, colvec = colvec, pointsize = 20, legend_ext = T)
  plotAverage(gpsa,xlab='Relative position [bp]', ylim=ylim, ylab='Signal', main = paste0("Plot profile \n SD Removed_",sd[1]," ",sd[2]), keepratio = KR,error.estimates = err, colvec = colvec, pointsize = 20, legend=F)
  # par(mfrow=c(2,2))
  # plot(density(Value_per_bin))

}



#####################################################################################-

split_GRanges_inList <- function(GR, NnamesSplit, Nsplit = NULL){
  # namesSplit is either a character vector (then  decreasinglyordered splitted with Nsplit)
  # either a list of character vectors (then a list of granges is created according to list of names )
  namesSplit = get(NnamesSplit)
  if(is.numeric(namesSplit)){
    namesSplit = names(namesSplit[order(namesSplit, decreasing=T)])
    GR = GR[namesSplit]
    GRList = split(GR, ceiling(seq_along(GR)/ceiling(length(GR)/Nsplit)))
    names(GRList) =  paste0("GR_",NnamesSplit, "_Q" ,seq(1,Nsplit,1))
  }
  if(is.character(namesSplit)){
    GR = GR[namesSplit]
    GRList = split(GR, ceiling(seq_along(GR)/ceiling(length(GR)/Nsplit)))
    names(GRList) =  paste0("GR_",NnamesSplit, "_Q" ,seq(1,Nsplit,1))
  }
  if(is.list(namesSplit)){
    GRList = list()
    GRList = unlist(lapply(namesSplit, function(subnames){GRList = c(GRList, GR[names(GR) %in% subnames])}))
    names(GRList) =  paste0("GR_",names(namesSplit))
  }
  return(GRList)
}

## Sort each list of genes in POSDOM according to a quantif
sortListGenes = function(GNList, Quantif){
## EXEMPLE ::  List_genes_DOM_pause_indice_ctrl = sortListGenes(GNList = List_genes_DOM, Quantif = pause_indice_ctrl)
  lapply(GNList, function(GN){
    QuantifGN = Quantif[names(Quantif) %in% GN]
    GN = names(QuantifGN[order(QuantifGN, decreasing=T)])
  })
}



#### GET GN LIST of a vector decile or else
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
#######################""

############################################################################################################
#     																									   #
#										BOXPLOT 														   #
#
############################################################################################################



Boxplot_wilcoxListFilter_REF = function(quantifWT, quantifKD, cond1 = "WT", cond2="KD", filterGNList, effMin = NULL,  SampleNorm = c(F, "NULL"), YLIM = NULL, bxplt_color=NULL, outlierTH = 0.01, logTrans =T, outdir, readQuantif = "readQuantif", Cond = "KD", select = "filterValue", info = "")
{
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

  # bxplt_color
  if(is.null(bxplt_color)){
    bxplt_color =  c("#aaaaaa", "#5e5e5e")
  }

#####################################

 #head(df_toPlot)
#readsCount select       condition
# 9.574139 Random Q_H2AV_GB_PW_R1
# 8.763884 Random Q_H2AV_GB_PW_R1
#10.442389 Random Q_H2AV_GB_PW_R1
# 8.651634 Random Q_H2AV_GB_PW_R1
# 9.517005 Random Q_H2AV_GB_PW_R1
# 9.160054 Random Q_H2AV_GB_PW_R1


# vecteur de position des p val dans le plot 
y_position_val= max(quantifWT)+1
y_position_val= trunc(y_position_val)
y_position_vect=rep(y_position_val,length.out= length(filterGNList))

### Paired wilcox_test
stat.test_1 <- df_toPlot %>%
group_by(select) %>%
wilcox_test(readsCount ~ condition,paired = TRUE) %>% # paire test 
adjust_pvalue(method = "bonferroni") %>%
add_significance("p")
stat.test_1
##### NO_ paired wilcox_test
stat.test_2 <- df_toPlot %>%
group_by(select) %>%
wilcox_test(readsCount ~ condition,paired = FALSE) %>% # NO paire test 
adjust_pvalue(method = "bonferroni") %>%
add_significance("p")
stat.test_2

# Créer  boxplot 1 
bxp_1 <- ggboxplot(df_toPlot, x = "select", y = "readsCount", color = "black", fill="condition", palette = bxplt_color,title="wilcox paire test")
# Ajoutez des p-values sur les graphiques en box plot
stat.test_1 <- stat.test_1 %>% add_xy_position(x = "select", dodge = 0.8)
# ajouter les crochets + des etoiles 
bx1=bxp_1 + stat_pvalue_manual(stat.test_1,  label = "{p}", tip.length = 0,remove.bracket = FALSE, y.position=y_position_vect)
#bx1=bxp_1 + stat_pvalue_manual(stat.test_1,  label = "{p}{p.signif}",step.increase=0.001 , y.position=y_position_vect)

# Créer box plot 2 
bxp_2 <- ggboxplot(df_toPlot, x = "select", y = "readsCount", color = "black", fill="condition",  palette = bxplt_color,
title= "wilcox NO paire test")
# Ajoutez des p-values sur les graphiques en box plot
stat.test_2 <- stat.test_2 %>% add_xy_position(x = "select", dodge = 0.8)
# ajouter les crochets
bx2=bxp_2 + stat_pvalue_manual(stat.test_2,  label = "{p}", tip.length = 0,remove.bracket = FALSE, y.position=y_position_vect)

print(bx1)
textplot(stat.test_1 )
print(bx2)
textplot(stat.test_2)
textplot(lapply(filterGNList, length))
textplot(lapply(filterGNList_effMin, length))

}
