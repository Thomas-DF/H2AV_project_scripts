
# Single Cell RNA seq project
# Cuvier''s Team
# Schaak - Heurteau - Depierre*
# 2017

#~ Heat map plot function
#~ To plot count matrix around TSS



library(gplots)
library(pheatmap)
'%ni%' = Negate('%in%')
require(GenomicRanges)
require(BiocGenerics)
require(parallel)
source(paste0(workdir,"David_M2_K9K27DE/functionR/get_isolated_genes.R"))
source(paste0(workdir,"David_M2_K9K27DE/functionR/bin_profile_mat.R"))


#~ cov1 = h3k27me3_ctrl
#~ cov2 = h3k27me3_hypb
#~ sortby = rownames(DErna_hypb[order(DErna_hypb[,"logFC"], decreasing=F),])
#~ start = -1000
#~ end = 1000
#~ outdir = "/home/ddepierre/work/David_M2_K9K27DE/heatmap_prof/"
#~ figure.name = "testheatmap2"


plot_heatmap_coverage = function(cov1, cov2 = NULL, sortby = NULL, start = -1000, end = 1000, outdir, figure.name, Main ="heatmap", Tlog = TRUE){
	if(is.null(sortby)){
		if(is.null(cov2)){
			# get matrix in choosen window between start and stop around TSS
			cov1 = cov1[,(round(dim(cov1)[2]/2,0)+start+1):(round(dim(cov1)[2]/2,0)+end)]
			if(Tlog ==TRUE){
				cov1 = log2(cov1+1)
			}
			# bin matrix
			cov1 = t(apply(cov1,1,binMean,xTend=ncol(cov1),binSize=10))

			# plot matrix
			bk = unique(c(seq(min(range(cov1)),3, length=100),seq(3,max(range(cov1)),length=100)))
			bicols<- colorRampPalette(c("white","red4"))(length(bk)-1)
			pdf(paste0(outdir, figure.name, ".pdf"), width=10, height = 30)     # width and height are in pixel
			par(cex.main=2)
			heatmap.2(cov1, col=bicols, breaks = bk, Rowv= T , Colv=F, dendrogram="row", useRaster = TRUE, symkey=FALSE, symm=F, symbreaks=T,
					scale="none", trace="none", labRow=NA, labCol=NA, main=Main)
			dev.off()
		}else{
			# get matrix in choosen window between start and stop around TSS
			cov = cov2 - cov1
			cov = cov1[,(round(dim(cov1)[2]/2,0)+start+1):(round(dim(cov1)[2]/2,0)+end)]
			if(Tlog ==TRUE){
				cov = log2(cov+1)
			}
			# bin matrix
			cov = t(apply(cov,1,binMean,xTend=ncol(cov),binSize=10))
			# plot heatmap
			bk = unique(c(seq(min(range(cov)),-0.01,length=100),seq(-0.01,0.01,length=100),seq(0.01,max(range(cov)),length=100)))
			tricols =  colorRampPalette(c("red4","white","green4"))(length(bk)-1)
			pdf(paste0(outdir, figure.name, ".pdf"), width=10, height = 30)     # width and height are in pixel
			par(cex.main=2)
			heatmap.2(cov, col=tricols, breaks = bk, Rowv= TRUE , Colv=FALSE, dendrogram="row", useRaster = TRUE, symkey=FALSE, symm=F, symbreaks=T,
					scale="none", trace="none", labRow=NA, labCol=NA, main=Main)
			dev.off()

		}
	}else{
		if(is.null(cov2)){
			comRownames = Reduce(intersect, list(rownames(cov1), sortby))
			cov1 = cov1[rownames(cov1) %in% comRownames,]
			sortby = sortby[sortby %in% comRownames]
			cov1 = cov1[sortby,]
			cov1 = cov1[,(round(dim(cov1)[2]/2,0)+start+1):(round(dim(cov1)[2]/2,0)+end)]
			if(Tlog ==TRUE){
				cov1 = log2(cov1+1)
			}
			# bin matrix

			cov1 = t(apply(cov1,1,binMedian,xTend=ncol(cov1),binSize=10))
			# cov1 = t(apply(cov1,1,binMean,xTend=ncol(cov1),binSize=10))

			# plot matrix

			# bk = unique(c(seq(0,8, length=1),seq(8,15,length=10),seq(15,20,length=1),seq(20,200,length=1)))
			bk = unique(c(seq(0,5, length=1),seq(5,8,length=100),seq(8,15,length=10),seq(15,20,length=1)))
			# bicols<- c("#FFFFFF", colorRampPalette(c("grey","black"))(length(bk)-3), "#000000")
			bicols<- c("#FFFFFF", rep("#000000",length(bk)-3), "#000000")
			pdf(paste0(outdir, figure.name, ".pdf"), width=10, height = 30)     # width and height are in pixel
			par(cex.main=2)
			heatmap.2(cov1, col=bicols, breaks = bk, Rowv= F , Colv=F, dendrogram="none", useRaster = TRUE, symkey=FALSE,
					symm=F, symbreaks=T, scale="none", trace="none", labRow=NA, labCol=NA, main=Main)
			# heatmap.2(cov1, Rowv= F , Colv=F, dendrogram="none", useRaster = TRUE, symkey=FALSE,
			# 		symm=F, symbreaks=T, scale="none", trace="none", labRow=NA, labCol=NA, main=Main)
			# heatmap.2(cov1, Rowv= F , Colv=F, dendrogram="none", useRaster = TRUE, symkey=FALSE,
			# 		symm=F, symbreaks=T, scale="row", trace="none", labRow=NA, labCol=NA, main=Main)
			# heatmap.2(cov1, Rowv= F , Colv=F, dendrogram="none", useRaster = TRUE, symkey=FALSE,
			# 		symm=F, symbreaks=T, scale="column", trace="none", labRow=NA, labCol=NA, main=Main)
			dev.off()
		}
		else{

			# get matrix in choosen window between start and stop around TSS
#~ 			cov1 = log2(cov1+1)
#~ 			cov2 = log2(cov2+1)
			comRownames = Reduce(intersect, list(rownames(cov1), rownames(cov2), sortby))
			cov1 = cov1[rownames(cov1) %in% comRownames,]
			cov2 = cov2[rownames(cov2) %in% comRownames,]
			sortby = sortby[sortby %in% comRownames]
			cov1 = cov1[,(round(dim(cov1)[2]/2,0)+start+1):(round(dim(cov1)[2]/2,0)+end)]
			cov2 = cov2[,(round(dim(cov2)[2]/2,0)+start+1):(round(dim(cov2)[2]/2,0)+end)]
			cov1 = cov1[rownames(cov1) %in% sortby,]
			cov2 = cov2[rownames(cov2) %in% sortby,]
			cov1 = cov1[sortby,]
			cov2 = cov2[sortby,]
			# bin matrix
			cov = cov2-cov1
			cov = t(apply(cov,1,binMean,xTend=ncol(cov),binSize=10))
			if(Tlog ==TRUE){
				cov = log2(cov+1)
			}
			#~ 			cov = t(apply(cov, 1, function(row){row = row/sqrt(mean(row))}))
			# plot heatmap
#~ 			bk = unique(c(seq(0,1, length=100),seq(1,3,length=100)))
#~ 			cols<- colorRampPalette(c("white","red4"))(length(bk)-1)
			bk = unique(c(seq(-5,-1,length=100),seq(-1,1,length=100),seq(1,5,length=100)))
			cols =  colorRampPalette(c("white","white","black"))(length(bk)-1)
			pdf(paste0(outdir, figure.name, ".pdf"), width=10, height = 30)     # width and height are in pixel
			par(cex.main=2)
			heatmap.2(cov, col=cols, breaks = bk, Rowv= FALSE , Colv=FALSE, dendrogram="none", useRaster = TRUE, symkey=FALSE, symm=F, symbreaks=T,
					scale="none", trace="none", labRow=NA, labCol=NA, main=Main)
			dev.off()
		}
	}
}

## END
