
## Get Arguments
args=(commandArgs(TRUE))
date<-Sys.Date()
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)< 12){
	stop("Not enough argument in command line ! See ReadMe.txt for required inputs.")
}else{
	for(i in 1:length(args)){
		eval(parse(text=args[[i]]))
	}
	data = data # data type called 'chipseq' or 'rnaseq' literally
	file_path = file_path # path to input file
	file_name = file_name # input file name
	window_start = w_start # FOR chipseq uniquely : start of window (distance from TSS) ex: -2000
	window_end = w_end # FOR chipseq uniquely : end of window (distance from TSS) ex: 2000
	design = design # all conditions according to reads coutn matrix ex: for ncol=6 matrix, design="control,control,factor1_KD,factor1_KD,factor2_KD,factor2_KD"
	cond1 = cond1 # first condition to compare in DE analysis, must be strictelly same name as in design
	cond2 = cond2 # second condition to compare to the first one in DE analysis, must be strictelly same name as in design
	filter = filt # filt='TRUE' if we want a HTSFilter to be apply onc ount matrix
	th_fc = FC # threshold foldchange pour volcano plot
	out_path = out_path # path to output files
	out_name = out_name # Out put prefix
}

#~ 	data = 'chipseq' # data type called 'chipseq' or 'rnaseq' literally
#~ 	file_path = "/home/ddepierre/work/David_M2_K9K27DE/DHMG/CountMatProf/profiles_H3k9me3" # path to input file
#~ 	file_name = "prof_k9Norm_list.rds" # input file name
#~ 	window_start = "-2000" # FOR chipseq uniquely : start of window (distance from TSS) ex: -2000
#~ 	window_end = "0" # FOR chipseq uniquely : end of window (distance from TSS) ex: 2000
#~ 	design = "control,control,HypB_KD,HypB_KD,Mes4_KD,Mes4_KD" # all conditions according to reads coutn matrix ex: for ncol=6 matrix, design="control,control,factor1_KD,factor1_KD,factor2_KD,factor2_KD"
#~ 	cond1 = "control" # first condition to compare in DE analysis, must be strictelly same name as in design
#~ 	cond2 = "Mes4_KD" # second condition to compare to the first one in DE analysis, must be strictelly same name as in design
#~ 	filter = "TRUE" # filt='TRUE' if we want a HTSFilter to be apply onc ount matrix
#~ 	th_fc = "1" # threshold foldchange pour volcano plot
#~ 	out_path = "/home/ddepierre/work/David_M2_K9K27DE/DHMG/DHMG_analysis/DHMG_limma_FROMm2kb_TOp0kb" # path to output files
#~ 	out_name = "FROMm2kb_TOp0kb" # Out put prefix


design = unlist(strsplit(design, ','))
th_fc = as.numeric(th_fc)

print(paste0("data type : ", data))
print(paste0("file_path : ", file_path))
print(paste0("file_name : ", file_name))
print(paste0("window_start : ", window_start))
print(paste0("window_end : ", window_end))
print(paste0("design : ", design))
print(paste0("cond1 : ", cond1))
print(paste0("cond2 : ", cond2))
th_fc
print(paste0("HTSFilter : ", filter))
print(paste0("out_path : ", out_path))
print(paste0("out_name : ", out_name))


library("limma")
library("edgeR")
library("cluster")
library(mixOmics)
library(gplots)
library(ggplot2)
library(HTSFilter)


if(data == 'chipseq'){
	# Get profiles list file
	prof = readRDS(file.path(file_path, file_name))
	# Get profiles for the given window
	zoomIn = function(profil, start, end){
		s = as.numeric(start) + round(dim(profil)[2]/2, 0) +1
		e = as.numeric(end) + round(dim(profil)[2]/2 ,0) +1
		profilZoomed = profil[,c( s : e )]
	}
	prof2 = lapply(prof, FUN=zoomIn, window_start, window_end)
	# Get sum of count by gene for the given window and put it in a matrix
	byGeneRowSum = lapply(prof2, rowSums)
	matByGene = do.call(cbind, byGeneRowSum)
	head(matByGene)
}else if (data == 'rnaseq'){
	matByGene = readRDS(file.path(file_path, file_name))
}else {
	stop("Argument data must be 'chipseq' or 'rnaseq' ! See ReadMe.txt for required inputs.")
}


## Data filtering HTSFilter
if(filter == 'TRUE'){
	condition = as.factor(design) 
	condition
	matByGene = HTSFilter(matByGene, condition, plot=FALSE)$filteredData
}else{
	matByGene = matByGene[rowSums(matByGene)>0,]
}

# get the splitted matrix
subMatByGene = matByGene[,design %in% c(cond1,cond2)]

## Create a data frame of conditions
matExpeDesign = matrix(ncol=1, nrow =4)
row.names(matExpeDesign) = colnames(subMatByGene)
colnames(matExpeDesign) = "condition"
matExpeDesign[c(1,2),"condition"] = cond1
matExpeDesign[c(3,4),"condition"] = cond2
matExpeDesign = as.data.frame(matExpeDesign)
matExpeDesign

dge1 = DGEList(subMatByGene)
dge1 = calcNormFactors(dge1) # default method  = TMM,  

# create a matrix condition used as design for the experimetn analysis
Design1 = model.matrix(~factor(matExpeDesign$condition)) 
Design1


Voom1 = voom(dge1, design=Design1, plot = F)
VoomedMat = Voom1$E
VoomedMat = cbind(rowMeans(VoomedMat[,c(1,2)]),
				rowMeans(VoomedMat[,c(3,4)]))
colnames(VoomedMat) = c(paste0("AveExpr_",cond1), paste0("AveExpr_",cond2))

fit1 = lmFit(Voom1, Design1, method = "ls") 

ebayes1 = eBayes(fit1)


outputLimma1 = topTable(ebayes1, coef=ncol(Design1), number = dim(ebayes1)[1])
VoomedMat = VoomedMat[rownames(outputLimma1),]
outputLimma1 = cbind(outputLimma1, VoomedMat)
dim(outputLimma1)
head(outputLimma1)

write.table(as.data.frame(outputLimma1),file=paste0(out_path, "/", "output_matrix_Limma_",data,"_",out_name,"_",cond1 ,"_",cond2 ,".txt"))
saveRDS(outputLimma1, file=paste0(out_path, "/", "output_matrix_Limma_",data,"_",out_name,"_",cond1 ,"_",cond2 ,".rds"))

##############################################################################
############################ PLOT  ##########################################
####################################################################################

### VOLCANO PLOT + QQplot


ggd.qqplot = function(pvector, main=NULL, ...) {
    o = -log10(sort(pvector,decreasing=F))
    e = -log10( 1:length(o)/length(o) )
    plot(e,o,pch=19,cex=1, main=main, ...,
        xlab=expression(Expected~~-log[10](italic(p))),
        ylab=expression(Observed~~-log[10](italic(p))),
        xlim=c(0,max(e)), ylim=c(0,max(o[is.finite(o)])))
    lines(e,e,col="red")
}


up = dim(outputLimma1[outputLimma1['logFC'] >log2(th_fc) & outputLimma1[,'P.Value'] < 0.05 ,])[1]
down = dim(outputLimma1[outputLimma1['logFC'] <(-log2(th_fc)) & outputLimma1[,'P.Value'] < 0.05 ,])[1]

up2 = dim(outputLimma1[outputLimma1['logFC'] >log2(th_fc) & outputLimma1[,'adj.P.Val'] < 0.05 ,])[1]
down2 = dim(outputLimma1[outputLimma1['logFC'] <(-log2(th_fc)) & outputLimma1[,'adj.P.Val'] < 0.05 ,])[1]

pdf(paste0(out_path,"/","Plot_",cond1 ,"-vs-",cond2,"_" ,out_name, "_FC>",th_fc ,".pdf"))
# Volcano plot
attach(outputLimma1)
ggplot(outputLimma1) +
	geom_point(aes(x = logFC, y = -log(P.Value), colour =logFC < -log2(th_fc) & P.Value < 0.05 | logFC >log2(th_fc) & P.Value < 0.05), shape=1) +
	scale_colour_manual(name = paste0('raw P.Val < 0.05 & |FoldChange| > ',th_fc), values = setNames(c('darkblue','black'),c(T, F))) +
	xlab(paste("log(FoldChange)",cond1 ,"vs",cond2, sep=" ")) + ylab('-log(P.Val)') +
	scale_shape(solid=FALSE) +xlim(-9,9) + ylim(0,31)+
	annotate("text", x = 6, y = 6, label = up) +
	annotate("text", x = -5, y = 6, label = down)
	ggtitle(paste0(cond1 ,"-vs-",cond2,"_" ,out_name, "_FC>",th_fc, "_rawPval"))

attach(outputLimma1)
ggplot(outputLimma1) +
	geom_point(aes(x = logFC, y = -log(adj.P.Val), colour =logFC < -log2(th_fc) & adj.P.Val < 0.05 | logFC >log2(th_fc) & adj.P.Val < 0.05), shape=1) +
	scale_colour_manual(name = paste0('adj.P.Val < 0.05 & |FoldChange| > ',th_fc), values = setNames(c('darkblue','black'),c(T, F))) +
	xlab(paste("log(FoldChange)",cond1 ,"vs",cond2, sep=" ")) + ylab('-log(P.Val)') +
	scale_shape(solid=FALSE) +xlim(-9,9) + ylim(0,31)+
	annotate("text", x = 6, y = 6, label = up2) +
	annotate("text", x = -5, y = 6, label = down2) + 
	ggtitle(paste0(cond1 ,"-vs-",cond2,"_" ,out_name, "_FC>",th_fc, "_adjPval"))

### QQ PLOT

	ggd.qqplot(outputLimma1[,'P.Value'], main=paste0("QQPlot_",cond1 ,"-vs-",cond2,"_" ,out_name, "_rawPval"))
	ggd.qqplot(outputLimma1[,'adj.P.Val'], main=paste0("QQPlot_",cond1 ,"-vs-",cond2,"_" ,out_name, "_adjPval"))
dev.off()


stop("End of processing")




