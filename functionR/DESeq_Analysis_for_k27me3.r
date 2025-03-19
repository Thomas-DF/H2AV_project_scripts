# Differential histone modification genes
# As we don't have biological replicates for this experimetn we use DESeq2 to apply an rlog Transformation
# and extract a rlogFC by simply subtract the control LogExpression on knockDown LogExpression
# !! No replicates = naive apporoach → no statistic evaluation, no pvalue

## Get Arguments
args=(commandArgs(TRUE))
date<-Sys.Date()
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)< 11){
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
	out_path = out_path # path to output files
	out_name = out_name # Out put prefix
}

# data = 'chipseq'
# file_path = "/home/ddepierre/work/David_M2_K9K27DE/MNase_seq/CountMatProf/NORM_byMNASEseq"
# file_name = "prof_k27_mnase_norm.rds"
# window_start = "-2000"
# window_end = "0"
# design = "control,Mes4_KD,HypB_KD"
# cond1 = "control"
# cond2 = "Mes4_KD"
# filter = "FALSE"
# out_path = "/home/ddepierre/work/David_M2_K9K27DE/DHMG_normONmnase/DHMG_analysis_k27me3/DHMG_DESeq2_FROMm2kb_TOp0kb"
# out_name = "FROMm2kb_TOp0kb"


design = unlist(strsplit(design, ','))

print(paste0("data type : ", data))
print(paste0("file_path : ", file_path))
print(paste0("file_name : ", file_name))
print(paste0("window_start : ", window_start))
print(paste0("window_end : ", window_end))
print(paste0("design : ", design))
print(paste0("cond1 : ", cond1))
print(paste0("cond2 : ", cond2))

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
library('DESeq2')


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


subMatByGene = matByGene[,design %in% c(cond1,cond2)]
head(subMatByGene)


colData = as.data.frame(cbind("condition"=c(cond1,cond2)))
rownames(colData) = colnames(subMatByGene)
subMatByGene = round(subMatByGene)

dds <- DESeqDataSetFromMatrix(countData = subMatByGene, colData = colData, design = ~ condition)

rld <- rlogTransformation( dds )

res <- data.frame(assay(rld), avgLogExpr = ( assay(rld)[,2] + assay(rld)[,1] ) / 2, rLogFC = assay(rld)[,2] - assay(rld)[,1] )

head(res[order(res$rLogFC, decreasing=T),])

outputDESeq = res[order(res$rLogFC, decreasing=T),]

pdf(paste0(out_path,"/","PlotDensity_rLogFC_",cond1 ,"-vs-",cond2,"_" ,out_name,".pdf"))

plot(density(outputDESeq$rLogFC))
dev.off()

#~ ## III. Descriptive analysis
#~ ### 1 - ACP
#~ ```{r ACP  + heatmap}
#~ pdf('PCA_Heatmap_for_DEanalysis_RNAseq_mes4-hypb-KD.pdf')
#~ plotIndiv(pca(t(dge01.cpm)), col = c(1,1,2,2,3,3))
#~ ### 2 - Heatmap
#~ heatmap.2(scale(t(dge01.cpm), scale =T), col =  greenred(10), scale = "column", trace = 'none', main='heatmap des data d expression en mes4 et hypb deplete')
#~ dev.off()
#~ ```


write.table(as.data.frame(outputDESeq),file=paste0(out_path, "/", "output_matrix_DESeq2_",data,"_",out_name,"_",cond1 ,"_",cond2 ,".txt"))
saveRDS(outputDESeq, file=paste0(out_path, "/", "output_matrix_DESeq2_",data,"_",out_name,"_",cond1 ,"_",cond2 ,".rds"))


stop("End of processing")
