# Single Cell RNA seq project
# Cuvier's Team
# Schaak - Heurteau - Depierre*
# 2017
# adaptation from Gael Micas script Dm541_pm2Kb.R


profcomp_tx=function(covr,gr)
#Compute profiles from a genome-wide coverage and a set of windows
#covr is a genome-wide coverage (typically obtained with the coverage function)
#gr is a GRanges object containing the windows over which to compute the profiles
#Note that names of covr and seqnames of gr must match
{
covV=RleViewsList()
for (i in names(covr))
	{covV[[i]]=Views(covr[[i]],ranges(gr[seqnames(gr)==i]))
	names(covV[[i]])=names(gr[seqnames(gr)==i])}
#Compute coverage over the windows
prof_chr=list()
for (i in names(covV[sapply(covV,length)!=0]))
	{
	prof_chr[[i]]=sapply(covV[[i]],function(x){Rle(as.vector(x))})
	prof_chr[[i]][as.vector(strand(gr[names(prof_chr[[i]])])=="-")]=lapply(prof_chr[[i]][as.vector(strand(gr[names(prof_chr[[i]])])=="-")],rev) #reverse the profiles for minus strand
	}
#Organize in a list
return(prof_chr)
}



covr = cov_H3K27me3$H3K27me3_ctrl
gr = dom_heterochr


profcomp_dom=function(covr,gr)
#Compute profiles from a genome-wide coverage and a set of windows
#covr is a genome-wide coverage (typically obtained with the coverage function)
#gr is a GRanges object containing the windows over which to compute the profiles
#Note that names of covr and seqnames of gr must match
{
covV=RleViewsList()
for (i in names(covr))
	{covV[[i]]=Views(covr[[i]],ranges(gr[seqnames(gr)==i]))
	names(covV[[i]])=names(gr[seqnames(gr)==i])}
#Compute coverage over the windows
prof_chr=list()
for (i in names(covV[sapply(covV,length)!=0]))
	{
	prof_chr[[i]]=sapply(covV[[i]],function(x){Rle(as.vector(x))})
	prof_chr[[i]][as.vector(strand(gr[names(prof_chr[[i]])])=="-")]=lapply(prof_chr[[i]][as.vector(strand(gr[names(prof_chr[[i]])])=="-")],rev) #reverse the profiles for minus strand
	}
	prof_dom = list()
	for(i in 1:length(names(prof_chr))){
		prof_dom = c(prof_dom, lapply(prof_chr[[i]], as.numeric))
	}
return(prof_dom)
}


profcomp=function(covr,gr)
#Compute profile matrix from a genome-wide coverage and a set of windows
#covr is a genome-wide coverage (typically obtained with the coverage function -use resize(granges(bam_aln),FragmentSize) for ChIP-seq data-)
#gr is a GRanges object containing the windows (all the same size) over which to compute the profiles (e.g. promoter regions)
#Note that names of covr and seqnames of gr must match
{
covV=RleViewsList()
for (i in names(covr))
	{covV[[i]]=Views(covr[[i]],ranges(gr[seqnames(gr)==i]))
	names(covV[[i]])=names(gr[seqnames(gr)==i])}
#Compute coverage over the windows
prof=matrix(rep(NA,max(width(gr))),nrow=1)
for (i in names(covV[sapply(covV,length)!=0]))
	{prof=rbind(prof,as.matrix(covV[[i]]))}
prof=prof[-1,]
#Organize in a matrix
prof[as.vector(strand(gr[rownames(prof)])=="-"),]=prof[as.vector(strand(gr[rownames(prof)])=="-"),ncol(prof):1]
#Reverse the data for the genes on the minus strand
return(prof)
}
