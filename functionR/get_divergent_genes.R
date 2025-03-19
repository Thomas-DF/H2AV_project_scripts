# Single Cell RNA seq project
# Cuvier's Team
# Schaak - Heurteau - Depierre* 
# 2017


# From genes dm6 GRanges object, get genes list that are isolated i.e no overlapping TSS region +/-1000 ? +1000
# The goal is to have a genes list for which we have only the TSS of 1 unique gene which have an influence

# function taking Granges, start and stop position around TSS of where we want to find overlapping

#~ genes_dm6 = readRDS("/home/belhocine/Bureau/work_ddepierre/David_M2_K9K27DE/DHMG/r6.13/genes_dm6_Granges.rds")
#~ genes_dm6 = readRDS("/home/ddepierre/work/David_M2_K9K27DE/DHMG/r6.13/genes_dm6_Granges.rds")

#~ gRange = genes_dm6
#~ start = -500
#~ end = 0

require(GenomicRanges)


get_divergent_genes = function(gRange, start = 0, end = 0, omit.str = T){
	# USAGE : get_divergent_genes(gRange, start=-500, end = 0, omit.str = T)
	# Output will be GRanges object of divergent genes i.e. that have a TSS from a gene from an other strand on their -500/0 around TSS
	# gRange = genes GRanges
	#
	# omit.str for ignore.strand param of findOverlaps, if TRUE, overlaps will be not strand specific, 
	#												, if FALSE, overlaps found only between same strand genes
	# Get Granges of TSS -start/+end
	gRange_start_end = gRange
	gRange_TSS = gRange
	# Get  granges of TSS 0/+1
		end(gRange_TSS[strand(gRange_TSS)=='+'])=start(gRange[strand(gRange)=='+'])+1
		start(gRange_TSS[strand(gRange_TSS)=='+'])=start(gRange[strand(gRange)=='+'])+0
		start(gRange_TSS[strand(gRange_TSS)=='-'])=end(gRange[strand(gRange)=='-'])-1
		end(gRange_TSS[strand(gRange_TSS)=='-'])=end(gRange[strand(gRange)=='-'])-0
	if(end >=0){
		end(gRange_start_end[strand(gRange_start_end)=='+'])=start(gRange[strand(gRange)=='+'])+end
		start(gRange_start_end[strand(gRange_start_end)=='+'])=start(gRange[strand(gRange)=='+'])+start
		start(gRange_start_end[strand(gRange_start_end)=='-'])=end(gRange[strand(gRange)=='-'])-end
		end(gRange_start_end[strand(gRange_start_end)=='-'])=end(gRange[strand(gRange)=='-'])-start
	}
	if(end < 0){
		start(gRange_start_end[strand(gRange_start_end)=='+'])=start(gRange[strand(gRange)=='+'])+start
		end(gRange_start_end[strand(gRange_start_end)=='+'])=start(gRange[strand(gRange)=='+'])+end
		end(gRange_start_end[strand(gRange_start_end)=='-'])=end(gRange[strand(gRange)=='-'])-start
		start(gRange_start_end[strand(gRange_start_end)=='-'])=end(gRange[strand(gRange)=='-'])-end
	}
	gene_Ovlp = findOverlaps(gRange_start_end, gRange_TSS, ignore.strand=omit.str)
	matOvlp = cbind(gene_Ovlp@queryHits, gene_Ovlp@subjectHits)
	cplOvlp = apply(matOvlp, 1, function(coupleOvlp){coupleOvlp[1] != coupleOvlp[2]})
	Ovpldiff = matOvlp[cplOvlp,]
	matStrand = cbind(as.vector(strand(gRange[Ovpldiff[,1]])), as.vector(strand(gRange[Ovpldiff[,2]])))
	cplStrand = apply(matStrand, 1, function(coupleOvlp){coupleOvlp[1] != coupleOvlp[2]})
	OvlpInvStr = Ovpldiff[cplStrand,]
	return(gRange[unique(OvlpInvStr[,1])])
}

























