# Single Cell RNA seq project
# Cuvier's Team
# Schaak - Heurteau - Depierre* 
# 2017


# From genes dm6 GRanges object, get genes liste that are isolated i.e no overlapping TSS region +/-1000 ? +1000
# The goal is to have a genes list for which we have only the TSS of 1 unique gene which have an influence

# function taking Granges, start and stop position around TSS of where we want to find overlapping

#~ genes_dm6 = readRDS("/home/belhocine/Bureau/work_ddepierre/David_M2_K9K27DE/DHMG/r6.13/genes_dm6_Granges.rds")
#~ genes_dm6 = readRDS("/home/ddepierre/work/David_M2_K9K27DE/DHMG/r6.13/genes_dm6_Granges.rds")

#~ gRange = genes_dm6
#~ start = -500
#~ end = 0

require(GenomicRanges)


getIsolatedGeneGRanges = function(gRange, start = 0, end = 0, omit.str = T){
	# USAGE : getSingleGeneGRanges(gRange, start=-500, end = 0, omit.str = T)
	# Output will be GRanges object of genes that have no overlapping with other gene body on their -500/0 around TSS
	# gRange = genes GRanges
	#
	# omit.str for ignore.strand param of findOverlaps, if TRUE, overlaps will be not strand specific, 
	#												, if FALSE, overlaps found only between same strand genes
	# Get Granges of TSS -start/+end
	gRange_start_end = gRange
	if(end >=0){
		end(gRange_start_end[strand(gRange_start_end)=='+'])=start(gRange[strand(gRange)=='+'])+end
		start(gRange_start_end[strand(gRange_start_end)=='+'])=start(gRange[strand(gRange)=='+'])+start
		start(gRange_start_end[strand(gRange_start_end)=='-'])=end(gRange[strand(gRange)=='-'])-end
		end(gRange_start_end[strand(gRange_start_end)=='-'])=end(gRange[strand(gRange)=='-'])-start
		# get gene with no overlapping 
		gene_Ovlp = findOverlaps(gRange_start_end, gRange, ignore.strand=omit.str)
		gene_Ovlp = rle(gene_Ovlp@queryHits)
		gRange_NoOvlp = gRange[gene_Ovlp$values[which(gene_Ovlp$lengths == 1)]]
		return(gRange_NoOvlp)
	}
	if(end < 0){
		start(gRange_start_end[strand(gRange_start_end)=='+'])=start(gRange[strand(gRange)=='+'])+start
		end(gRange_start_end[strand(gRange_start_end)=='+'])=start(gRange[strand(gRange)=='+'])+end
		end(gRange_start_end[strand(gRange_start_end)=='-'])=end(gRange[strand(gRange)=='-'])-start
		start(gRange_start_end[strand(gRange_start_end)=='-'])=end(gRange[strand(gRange)=='-'])-end
		# get gene with no overlapping 
		gene_Ovlp = findOverlaps(gRange_start_end, gRange, ignore.strand=omit.str)
		gRange_NoOvlp = gRange[-gene_Ovlp@queryHits]
		return(gRange_NoOvlp)
	}
}


#~ length(getIsolatedGeneGRanges(genes_dm6, start = -1000, end =0, omit.str=T))
#~ [1] 4455
