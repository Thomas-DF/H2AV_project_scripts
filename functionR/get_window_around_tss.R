# Single Cell RNA seq / k9me3.k27me3 project
# Cuvier's Team
# Depierre
# 2017

require(GenomicRanges)

getWindowAroundTSS = function(gRange, start = 0, end = 0){
	# USAGE : getWindowAroundTSS(gRange, start=-500, end = 0)
	# Output will be GRanges object of genes on their -500/0 around TSS
	# gRange = genes GRanges
	#
	# Get Granges of TSS -start/+end
	gRange_start_end = gRange
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
return(gRange_start_end)
}


