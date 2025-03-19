# Single Cell RNA seq project
# Cuvier's Team
# Schaak - Heurteau - Depierre
# 2017
library(GenomicRanges)
library(Rsamtools)
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library(gsubfn)



bam2BigWig = function(bam_path, out_path, bam_pattern){
	######################################################
	# Single Cell RNA seq project
	# Cuvier's Team
	# Schaak - Heurteau - Depierre
	# 2017
	#
	#     USAGE :
	#> bam2BigWig(path/to/bamfiles, path/to/output, .bam)
	#
	#   	Function which takes bam files path and pattern file name
	#   	And computes BigWig file using GRanges in order to visualization
	#
	#   	out_path must already exists
	#   	bam_pattern correpsond to consensus end of your bamfiles
	#   	ex : for file "myBam_filtered.bam", pattern could be "_filtered.bam"



	######################################################
	pat=paste0("(.*)(", as.symbol(bam_pattern), ")")
	bamfile = dir(bam_path, pattern=pat)
	bamfile = paste0(bam_path, "/", bamfile)
	name = strapplyc(bamfile, pat, simplify = TRUE)[1,]
	bam = lapply(bamfile, readGAlignments)
	bam_gr = lapply(bam, granges)
	names(bam_gr) = name
  # IGV is in UCSC style :
	bam_gr_UCSC = lapply(bam_gr, function(GR){seqlevelsStyle(GR)="UCSC"})
		myExport = function(granges, NameGR){
			export.bw(granges, format="BigWig", con=paste0(out_path,"/", NameGR ,".BigWig"))
		}
	mapply(myExport, bam_gr_UCSC, name)

  }
