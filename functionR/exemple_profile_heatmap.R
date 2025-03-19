# Author: Heurteau Alexandre
# Date: Tue Nov  6 15:56:03 2018
#####################################################################################--
#          DESCRIPTION  ----
#####################################################################################--

# THIS SCRIPT AIMS TO


#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

source("/media/alexandre/Data/Recherche/LBME/Database/Scripts_DB/R_script_DB/lrc_lib.R")
library(rtracklayer) # Pour genomicRanges
library(GenomicRanges) # Idem
library("seqplots") # Pour les plots
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library("TxDb.Dmelanogaster.UCSC.dm3.ensGene") # pour les bases de données droso
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")  # pour les bases de données droso

#####################################################################################-
#          INSIDE FUNCTIONS  ----
#####################################################################################-
seqPlotSDoutliers <- function(bw.l,tmp,gr.v,ylim,xlim=xlim,bin=bin,sd=3,err=F,type="mf",smooth=F,spar=0.35,gnme="hg38",ignore.strand=F){
  bw.n <- NULL
  o.tmp <- NULL
  for(n in 1:length(bw.l)){
    bw.c <- bw.l[n]
    bw.n[n] <- gsub("(.*).bw","\\1",basename(bw.c))
  }
  for(mygr in gr.v){
    sze <- length(get(mygr))
    print(mygr)
    o.tmp <- c(o.tmp,toString(export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed"))))
  }
  gpsa <- getPlotSetArray(bw.l,o.tmp,gnme,bin = bin,ignore_strand = ignore.strand,xmin = xlim[1],xmax=xlim[2],rm0 = F,type=type)
  gpsa.data <- gpsa$data
  for(mygr in gr.v){
    for(my.bw in bw.n){
      sze <- length(get(mygr))
      gpsa.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["heatmap"]])
      gpsa.scl.mtx <- gpsa.mtx %>% mutate_all(scale) # scale the data (center reduce)
      gpsa.scl.mtx[abs(gpsa.scl.mtx) > sd] <- NA # Remove value X SD away (sd = 3 by default ~ 98% of the data)
      means <- colMeans(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
      if(smooth){
        means = smooth.spline(1:(length(means)), means, spar=spar)$y
      }
      stderror <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,2,function(n){
        sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
      })
      conint <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx, 2, function(n) {
        qt(0.95, sum(!is.na(n))) * sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
      })
      stderror[is.na(stderror)] <- 0
      conint[is.na(conint)] <- 0
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["means"]] <- means # change the means vector from getPlotSetArray object
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["stderror"]] <- stderror # change the means vector from getPlotSetArray object
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["conint"]] <- conint  # change the means vector from getPlotSetArray object
    }
    
  }
  file.remove(paste0(o.tmp))
  plotAverage(gpsa,xlab='Relative position [bp]', 
              ylim=ylim, ylab='Signal',
              main = paste0("Plot profile \n",sd," SD Removed"), 
              keepratio = F,error.estimates = err,
              frame.plot=F,
              cex.legend = 8,
              legend_pos = "topleft",
              pointsize = 16,
              colvec =brewer.pal(n = length(bw.l)*length(gr.v)+1 , name = "Paired"))
}

myCHR <- c("2L","2R","3L","3R","X")
kc.chr <- c("2L","2R","3L","3R","X")
s2.chr <- c("2L","2R","3L","3R","X")
sqfDM3 <- Seqinfo(genome="dm3")
# sqfDM6 <- Seqinfo(genome="dm6")
dm3.gr <- GRanges(sqfDM3)



# HG19
hela.chr <- c(seq(1,22,1),"X")
sqfHG19 <- Seqinfo(genome="hg19")
sqfHG38 <- Seqinfo(genome="hg38")

genome <- list(dm3=c(txdb=TxDb.Dmelanogaster.UCSC.dm3.ensGene,kc.chr=list(kc.chr),s2.chr=list(s2.chr),sqf=sqfDM3),
               # dm6=c(txdb=TxDb.Dmelanogaster.UCSC.dm6.ensGene,kc.chr=list(kc.chr),s2.chr=list(s2.chr),sqf=sqfDM6),
               hg19=c(txdb=TxDb.Hsapiens.UCSC.hg19.knownGene,hela.chr=list(hela.chr),sqf=sqfHG19),
               hg38=c(txdb=TxDb.Hsapiens.UCSC.hg38.knownGene,hela.chr=list(hela.chr),sqf=sqfHG38)
)
create <- function(path){
  for(idx in 1:length(path)){
    dir.create(path[idx],showWarnings = F,recursive = T)
  }
  return(path)
}
addSeqinfo <- function(my.gr,g_name="dm3",celltype="kc",keepSeq=T){
  sqf <- genome[[g_name]]$sqf
  chr <- as.character(unlist(genome[[g_name]][paste0(celltype,".chr")]))
  seqlevelsStyle(my.gr) <- "UCSC"
  seqinfo(my.gr) <- sqf[seqlevels(my.gr),]
  seqlevelsStyle(my.gr) <- "Ensembl"
  print(my.gr)
  print(chr)
  print( keepSeqlevels(my.gr,chr,pruning.mode="coarse"))
  if(keepSeq) my.gr <- keepSeqlevels(my.gr,chr,pruning.mode="coarse")
  my.gr <- sortSeqlevels(my.gr)
  my.gr <- sort(my.gr,ignore.strand=T)
}


genes.dm3.gr <- addSeqinfo(genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene),"dm3","kc")
dm3_mychr.gr <- addSeqinfo(dm3.gr,"dm3")
tss.dm3.gr <- resize(genes.dm3.gr,1,"start")

#####################################################################################-
#          LOAD DATA  ----
#####################################################################################-

h3k27wt_2013_inpNorm_log2_shft.bwp <- "/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/H3K27_WT_KD_2013/BOWTIE/BW/K27_WT_2013_Q10_sorted_dedup_SHFT135_log2_inputReadCountNorm.bw"
h3k27kd_2013_inpNorm_log2_shft.bwp <- "/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/H3K27_WT_KD_2013/BOWTIE/BW/K27_KDbeaf_2013_Q10_sorted_dedup_SHFT135_log2_inputReadCountNorm.bw"


#####################################################################################-
#
#          RUN  ----
#
#####################################################################################-

# Profile
SD <- 2
bw.l <- c(h3k27wt_2013_inpNorm_log2_shft.bwp,
          h3k27kd_2013_inpNorm_log2_shft.bwp,
          h3k27bi1_2018_inpNorm_log2_shft.bwp,
          h3k27bi2_2018_inpNorm_log2_shft.bwp,
          h3k27mbi1_2018_inpNorm_log2_shft.bwp,
          h3k27mbi2_2018_inpNorm_log2_shft.bwp,
          h3k27mbi_2013_inpNorm_log2_shft.bwp)
out_PRF.d <- create(paste0(out.d,"/PROFILE_TSS_CENTERED/"))


pdf(paste0(out_PRF.d,"PROFILE_TSS_Border_K27_multiK27_chipseq_WT_KD_normINP_shift_noshift_LOG2.pdf"))
seqPlotSDoutliers(bw.l,tmp,c("TSSborderbgLeft.gr"),ylim=c(-5,5),xlim=c(2000,2500),bin=10L,sd=SD,err = F,smooth = T,type="mf",gnme = "dm3",ignore.strand = F)
seqPlotSDoutliers(bw.l,tmp,c("TSSborderbgRight.gr"),ylim=c(-5,5),xlim=c(2000,2500),bin=10L,sd=SD,err = F,smooth = T,type="mf",gnme = "dm3",ignore.strand = F)
seqPlotSDoutliers(bw.l,tmp,c("ctl.gr"),ylim=c(-5,5),xlim=c(2000,2500),bin=10L,sd=SD,err = F,smooth = T,type="mf",gnme = "dm3",ignore.strand = F)
dev.off()

