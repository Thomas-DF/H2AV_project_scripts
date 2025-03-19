#####################################################################################
#      venn_diagram
#####################################################################################
#
# Author(s): David Depierre, Tachibana Lab
# depierre@mpg.biochem.de | dav.depierre@gmail.com
# whoami: depierre
#
# Date: 2023-03-21
# Last update: 2023-03-21
#
#####################################################################################
#
# Project: plot functions R
# petit cadeau de David, meme s'il est parti, il continue d'updater ses super fonctions de plot, dont cette magnifique fonction pour faire des venn proportinel avec un test de fisher et meme un test d'intersection multiple (ca c'est de la balle!)
# 
#
# Description: venn diagram based on vennerable
#
# R version: 
#
#####################################################################################
# LOAD LIBRARY & SOURCE FUNCTIONS
#####################################################################################

	# Set workdir 
		# workdir=''
		# setwd('')
	# Source functions 
		# funcdir='/fs/home/depierre/work/SCRIPTS/'
	    # source(paste0(funcdir, 'script_repo/r_scripts/plot_function/avg_profile/seqplots_avg_prof_on_R_4.1.3_env.R'))
	    # source(paste0(funcdir, 'script_repo/r_scripts/...'))
	# Set personnal library path 
		# .libPaths(c(.libPaths(), '/path/to/'))
		# envdir='' # if renv library is used
	# Load library 
	    # library('')

#####################################################################################
# LOAD DATA
#####################################################################################



######################################################################################################################################################
## LIBRARY
######################################################################################################################################################
library(gplots)
library(Vennerable)
library(SuperExactTest)
library(RColorBrewer)
######################################################################################################################################################
## FUNCTION
######################################################################################################################################################


plotVENNerable = function(
    input.lst,
    col.vec = c("#0000ff", "#6d0000", "#cc9900"),
    fisher.eff.tot = NULL,
    doWeights = TRUE,
    info=NULL,
    outdir
    ){


    vennNAME = paste(names(input.lst), collapse="_")
    Venn_ovlp = Vennerable::Venn(input.lst)
    if(is.null(fisher.eff.tot)){
        fisher.eff.tot = length(unique(do.call(c,input.lst)))
    }

    catched.msg = evaluate::try_capture_stack(Vennerable::VennThemes(Vennerable::compute.Venn(Venn_ovlp, doWeights = doWeights)), env = parent.frame())
    if(grepl('Error', as.character(catched.msg))){
        doWeights = FALSE
    }
    color_Venn = Vennerable::VennThemes(Vennerable::compute.Venn(Venn_ovlp, doWeights = doWeights))

    #prepare color sets
    while(length(col.vec)<3){
        col.vec=c(col.vec, "#000000")
    }
    col.vec2 = lighten_copied(col.vec, 0.4)
    combn.mtx = combn(seq(col.vec2), m=2)
    col.combn = c()
    for(i in seq(col.vec2)){
        col.combn = c(col.combn,
            colorspace::mixcolor(alpha = 0.5, colorspace::sRGB(c(col2rgb(col.vec2[combn.mtx[1,i]]))[1], c(col2rgb(col.vec2[combn.mtx[1,i]]))[2], c(col2rgb(col.vec2[combn.mtx[1,i]]))[3]),
                                        colorspace::sRGB(c(col2rgb(col.vec2[combn.mtx[2,i]]))[1], c(col2rgb(col.vec2[combn.mtx[2,i]]))[2], c(col2rgb(col.vec2[combn.mtx[2,i]]))[3])))
    }
    col.vec3 = c(
        rgb(c(col.combn[[1]]@coords)[1]/255, c(col.combn[[1]]@coords)[2]/255, c(col.combn[[1]]@coords)[3]/255),
        rgb(c(col.combn[[2]]@coords)[1]/255, c(col.combn[[2]]@coords)[2]/255, c(col.combn[[2]]@coords)[3]/255),
        rgb(c(col.combn[[3]]@coords)[1]/255, c(col.combn[[3]]@coords)[2]/255, c(col.combn[[3]]@coords)[3]/255)
    )
    if(length(input.lst)>2){
            col.vec4 = colorspace::mixcolor(alpha = 0.5, colorspace::sRGB(c(col2rgb(col.vec3[1]))[1], c(col2rgb(col.vec3[1]))[2], c(col2rgb(col.vec3[1]))[3]),
                                        colorspace::sRGB(c(col2rgb(col.vec3[3]))[1], c(col2rgb(col.vec3[3]))[2], c(col2rgb(col.vec3[3]))[3]))
            col.vec4 = rgb(c(col.vec4@coords)[1]/255, c(col.vec4@coords)[2]/255, c(col.vec4@coords)[3]/255)
    }
    if(length(input.lst)==2){
        #set1
        color_Venn[["Set"]][["Set1"]]$col = col.vec[1]
        color_Venn[["SetText"]][["Set1"]]$col = col.vec[1]
        color_Venn[["Face"]][["10"]]$fill = col.vec2[1]
        color_Venn[["Face"]][["10-1"]]$fill = col.vec2[1]
        #set2
        color_Venn[["Set"]][["Set2"]]$col = col.vec[2]
        color_Venn[["SetText"]][["Set2"]]$col = col.vec[2]
        color_Venn[["Face"]][["01"]]$fill = col.vec2[2]
        color_Venn[["Face"]][["01-1"]]$fill = col.vec2[2]
        #set1 n set2
        color_Venn[["Face"]][["11"]]$fill = col.vec3[1]
        color_Venn[["Face"]][["11-1"]]$fill = col.vec3[1]

    }
    if(length(input.lst)==3){
        #set1
        color_Venn[["Set"]][["Set1"]]$col = col.vec[1]
        color_Venn[["SetText"]][["Set1"]]$col = col.vec[1]
        color_Venn[["Face"]][["100"]]$fill = col.vec2[1]
        color_Venn[["Face"]][["100-1"]]$fill = col.vec2[1]
        #set2
        color_Venn[["Set"]][["Set2"]]$col = col.vec[2]
        color_Venn[["SetText"]][["Set2"]]$col = col.vec[2]
        color_Venn[["Face"]][["010"]]$fill = col.vec2[2]
        color_Venn[["Face"]][["010-1"]]$fill = col.vec2[2]
        #set3
        color_Venn[["Set"]][["Set3"]]$col = col.vec[3]
        color_Venn[["SetText"]][["Set3"]]$col = col.vec[3]
        color_Venn[["Face"]][["001"]]$fill = col.vec2[3]
        color_Venn[["Face"]][["001-1"]]$fill = col.vec2[3]
        #set1 n set2
        color_Venn[["Face"]][["110"]]$fill = col.vec3[1]
        color_Venn[["Face"]][["110-1"]]$fill = col.vec3[1]
        #set1 n set3
        color_Venn[["Face"]][["101"]]$fill = col.vec3[2]
        color_Venn[["Face"]][["101-1"]]$fill = col.vec3[2]
        #set2 n set3
        color_Venn[["Face"]][["011"]]$fill = col.vec3[3]
        color_Venn[["Face"]][["011-1"]]$fill = col.vec3[3]
        color_Venn[["Face"]][["111"]]$fill = col.vec4
        color_Venn[["Face"]][["111-1"]]$fill = col.vec4 
    }

    #fisher exact test
    fisher.res.lst = fisherExactTest_onVecList(input.lst,input.lst,fisher.eff.tot)
    # Supertest for triple intersection
    ResSuperTest = SuperExactTest::supertest(input.lst, n = fisher.eff.tot)
    SupetTestToPrint = as.matrix(ResSuperTest$P.value)
    colnames(SupetTestToPrint) = "intersect_pval"

    ########### PLOT #########################
    pdf(paste0(outdir, "PLOTvenn_",vennNAME,"_",info,".pdf"))
        plot(Venn_ovlp, gp = color_Venn, doWeights = doWeights)
        par(mfrow=c(2,1))
        textplot(fisher.res.lst$p.value,  valign="top")
        textplot(fisher.res.lst$odds.ratio)
        textplot(SupetTestToPrint)
        if(!doWeights){
            textplot("Can't do proportional...")
        }
        plot(Venn_ovlp, gp = color_Venn, doWeights = doWeights, show = list(SetLabels = FALSE, FaceText = ""))
    dev.off()

}

#####################################################################################
# LOAD DATA
#####################################################################################

fisherExactTest_onVecList=function(list1, list2, total, test.side = "greater"){
    ## init output results list
        res.lst = list()
        res.lst$p.value=matrix(nrow=length(list1),ncol=length(list2))
        res.lst$odds.ratio=matrix(nrow=length(list1),ncol=length(list2))
        res.lst$effectives=matrix(nrow=length(list1),ncol=length(list2))
        res.lst$pct_list1=matrix(nrow=length(list1),ncol=length(list2))
        res.lst$pct_list2=matrix(nrow=length(list1),ncol=length(list2))
        colnames(res.lst$p.value)=names(list2)
        rownames(res.lst$p.value)=names(list1)
        colnames(res.lst$odds.ratio)=names(list2)
        rownames(res.lst$odds.ratio)=names(list1)
        colnames(res.lst$effectives)=names(list2)
        rownames(res.lst$effectives)=names(list1)
        colnames(res.lst$pct_list1)=names(list2)
        rownames(res.lst$pct_list1)=names(list1)
        colnames(res.lst$pct_list2)=names(list2)
        rownames(res.lst$pct_list2)=names(list1)

    ## exec fisherexact test and fill matrices with results, effectives or percentage
        for (i in 1:length(list1)){
            liste1=list1[[i]]
            for (j in 1:length(list2)){
                liste2=list2[[j]]
                inter=length(Reduce(intersect, list(liste1,liste2)))
                res.lst$p.value[i,j]=signif(fisher.test(matrix(c(inter,length(liste1)-inter,length(liste2)-inter,max(total-length(liste1)-length(liste2)+inter,0)),ncol=2),alternative=test.side)$p.value,2)
                res.lst$odds.ratio[i,j]=signif(fisher.test(matrix(c(inter,length(liste1)-inter,length(liste2)-inter,max(total-length(liste1)-length(liste2)+inter,0)),ncol=2),alternative=test.side)$estimate,2)
                res.lst$effectives[i,j]=inter
                res.lst$pct_list1[i,j] = inter/length(liste1)
                res.lst$pct_list2[i,j] = inter/length(liste2)
            }
        }
        return(res.lst)
}

# https://rdrr.io/rforge/colorspace/src/R/lighten.R#sym-lighten
lighten_copied <- function(col, amount = 0.1,
                    method = c("relative", "absolute"), space = c("HCL", "HLS", "combined"), fixup = TRUE)
{
  ## method
  space <- match.arg(space, c("HCL", "HLS", "combined"))
  method <- match.arg(method, c("relative", "absolute"))
  
  ## number of colors
  n <- max(c(length(col), length(amount)))
  col <- rep_len(col, length.out = n)
  amount <- rep_len(amount, length.out = n)
  
  ## save original colors for later, to substitute any cases with amount == 0
  col_orig <- col
  
  ## col has to be hex code, otherwise col2rgb is used
  if(is.character(col) &&
     (all(substr(col, 1L, 1L) == "#") & all(nchar(col) %in% c(7L, 9L))))
  {
    ## extract alpha from hex (if any)
    alpha <- substr(col, 8L, 9L)
    ## retain only RGB in hex
    col <- substr(col, 1L, 7L)
    ## convert to colorspace::RGB
    col <- hex2RGB(col)
  } else {
    col <- col2rgb(col, alpha = TRUE)
    # save original colors in hex format, in case some were specified as 
    # named colors or as palette entries
    col_orig <- rgb(t(col), maxColorValue = 255)
    ## extract alpha values (if non-FF)
    alpha <- format(as.hexmode(col[4L, ]), width = 2L, upper.case = TRUE)
    alpha[alpha == "FF"] <- ""
    ## retain only sRGB
    col <- sRGB(t(col[1L:3L, ])/255)
  }
  
  if (space == "HCL") {
    ## *** darkening/lightening in HCL space ***
    
    ## convert to HCL
    col <- as(col, "polarLUV")
    ## fix-up extreme luminance cases
    col@coords[, "L"] <- pmin(100, pmax(0, col@coords[, "L"]))

    ## adjust luminance
    Lold <- col@coords[, "L"]
    col@coords[, "L"] <- if(method == "relative") {
      (amount >= 0) * (100 - (100 - Lold) * (1 - amount)) +
        (amount < 0) * Lold * (1 + amount)
    } else {
      Lold + amount * 100
    }
    col@coords[, "L"] <- pmin(100, pmax(0, col@coords[, "L"]))
    
    ## transform chroma correspondingly (relative to maximum chroma possible)
    ## It seems better to not apply this adjustment here. Lighened colors look better without it, 
    ## and darkened colors look better under the combined model.
    #col@coords[, "C"] <- col@coords[, "C"]/ceiling(max_chroma(col@coords[, "H"], Lold) + 1e-8) *
    #  max_chroma(col@coords[, "H"], col@coords[, "L"], floor = TRUE)
  
    ## check that resulting chroma is within appropriate bounds
    col@coords[, "C"] <- pmin(max_chroma(col@coords[, "H"], col@coords[, "L"], floor = TRUE),
                              pmax(0, col@coords[, "C"]))
  } 
  else if (space == "HLS") {
    ## *** darkening/lightening in HLS space ***
    
    col <- as(col, "HLS")
    col@coords[, "L"] <- if(method == "relative") {
      (amount >= 0) * (1 - (1 - col@coords[, "L"]) * (1 - amount)) +
        (amount < 0) * col@coords[, "L"] * (1 + amount)
    } else {
      col@coords[, "L"] + amount
    }
    col@coords[, "L"] <- pmin(1, pmax(0, col@coords[, "L"]))
  } else {
    ## *** darkening/lightening in combined space ***
    
    ## first do adjustment in HLS space
    colHLS <- as(col, "HLS")
    colHLS@coords[, "L"] <- if(method == "relative") {
      (amount >= 0) * (1 - (1 - colHLS@coords[, "L"]) * (1 - amount)) +
        (amount < 0) * colHLS@coords[, "L"] * (1 + amount)
    } else {
      colHLS@coords[, "L"] + amount
    }
    colHLS@coords[, "L"] <- pmin(1, pmax(0, colHLS@coords[, "L"]))
    
    colHLSHCL <- as(as(colHLS, "RGB"), "polarLUV")
    
    ## now do adjustment in HCL space
    col <- as(col, "polarLUV")
    ## fix-up extreme luminance cases
    col@coords[, "L"] <- pmin(100, pmax(0, col@coords[, "L"]))

    ## transform luminance
    Lold <- col@coords[, "L"]
    col@coords[, "L"] <- if(method == "relative") {
      (amount >= 0) * (100 - (100 - Lold) * (1 - amount)) +
        (amount < 0) * Lold * (1 + amount)
    } else {
      Lold + amount * 100
    }
    
    ## fix-up L and copy C over from HLS-converted color
    col@coords[, "L"] <- pmin(100, pmax(0, col@coords[, "L"]))
    #col@coords[, "H"] <- colHLSHCL@coords[, "H"]
    col@coords[, "C"] <- colHLSHCL@coords[, "C"]
    
    ## make sure chroma is in allowed range
    col@coords[, "C"] <- pmin(max_chroma(col@coords[, "H"], col@coords[, "L"], floor = TRUE), col@coords[, "C"])
  }
  
  ## convert back to hex and add alpha again (if any)
  col <- hex(col, fixup = fixup)
  col[!is.na(col)] <- paste(col[!is.na(col)], alpha[!is.na(col)], sep = "")
  
  ## return original colors whenever amount == 0
  col[amount == 0] <- col_orig[amount == 0]
  
  return(col)
}
# https://rdrr.io/rforge/colorspace/src/R/max_chroma.R
max_chroma <- function(h, l, floor = FALSE) {
  ## align h and l
  n <- max(c(length(h), length(l)))
  h <- rep_len(h, n)
  l <- rep_len(l, n)

  ## assure h in [0, 360]
  while(any(h < 0)) h[h < 0] <- h[h < 0] + 360
  while(any(h >= 360)) h[h >= 360] <- h[h >= 360] - 360

  ## assure l in [0, 100]
  l <- pmin(100, pmax(0, l))

  ## obtain surrounding h/l coordinates
  hmin <- floor(h + 1e-8)
  hmax <- ceiling(h + 1e-8)
  lmin <- floor(l + 1e-8)
  lmax <- ceiling(l + 1e-8)

  ## average
  c <- (hmax - h) * (lmax - l) * colorspace::max_chroma_table[paste(hmin, lmin, sep = "-")] + 
       (hmax - h) * (l - lmin) * colorspace::max_chroma_table[paste(hmin, lmax, sep = "-")] + 
       (h - hmin) * (lmax - l) * colorspace::max_chroma_table[paste(hmax, lmin, sep = "-")] + 
       (h - hmin) * (l - lmin) * colorspace::max_chroma_table[paste(hmax, lmax, sep = "-")]

  ## catch border cases
  c <- as.numeric(c) # pmin(c, 100)
  c[l <= 0 | l >= 100] <- 0
  
  ## take floor to be "on the safe side"
  if(floor) c <- floor(c)
  
  return(c)
}


#####################################################################################
# END
#####################################################################################
