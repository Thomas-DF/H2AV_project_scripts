################################################################
# PLOT FUNCTION LIBRARY
# David Depierre
# dav.depierre@gmail.com
################################################################
# BOXPLOT
################################################################

################################################################
# library
################################################################
library(rstatix)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gplots)
library(viridis)
'%ni%' = Negate('%in%')



# source(paste0("/home/ddepierre/work/script_repo/utils_function/utils_function.R"))
################################################################
# function
################################################################

                # get_density for scatter plot
                    get_density <- function(x, y, n = 100) {
                        dens <- MASS::kde2d(x = x, y = y, n = n)
                        ix <- findInterval(x, dens$x)
                        iy <- findInterval(y, dens$y)
                        ii <- cbind(ix, iy)
                        return(dens$z[ii])
                    }



                        feat1.vec = Q_Zj12_shC_RBM26_TSS
                        feat2.vec = Q_Zj11_shC_PapolG_TSS
                        namefeat1 = "Q_Zj12_shC_RBM26_TSS"
                        namefeat2 = "Q_Zj11_shC_PapolG_TSS"
                        subgroup.lst = list(
                            topPol2 = getNameList(Q_Zj7_shC_PolII_TSS, topdown = "top", prct = 50),
                            downPol2 = getNameList(Q_Zj7_shC_PolII_TSS, topdown = "down", prct = 50)
                        )
                        colorsubgroup = c("#f5a623ff", "#8f8f8f")
                        logT=T
                        yl=c(0,15)
                        xl=c(0,15)
                        abline.type="diagonal"
                        info = NULL
                        outDir = "/home/cperrois/work/PROJET_MTREC_IGH/FIGURE/SCATTER/SCATTER_rev/SCATTER_highlow_pol2/"



                # scatterPlot.func
                    scatterPlot.v230409 = function(
                        feat1.vec,
                        feat2.vec,
                        namefeat1,
                        namefeat2,
                        subgroup.lst = NULL,
                        colorsubgroup = NULL,
                        logT=F,
                        yl=c(-1,1),
                        xl=c(-1,1),
                        abline.type="diagonal",
                        info = NULL,
                        outDir
                    ){
                        
                        # Helper
                            #Na. scatterPlot.func
                            #De. Function that scatterPlot.func
                            #Us. scatterPlot.func(feat1.vec, feat2.vec,  namefeat1, namefeat2, logT=F, info = NULL, outDir, yl=c(-1,1), xl=c(-1,1)
                            #Ar. abline.type= "diagonal" or "horizontal" or "vertical" or NULL
                            #Va. None

                            #Ex. 
                            #Ex. # EXEMPLE 1 : 
                            #Ex. A = rnorm(2000, 0, 1)
                            #Ex. names(A) = paste0("quantif", 1:length(A))
                            #Ex. B=jitter(A, factor = 2, amount = 1)
                            #Ex. 
                            #Ex. scatterPlot.v220622(feat1.vec = A, feat2.vec = B,  namefeat1 = "someQuantifA", namefeat2 = "someQuantifB",
                            #Ex.                     subgroup.lst = NULL, colorsubgroup = NULL, logT=F, abline.type="diagonal", yl=c(-5,5), xl=c(-5,5),
                            #Ex.                     outDir = "/home/ddepierre/work/script_repo/plot_function/scatterplot/", info = "exemple_1")
                            #Ex. 
                            #Ex. # EXEMPLE 2 : 
                            #Ex. A = rnorm(2000, 0, 1)
                            #Ex. names(A) = paste0("quantif", 1:length(A))
                            #Ex. B=jitter(A, factor = 2, amount = 2)
                            #Ex. grp1 = sample(names(B), 500)
                            #Ex. B[grp1] = B[grp1]*2
                            #Ex. 
                            #Ex. subgroup.lst = list(
                            #Ex.   grp1 = grp1,
                            #Ex.   grp2 = names(B)[names(B) %ni% grp1]
                            #Ex. )
                            #Ex. 
                            #Ex. scatterPlot.v220622(feat1.vec = A, feat2.vec = B,  namefeat1 = "someQuantifA", namefeat2 = "someQuantifB",
                            #Ex.                     subgroup.lst = subgroup.lst, colorsubgroup = NULL, logT=F, abline.type="diagonal", yl=c(-5,5), xl=c(-5,5),
                            #Ex.                     outDir = "/home/ddepierre/work/script_repo/plot_function/scatterplot/", info = "exemple_2")
                            #Ex. 

                        if(is.null(names(feat1.vec))){
                            names(feat1.vec) = paste0("row_", seq(1,length(feat1.vec),1))
                        }
                        if(is.null(names(feat2.vec))){
                            names(feat2.vec) = paste0("row_", seq(1,length(feat2.vec),1))
                        }
                        if(!is.null(subgroup.lst)){
                            if(sum(duplicated(unlist(subgroup.lst))) > 0){ # remove duplicated elements
                                print((paste0("WARNING : ",sum(duplicated(unlist(subgroup.lst)))," duplicated elements of subgroup.lst were removed")))
                                dupElement = unlist(subgroup.lst)[duplicated(unlist(subgroup.lst))]
                                subgroup.lst = lapply(subgroup.lst, function(vec){vec = vec[vec %ni% dupElement]})
                            }
                            if(is.null(colorsubgroup) | length(colorsubgroup)<length(subgroup.lst)){ # if no color subgroup defined then define color
                                colorsubgroup = rainbow(length(subgroup.lst))
                                print((paste0("colorsubgroup created : ",paste0(colorsubgroup, collapse=", "))))

                            }
                        }

                        attribtoPlot = cbind(attributes(feat1.vec)[unlist(lapply(attributes(feat1.vec), length))<4],
                                             attributes(feat2.vec)[unlist(lapply(attributes(feat2.vec), length))<4])

                        feat1.vec = feat1.vec[!is.na(feat1.vec)]
                        feat2.vec = feat2.vec[!is.na(feat2.vec)]


                        # if logT %in% T, value in .vec will be log2 transformed (usefull when plotting -seq quantif)
                        if(logT %in% T){
                          feat1.vec = log2(feat1.vec+1)
                          feat2.vec = log2(feat2.vec+1)
                        }
                        # transform .vec en matrix
                        feat1 = as.matrix(feat1.vec)
                        feat2 = as.matrix(feat2.vec)

                        # define a set of common rownames + filter by this set of names (i.e. only rows with values in feat1 and feat2 will be plotted)
                    	comRownames = Reduce(intersect, list(rownames(feat1),rownames(feat2)))
                    	feat1 = feat1[comRownames,1]
                    	feat2 = feat2[comRownames,1]
                        if(!is.null(subgroup.lst)){
                            subgroup = rep(NA, length(feat1))
                            lapply(seq_along(subgroup.lst), function(i){
                                subgroup[names(feat1) %in% subgroup.lst[[i]]] <<- names(subgroup.lst)[i]
                            })
                            subgroup = factor(subgroup, levels = names(subgroup.lst))

                             # create a data.frame
                            data=data.frame(feat1=feat1, feat2=feat2, subgroup=subgroup)
                        }else{
                            data=data.frame(feat1=feat1, feat2=feat2)
                        }
                        data.filt = data[rownames(data) %in% unlist(subgroup.lst),]

                        ## STAT TEST FOR CORRELATION
                        ## PEARSON TEST
                    	pearson = cor.test(data$feat1, data$feat2, na.rm=T, method = "pearson")
                    	pearsonT = round(pearson$estimate, 2)
                        pearsonP = pearson$p.value

                        # PEARSON for each subgroup
                        if(!is.null(subgroup.lst)){
                            PEARSON.subgroup = list()
                            for(i in seq_along(subgroup.lst)){
                                data.sub = data[data$subgroup %in% names(subgroup.lst)[i],]
                                PEARSON.subgroup[[i]] = data.frame(pearson.estimate = round(cor.test(data.sub$feat1, data.sub$feat2, na.rm=T, method = "pearson")$estimate,2),
                                            pearson.pval = cor.test(data.sub$feat1, data.sub$feat2, na.rm=T, method = "pearson")$p.value)
                                names(PEARSON.subgroup)[i] = names(subgroup.lst)[i]
                            }
                        }



                        ## SPEARMAN TEST
                    	spearman = cor.test(data$feat1, data$feat2, na.rm=T, method = "spearman")
                    	spearmanT = round(spearman$estimate, 2)
                        spearmanP = spearman$p.value
                        ## fit a classic linear model
                        res.lm = lm(feat2~feat1)
                        print(summary(res.lm))
                        summary(res.lm)$r.squared

                        ## PLOT scatter with R ggplot functions (harder to code but fancy and beautifull plots)
                        # get data density for ggplot density option, data is the data.frame create with the 2 vetors, that is used in the ggplot functions
                        data$density <- get_density(data$feat1,data$feat2, n=500)



                        ## Basic ggplot scatter with density
                        GGPLOT = ggplot(data = data) + geom_point(aes(x=feat1, y=feat2, color = density)) +
                                 scale_color_viridis() + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl) +
                                 theme_classic()
                        ## Basic ggplot scatter with density + fit linear model
                        GGPLOTlm = ggplot(data = data, aes(x=feat1, y=feat2, color = density)) + geom_point() +
                                 scale_color_viridis() + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl) +
                                 stat_smooth() + theme_classic()


                        ## Basic ggplot scatter with density + fit none linear model
                        GGPLOT_nlm = ggplot(data = data, aes(x=feat1, y=feat2, color = density)) + geom_point() +
                                 scale_color_viridis() + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl) +
                                 stat_smooth(method="lm") + theme_classic()

                        if(!is.null(subgroup.lst)){
                            GGPLOT_colorSubGrp <- ggplot(data = data, aes(x=feat1, y=feat2)) + geom_point(aes(color = subgroup),alpha = 0.5) + 
                                scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl) +
                                scale_color_manual(values = c(colorsubgroup, "#8f8f8f")) + theme_classic()
                            GGPLOT_colorSubGrpfilt <- ggplot(data = data.filt, aes(x=feat1, y=feat2)) + geom_point(aes(color = subgroup),alpha = 0.5) + 
                                scale_color_manual(values = c(colorsubgroup)) + scale_fill_manual(values = c(colorsubgroup)) + 
                                scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl) + theme_classic()
                            GGPLOT_colorSubGrpfilt_lm <- ggplot(data = data.filt, aes(x=feat1, y=feat2)) + geom_point(aes(color = subgroup),alpha = 0.5) + 
                                stat_smooth(aes(color=subgroup, fill=subgroup), method="lm", alpha = 0.2) +
                                scale_color_manual(values = c(colorsubgroup)) + scale_fill_manual(values = c(colorsubgroup)) + 
                                scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl) + theme_classic()
                            
                            GGPLOT_SubGrpfilt.lst = list()
                            lapply(seq_along(subgroup.lst), function(i){
                                data.filt2 = data[rownames(data) %in% subgroup.lst[[i]],]
                            GGPLOT_SubGrpfilt.lst[[i]] <<- ggplot(data = data.filt2, aes(x=feat1, y=feat2)) + geom_point(color = colorsubgroup[i]) +
                                scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl) + theme_classic()
                            GGPLOT_SubGrpfilt.lst[[i+length(subgroup.lst)]] <<- ggplot(data = data.filt2, aes(x=feat1, y=feat2)) + geom_point(color = colorsubgroup[i]) +
                                scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl) + theme_classic()+ 
                                stat_smooth(color = colorsubgroup[i], method="lm", alpha = 0.4)
                            })
                        
                        
                        
                        }


                        if(!is.null(abline.type)){
                            if(abline.type %in% "diagonal"){
                            GGPLOT = GGPLOT + geom_abline(intercept = 0, slope = 1, alpha=0.5)
                            GGPLOTlm = GGPLOTlm + geom_abline(intercept = 0, slope = 1, alpha=0.5)
                            GGPLOT_nlm = GGPLOT_nlm + geom_abline(intercept = 0, slope = 1, alpha=0.5)
                            }else if(abline.type %in% "horizontal"){
                            GGPLOT = GGPLOT + geom_hline(yintercept = 0)
                            GGPLOTlm = GGPLOTlm + geom_hline(yintercept = 0)
                            GGPLOT_nlm = GGPLOT_nlm + geom_hline(yintercept = 0)
                            }else if(abline.type %in% "vertical"){
                            GGPLOT = GGPLOT + geom_vline(xintercept = 0)
                            GGPLOTlm = GGPLOTlm + geom_vline(xintercept = 0)
                            GGPLOT_nlm = GGPLOT_nlm + geom_vline(xintercept = 0)
                            }
                        }

                        if(!is.null(subgroup.lst)){
                            if(!is.null(abline.type)){
                                if(sum(abline.type %in% "diagonal")>=1){
                                GGPLOT_colorSubGrp = GGPLOT_colorSubGrp + geom_abline(intercept = 0, slope = 1, alpha=0.5)
                                GGPLOT_colorSubGrpfilt = GGPLOT_colorSubGrpfilt + geom_abline(intercept = 0, slope = 1, alpha=0.5)
                                GGPLOT_colorSubGrpfilt_lm = GGPLOT_colorSubGrpfilt_lm + geom_abline(intercept = 0, slope = 1, alpha=0.5)
                                }
                                if(sum(abline.type %in% "horizontal")>=1){
                                GGPLOT_colorSubGrp = GGPLOT_colorSubGrp + geom_hline(yintercept = 0)
                                GGPLOT_colorSubGrpfilt = GGPLOT_colorSubGrpfilt + geom_hline(yintercept = 0)
                                GGPLOT_colorSubGrpfilt_lm = GGPLOT_colorSubGrpfilt_lm + geom_hline(yintercept = 0)
                                }
                                if(sum(abline.type %in% "vertical")>=1){
                                GGPLOT_colorSubGrp = GGPLOT_colorSubGrp + geom_vline(xintercept = 0)
                                GGPLOT_colorSubGrpfilt = GGPLOT_colorSubGrpfilt + geom_vline(xintercept = 0)
                                GGPLOT_colorSubGrpfilt_lm = GGPLOT_colorSubGrpfilt_lm + geom_vline(xintercept = 0)
                                }
                                lapply(seq_along(GGPLOT_SubGrpfilt.lst), function(i){
                                    if(sum(abline.type %in% "diagonal")>=1){
                                    GGPLOT_SubGrpfilt.lst[[i]] <<- GGPLOT_SubGrpfilt.lst[[i]] + geom_abline(intercept = 0, slope = 1, alpha=0.5)
                                    }
                                    if(sum(abline.type %in% "horizontal")>=1){
                                    GGPLOT_SubGrpfilt.lst[[i]] <<- GGPLOT_SubGrpfilt.lst[[i]] + geom_hline(yintercept = 0)
                                    }
                                    if(sum(abline.type %in% "vertical")>=1){
                                    GGPLOT_SubGrpfilt.lst[[i]] <<- GGPLOT_SubGrpfilt.lst[[i]] + geom_vline(xintercept = 0)
                                    }
                                })
                            }
                        }

                        ## open the pdf to plot inside
                       	pdf(paste0(outDir,"scatter_",namefeat1,"_VS_", namefeat2, "_",info,".pdf"))
                                # PLOT scatter with R "base" function
                                # free axis scales
                      	        plot(feat1, feat2, pch=20, xlab = namefeat1, ylab = namefeat2,
                                        main=paste0(namefeat1 , " & ", namefeat2, " free axis scales"))
                                        abline(lm(feat2~feat1), col="red") # fit linear model on dots and plot it in red
                                        lines(x = c(-1000,1000), y = c(-1000,1000)# plot diagonal in black
                                )

                                # defined axis scales
                                plot(feat1, feat2, pch=20,ylim=yl, xlim=xl, xlab = namefeat1, ylab = namefeat2,
                                        main=paste0(namefeat1 , " & ", namefeat2, " defined axis scales"))
                                        abline(lm(feat2~feat1), col="red") # fit linear model on dots and plot it in red
                                        lines(x = c(-1000,1000), y = c(-1000,1000)# plot diagonal in black
                                ) 

                                # PRINT ggplot 
                                print(GGPLOT)
                                print(GGPLOTlm)
                                print(GGPLOT_nlm)
                                if(!is.null(subgroup.lst)){
                                    print(GGPLOT_colorSubGrp)
                                    print(GGPLOT_colorSubGrpfilt)
                                    print(GGPLOT_colorSubGrpfilt_lm)
                                    lapply(GGPLOT_SubGrpfilt.lst, print)
                                }

 
                                ## WRITE TEST RESULTS IN THE PDF
                                textplot(pearson, cex=1)
                                # PEARSON for each subgroup
                                if(!is.null(subgroup.lst)){
                                    textplot(PEARSON.subgroup)
                                }


                                textplot(spearman, cex=1)

                                textplot(paste0("R squared de lm() :  ", summary(res.lm)$r.squared), cex=1)

                                # textplot(cbind(attributes(apa.mtx)[unlist(lapply(attributes(apa.mtx), class)) %ni% c("matrix", "data.frame")]))

                        # close the pdf
                    	dev.off()

                    }





################################################################
# Ntile_split
################################################################

        Ntile_split = function(vec, n, decreasing=T, return="values"){
                # Helper
                    #Na. Ntile_split
                    #De. split a vector into a list of vector
                    #Us. Ntile_split(A, n=5, return="values")
                    #Ar. vec: a numeric vector 
                    #Ar. n: an integer, number of groups wanted
                    #Ar. decreasing: should the numeric vector be ordered degreasingly or not
                    #Ar. return; "values" or "names" are returned
                    #Va. list of vector either values or names
                    #Ex. 
                    #Ex. A = rnorm(200, 0, 1)
                    #Ex. names(A) = paste0("quantif", 1:length(A))
                    #Ex. 
                    #Ex. Ntile_split(vec = A, n=5, return="values")
                    #Ex. 
                    #Ex. Ntile_split(vec = A, n=10, decreasing=F, return="names")
                    #Ex. 
        
                if(return %ni% "values" & return %ni% "names"){
                    stop("ERROR: 'return' should be either 'values' or 'names' ")
                }
        
                vec = sort(vec, decreasing=decreasing)
        
                if(return == "values"){
                    vec.lst = split(vec, ceiling(seq_along(vec)/ceiling(length(vec)/n)))
                }else if(return == "names"){
                    vec.lst = split(names(vec), ceiling(seq_along(names(vec))/ceiling(length(names(vec))/n)))
                }
                return(vec.lst)
        }


################################################################
# getNameList
################################################################

        getNameList = function(vec, topdown = "top", prct = 10){
            # Helper
                #Na. getNameList
                #De. get names of element in a vector fitlered on top, mid or down percentage
                #Us. Ntile_split(A, n=5, return="values")
                #Ar. vec: a numeric vector 
                #Ar. topdown: "top", "down", "mid" for topn down or middle part of the vector
                #Ar. prct: numeric ; which pecentage of the vector
                #Va. vector of names
                #Ex. 
                #Ex. A = rnorm(200, 0, 1)
                #Ex. names(A) = paste0("quantif", 1:length(A))
                #Ex. 
                #Ex. getNameList(vec = A, topdown = "top", prct = 10 )
                #Ex. 
            if(topdown %ni% "top" & topdown %ni% "down" & topdown %ni% "mid"){
                stop("ERROR: 'topdown' should be either 'top', 'down' or 'mid' ")
            }
            if(sum(is.na(vec))>0){
              cat(paste0("WARNING: ",sum(is.na(vec))," NA &/or NaN values are not considered \n"))
              vec = vec[!is.na(vec)]
            }
            vec = sort(vec, decreasing=T)
            if(topdown %in% "top"){
              GN = names(vec[vec > quantile(vec, (100-prct)/100, na.rm=T)])
            }
            if(topdown %in% "down"){
               GN = names(vec[vec < quantile(vec, (prct)/100, na.rm=T)])
            }
            if(topdown %in% "mid"){
               tmp1 = names(vec[vec < quantile(vec, (100/2-prct/2)/100, na.rm=T)])
               tmp2 = names(vec[vec < quantile(vec, (100/2-prct/2+prct)/100, na.rm=T)])
               GN = tmp2[tmp2 %ni% tmp1]
            }
            return(GN)
        }





# boxplot.v220617 = function(data.veclst, data.name = "nameVec", subgroup.lst = NULL, subgroup.name = "nameVec_GN", effMax=F, colorsubgroup = NULL, ylim = c(-1,1), outlier.pct=0.01, test.side = "two.sided", outDir = paste0(workdir, "path/to/out/"), info = NULL){
#     # Helper
#       #Na. histogramdensityPlot
#       #De. Function that do histogram and density plot from a vector or a list of vector
#       #De. Can do multiple histogram from different subgroup
#       #Us. scatterPlot.func(feat1.vec, feat2.vec,  namefeat1, namefeat2, logT=F, info = NULL, outDir, yl=c(-1,1), xl=c(-1,1)
#       #Ar. data = either a vector or a list of vector
#       #Ar. data.name = name of the data / metric. sue for generating plot name
#       #Ar. subgroup.lst = used only if data.veclst is a vector // list of vector of names, used to create subgroup to plot multiple histogram
#       #Ar. colorsubgroup = c("#3d3d3d", "#820002"), if null, defined from rainbow palette
#       #Ar. xlim = c(-10,10)
#       #Ar. filt_xlim [T or F], is the density plot made with values filtered by xlim or not (could affect the density if some values are cropped)
#       #Ar. ylim = c(0,0.1)
#       #Ar. ylimDiff ?
#       #Ar. binW: bandwith/binWitdh for density 
#       #Ar. alphaTransp: transparency of colors when multiple elements
#       #Ar. diffSpanSmth : Span of smooth o differential density / only in diff density on 2 elements
#       #Ar. info: additional info on plot name 
#       #Ar. outDir: path to output pdf
#       #Va. None
#       #Ex. 
#       #Ex. 
#       #Ex. # EXEMPLE 1 : 
#       #Ex. A = rnorm(200, 0, 1)
#       #Ex. names(A) = paste0("quantif", 1:length(A))
#       #Ex. data_vector = A
#       #Ex. subgroup.lst = Ntile_split(A, n=5, decreasing=T, return="names")
#       #Ex. 
#       #Ex. boxplot.v220617(data.veclst = data_vector, data.name = "some_quantifications", subgroup.lst = subgroup.lst, subgroup.name = "some_groups",
#       #Ex.               effMax=F, colorsubgroup = NULL, ylim = c(-7,7),
#       #Ex.               outDir = "/home/ddepierre/work/script_repo/plot_function/boxplot/", info = "exemple_1")
#       #Ex. 
#       #Ex. # EXEMPLE 2 : 
#       #Ex. B = rnorm(200, 2, 1.5)
#       #Ex. names(B) = paste0("quantif", 1:length(B))
#       #Ex. data2groups = list(A=A, B=B)
#       #Ex. 
#       #Ex. boxplot.v220617(data.veclst = data2groups, data.name = "some_quantifications", subgroup.lst = NULL, subgroup.name = "no_subgroups",
#       #Ex.                 effMax=F, colorsubgroup = NULL, ylim = c(-7,7),
#       #Ex.                 outDir = "/home/ddepierre/work/script_repo/plot_function/boxplot/", info = "exemple_2")
#       #Ex. 
#       #Ex. # EXEMPLE 3 : 
#       #Ex. subgroup.lst = list(
#       #Ex.   grp1 = intersect(names(A[A>quantile(A, 0.5)]), names(B[B>quantile(B, 0.5)])),
#       #Ex.   grp2 = names(A[A<quantile(A, 0.1)]),
#       #Ex.   grp3 = paste0("quantif", 1:length(A))[paste0("quantif", 1:length(A)) %ni% c(intersect(names(A[A>quantile(A, 0.5)]), names(B[B>quantile(B, 0.5)])), names(A[A<quantile(A, 0.1)]))]
#       #Ex. )
#       #Ex. 
#       #Ex. boxplot.v220617(data.veclst = data2groups, data.name = "some_quantifications", subgroup.lst = subgroup.lst, subgroup.name = "some_groups",
#       #Ex.                 effMax=50, colorsubgroup = NULL, ylim = c(-7,7),
#       #Ex.                 outDir = "/home/ddepierre/work/script_repo/plot_function/boxplot/", info = "exemple_3_effMax_50")
#       #Ex. 
#       #Ex. boxplot.v220617(data.veclst = data2groups, data.name = "some_quantifications", subgroup.lst = subgroup.lst, subgroup.name = "some_groups",
#       #Ex.                 effMax=T, colorsubgroup = NULL, ylim = c(-7,7),
#       #Ex.                 outDir = "/home/ddepierre/work/script_repo/plot_function/boxplot/", info = "exemple_3_effMax_T")
#       #Ex. 
#       #Ex. # EXEMPLE 4 : 
#       #Ex. C = rnorm(200, -4, 0.5)
#       #Ex. names(C) = paste0("quantif", 1:length(C))
#       #Ex. data3groups = list(A=A, B=B, C=C)
#       #Ex. 
#       #Ex. subgroup.lst = list(
#       #Ex.   grp1 = intersect(names(A[A>quantile(A, 0.5)]), names(B[B>quantile(B, 0.5)])),
#       #Ex.   grp2 = intersect(names(A[A<quantile(A, 0.5)]), names(C[C<quantile(C, 0.5)])),
#       #Ex.   grp3 = paste0("quantif", 1:length(C))[paste0("quantif", 1:length(C)) %ni% c(intersect(names(A[A>quantile(A, 0.5)]), names(B[B>quantile(B, 0.5)])), intersect(names(A[A<quantile(A, 0.5)]), names(C[C<quantile(C, 0.5)])))]
#       #Ex. )
#       #Ex. 
#       #Ex. boxplot.v220617(data.veclst = data3groups, data.name = "some_quantifications", subgroup.lst = subgroup.lst, subgroup.name = "some_groups",
#       #Ex.                 effMax=F, colorsubgroup = NULL, ylim = c(-7,7),
#       #Ex.                 outDir = "/home/ddepierre/work/script_repo/plot_function/boxplot/", info = "exemple_4")
#       #Ex. 
#       #Ex. 

#         # prepare data.lst
#         if(!is.list(data.veclst)){ # if 1 vector, transform in list according to subgroup.lst
#           data.vec = data.veclst
#           if(!is.null(subgroup.lst)){
#             data.lst = lapply(subgroup.lst, function(elem){
#               data.vec[elem]
#             })
#           }else{
#             data.lst = list()
#             data.lst$elem_1 = data.vec
#           }
#         }else if(is.list(data.veclst)){
#           data.lst = data.veclst
#         }
#         # name elements
#         if(is.null(names(data.lst))){ # if no names of list element, create names
#           names(data.lst) = paste0("elem_", 1:length(data.lst))
#         }
#         if(is.null(subgroup.name)){
#           subgroup.name = paste0(names(data.lst), collapse="_")
#         }
#         # prepare colorsubgroup if NULL
#         if(is.null(colorsubgroup) | length(colorsubgroup)<length(data.lst)){
#           colorsubgroup = rainbow(length(data.lst))
#           print((paste0("colorsubgroup created : ",paste0(colorsubgroup, collapse=", "))))
#         }

#     # apply maximum effectif on group / usefull when comparing pvalues to avoid sample-size
#       if(is.null(subgroup.lst)){
#         eff.df = sapply(data.lst, length)
#         rownames(eff.df) = c("eff in groups")
#         if(isTRUE(effMax)){
#           effMax =  min(unlist(lapply(data.lst, length)))
#         }
#         if(is.numeric(effMax)){
#           data.lst.effMax = lapply(data.lst, function(elem){
#             if(length(elem) >= effMax){
#               elem = sample(elem, effMax)
#             }else{
#               elem = elem
#             }
#           })
#           eff.df = rbind(sapply(data.lst, length), sapply(data.lst.effMax, length))
#           eff.df = sapply(data.lst, length)
#           rownames(eff.df) = c("eff in groups", "eff used after effMax")
#           data.lst = data.lst.effMax
#         }
#       }else{
#         eff.df = sapply(subgroup.lst, length)
#         rownames(eff.df) = c("eff in groups")
#         if(isTRUE(effMax)){
#           effMax =  min(unlist(lapply(subgroup.lst, length)))
#         }
#         if(is.numeric(effMax)){
#           subgroup.lst.effMax = lapply(subgroup.lst, function(elem){
#             if(length(elem) >= effMax){
#               elem = sample(elem, effMax)
#             }else{
#               elem = elem
#             }
#           })
#         eff.df = rbind(sapply(subgroup.lst, length), sapply(subgroup.lst.effMax, length))
#         rownames(eff.df) = c("eff in groups", "eff used after effMax")
#         subgroup.lst = subgroup.lst.effMax
#         }
#       }

#     ## create data.frame to be used in ggplot
#       data_ggplot = data.table(unlist(data.lst))
#       rownames(data_ggplot) = names(unlist(data.lst))
#       data_ggplot$elem = rep(names(data.lst), unlist(lapply(data.lst, length)))
#       data_ggplot$elem = factor(data_ggplot$elem, levels = c(names(data.lst)))
#       if(!is.null(subgroup.lst)){
#         data_ggplot$subgroup = NA
#         lapply(1:length(subgroup.lst), function(ndx){
#           data_ggplot$subgroup[unlist(lapply(strsplit(rownames(data_ggplot), "[.]"), "[[", 2)) %in% subgroup.lst[[ndx]]] <<- names(subgroup.lst)[ndx]
#         })
#         data_ggplot$subgroup = factor(data_ggplot$subgroup, levels = c(names(subgroup.lst)))
#         data_ggplot = data_ggplot[!is.na(data_ggplot$subgroup),]
#       }

#     # plotting
#       if(is.null(subgroup.lst) | !is.list(data.veclst)){ # If compare is on subgroup.lst
#             comp_list = combn(names(data.lst),2, simplify=FALSE)

#             boxplot_wilcox = ggplot(data_ggplot, aes(x=elem, y=V1)) + geom_boxplot(outlier.shape = 4, lwd=1,col="#000000", fill = colorsubgroup) + 
#                         stat_compare_means(method = "wilcox.test",comparisons = comp_list) + theme_classic(base_line_size = 1, base_rect_size=1) +
#                         labs(x = subgroup.name, y = data.name) + ggtitle("auto ylim, outliers plotted and ketp for wilcox")

#             boxplot_ttest = ggplot(data_ggplot, aes(x=elem, y=V1)) + geom_boxplot(outlier.shape = 4, lwd=1,col="#000000", fill = colorsubgroup) + 
#                         stat_compare_means(method = "t.test",comparisons = comp_list) + theme_classic(base_line_size = 1, base_rect_size=1) +
#                         labs(x = subgroup.name, y = data.name) + ggtitle("auto ylim, outliers plotted and ketp for ttest")

#             boxplot_NOoutlier = ggplot(data_ggplot, aes(x=elem, y=V1)) + geom_boxplot(outlier.shape = 4, lwd=1,col="#000000", fill = colorsubgroup) + 
#                         scale_y_continuous(limits = quantile(data_ggplot$V1, c(outlier.pct, 1-outlier.pct))) + # remove outlier 
#                         stat_compare_means(method = "wilcox.test",comparisons = comp_list) + theme_classic(base_line_size = 1, base_rect_size=1) +
#                         labs(x = subgroup.name, y = data.name) + ggtitle(paste0("auto ylim ",outlier.pct*100 ,"% outliers removed for plot and wilcox"))

#             boxplot_ylim = ggplot(data_ggplot, aes(x=elem, y=V1)) + geom_boxplot(outlier.shape = 4, lwd=1,col="#000000", fill = colorsubgroup) + 
#                         ylim(ylim) + theme_classic(base_line_size = 1, base_rect_size=1) +
#                         labs(x = subgroup.name, y = data.name) + ggtitle(paste0("figure boxplot"))

#             boxplotviolin_ylim = ggplot(data_ggplot, aes(x=elem, y=V1))+ geom_violin(aes(fill = elem), trim = FALSE) + 
#                         geom_boxplot(outlier.shape = 4, lwd=1,col="#000000", width = 0.2) + scale_fill_manual(values = colorsubgroup) +
#                         ylim(ylim) + theme_classic(base_line_size = 1, base_rect_size=1) +
#                         labs(x = subgroup.name, y = data.name) + ggtitle(paste0("figure violin"))
#       }else if(length(data.lst)==2){
#              comp_list = combn(names(data.lst),2, simplify=FALSE)

#             boxplot_wilcox = ggplot(data_ggplot, aes(x=subgroup, y=V1, fill=elem)) + geom_boxplot(outlier.shape = 4, lwd=1,col="#000000") + 
#                         theme_classic(base_line_size = 1, base_rect_size=1) + scale_fill_manual(values=colorsubgroup) + 
#                         stat_compare_means(method = "wilcox.test", paired=T, method.args = list(alternative = test.side)) +
#                         labs(x = subgroup.name, y = data.name) + ggtitle("auto ylim, outliers plotted and ketp for wilcox")

#             boxplot_ttest = ggplot(data_ggplot, aes(x=subgroup, y=V1, fill=elem)) + geom_boxplot(outlier.shape = 4, lwd=1,col="#000000") + 
#                         theme_classic(base_line_size = 1, base_rect_size=1) + scale_fill_manual(values=colorsubgroup) + 
#                         stat_compare_means(method = "t.test", paired=T, method.args = list(alternative = test.side)) +
#                         labs(x = subgroup.name, y = data.name) + ggtitle("auto ylim, outliers plotted and ketp for ttest")

#             boxplot_NOoutlier = ggplot(data_ggplot, aes(x=subgroup, y=V1, fill=elem)) + geom_boxplot(outlier.shape = 4, lwd=1,col="#000000") + 
#                         scale_y_continuous(limits = quantile(data_ggplot$V1, c(outlier.pct, 1-outlier.pct))) + # remove outlier 
#                         theme_classic(base_line_size = 1, base_rect_size=1) + scale_fill_manual(values=colorsubgroup) + 
#                         labs(x = subgroup.name, y = data.name) + ggtitle(paste0("auto ylim ",outlier.pct*100 ,"% outliers removed for plot and wilcox"))

#             boxplot_ylim = ggplot(data_ggplot, aes(x=subgroup, y=V1, fill=elem)) + geom_boxplot(outlier.shape = 4, lwd=1,col="#000000") + 
#                         ylim(ylim) + theme_classic(base_line_size = 1, base_rect_size=1) + scale_fill_manual(values=colorsubgroup) + 
#                         labs(x = subgroup.name, y = data.name) + ggtitle("figure boxplot")

#             boxplotviolin_ylim = ggplot(data_ggplot, aes(x=subgroup, y=V1, fill=elem)) +
#                         geom_violin(aes(fill = elem), trim = FALSE, lwd=1) + scale_fill_manual(values=colorsubgroup) +
#                         # geom_boxplot(outlier.shape = 4, lwd=1,col="#000000", width = 0.2) + scale_fill_manual(values="white") +
#                         ylim(ylim) + theme_classic(base_line_size = 1, base_rect_size=1)  + 
#                         labs(x = subgroup.name, y = data.name) + ggtitle("figure violin")

#       }else if(length(data.lst)>2){
#             stat.ttest <- data_ggplot %>%
#               group_by(subgroup) %>%
#               pairwise_t_test(
#                   V1 ~ elem, paired = TRUE, 
#                   p.adjust.method = "bonferroni"
#                 ) %>% select(-df, -statistic, -p) # Supprimer les détails
#             stat.ttest <- stat.ttest %>% add_xy_position(x = "subgroup") # add p-values des tests statistiques

#             stat.wilcox <- data_ggplot %>%
#               group_by(subgroup) %>%
#               pairwise_wilcox_test(
#                   V1 ~ elem, paired = TRUE, 
#                   p.adjust.method = "bonferroni"
#                 ) %>% select(-statistic, -p) # Supprimer les détails
#             stat.wilcox <- stat.wilcox %>% add_xy_position(x = "subgroup") # add p-values des tests statistiques


#             boxplot_wilcox = ggboxplot(data_ggplot, x = "subgroup", y = "V1", fill = "elem", palette = colorsubgroup, lwd=1,col="#000000") + 
#                             stat_pvalue_manual(stat.ttest, label = "p.adj", step.increase = 0.08) + 
#                             theme_classic(base_line_size = 1, base_rect_size=1) +
#                             labs(x = subgroup.name, y = data.name) + ggtitle("auto ylim, outliers plotted and ketp for wilcox")

#             boxplot_ttest = ggboxplot(data_ggplot, x = "subgroup", y = "V1", fill = "elem", palette = colorsubgroup, lwd=1,col="#000000") + 
#                             stat_pvalue_manual(stat.ttest, label = "p.adj", step.increase = 0.08) + 
#                             theme_classic(base_line_size = 1, base_rect_size=1) +
#                             labs(x = subgroup.name, y = data.name) + ggtitle("auto ylim, outliers plotted and ketp for t.test")

#             boxplot_NOoutlier = ggboxplot(data_ggplot, x = "subgroup", y = "V1", fill = "elem", palette = colorsubgroup, lwd=1,col="#000000") + 
#                             scale_y_continuous(limits = quantile(data_ggplot$V1, c(outlier.pct, 1-outlier.pct))) + # remove outlier 
#                             theme_classic(base_line_size = 1, base_rect_size=1) +
#                             labs(x = subgroup.name, y = data.name) + ggtitle(paste0("auto ylim ",outlier.pct*100 ,"% outliers removed for plot and wilcox"))

#             boxplot_ylim = ggboxplot(data_ggplot, x = "subgroup", y = "V1", fill = "elem", palette = colorsubgroup, lwd=1,col="#000000") + 
#                             ylim(ylim) +
#                             theme_classic(base_line_size = 1, base_rect_size=1) +
#                             labs(x = subgroup.name, y = data.name) + ggtitle("figure boxplot")

#             boxplotviolin_ylim = ggplot(data_ggplot, aes(x=subgroup, y=V1, fill=elem)) +
#                         geom_violin(aes(fill = elem), trim = FALSE, lwd=1) + scale_fill_manual(values=colorsubgroup) +
#                         # geom_boxplot(outlier.shape = 4, lwd=1,col="#000000", width = 0.2) + scale_fill_manual(values="white") +
#                         ylim(ylim) + theme_classic(base_line_size = 1, base_rect_size=1)  + 
#                         labs(x = subgroup.name, y = data.name) + ggtitle("figure violin")
#       }

#     # plot the plot
#       pdf(paste0(outDir, "BOXPLOT_",data.name ,"_compare_", subgroup.name,"_",info,".pdf"), width = 16, height = 10 )
#             print(boxplot_wilcox)
#             print(boxplot_ttest)
#             print(boxplot_NOoutlier)
#             textplot(eff.df)
#             print(boxplot_ylim)
#             print(boxplotviolin_ylim)
#       dev.off()
# }




































#end