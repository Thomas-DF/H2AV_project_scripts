#
# Stéphane Schaak - nov2016
# David Dépierre - mai 2017
#
# plot2profilesNew=function(dataToPlot,xRange=-2000:2000,yLims=NULL,dataNames=names(dataToPlot),colors=c("darkblue","red","darkgrey"),outputFile=NULL,fontSizes=list(axisNames=1,axisGrads=0.8,legend=0.8)) {
#   if (!is.null(outputFile)) {
#     pdf(file=outputFile,paper="a4")
#   }
#   if (is.null(yLims)) {
#     yMax=ceiling(max(c(dataToPlot[[1]],dataToPlot[[2]]))/5)*5
#     yMin=floor(min(c(dataToPlot[[1]],dataToPlot[[2]]),0)/5)*5
#     yLims=c(yMin,yMax)
#   }
#   if (is.null(dataNames)) {
#     dataNames=c("sample01","sample02")
#   }
#   plot(smooth.spline(dataToPlot[[1]],spar=0,x=xRange),col=colors[[1]],type='l',ylim=yLims,xlab="Distance from TSS (bp)",ylab="",cex.lab=fontSizes$axisNames,cex.axis=fontSizes$axisGrads)
#   abline(v=0,lty=5,col=colors[[length(colors)]])
#   lines(smooth.spline(dataToPlot[[2]],spar=0,x=xRange),col=colors[[2]],type='l')
#   legend('topright',legend=dataNames,lty=1,col=colors[1:(length(colors)-1)],text.col=colors[1:(length(colors)-1)],cex=fontSizes$legend)
#   if (!is.null(outputFile)) {
#     dev.off()
#   }
# }
# #
#
#
# plot1profileNew=function(dataToPlot,xRange=-2000:2000,yLims=NULL,dataNames=names(dataToPlot),colors=c("red"),outputFile=NULL,fontSizes=list(axisNames=1,axisGrads=0.8,legend=0.8)) {
#   if(length(xRange) != length(dataToPlot)){
# 	  dataToPlot = dataToPlot[(((length(dataToPlot)+1)/2)+xRange[1]):(((length(dataToPlot)+1)/2)+xRange[length(xRange)])]
#   }
#   if (!is.null(outputFile)) {
#     pdf(file=outputFile,12,6)
#   }
#   if (is.null(yLims)) {
#     yMax=ceiling(max(dataToPlot)/5)*5
#     yMin=floor(min(dataToPlot,0)/5)*5
#     yLims=c(yMin,yMax)
#   }
#   if (is.null(dataNames)) {
#     dataNames="sample01"
#   }
#   plot(smooth.spline(dataToPlot,spar=0,x=xRange),col=colors,type='l',ylim=yLims,xlab="Distance from TSS (bp)",ylab="",cex.lab=fontSizes$axisNames,cex.axis=fontSizes$axisGrads)
#   abline(v=0,lty=5,col=colors[[length(colors)]])
# #~   lines(smooth.spline(dataToPlot[[2]],spar=0,x=xRange),col=colors[[2]],type='l')
#   legend('topright',legend=dataNames,lty=1,col=colors,text.col=colors,cex=fontSizes$legend)
#   if (!is.null(outputFile)) {
#     dev.off()
#   }
# }
#
# plot10profilesNew=function(dataToPlot,xRange=-2000:2000,yLims=NULL,dataNames=names(dataToPlot),col=c("darkblue","red4","darkgrey"),outputFile=NULL,fontSizes=list(axisNames=1,axisGrads=0.8,legend=0.8)) {
#   if (!is.null(outputFile)) {
#     pdf(file=outputFile,paper="a4")
#   }
#   if (is.null(yLims)) {
#     yMax=ceiling(max(c(dataToPlot[[1]],dataToPlot[[10]]))/5)*5
#     yMin=floor(min(c(dataToPlot[[1]],dataToPlot[[10]]),0)/5)*5
#     yLims=c(yMin,yMax)
#   }
#   if (is.null(dataNames)) {
#     dataNames=names(dataToPlot)
# #~     dataNames=c("sample01","sample02", "sample03", "sample4", "sample5", "sample6", "sample7", "sample8", "sample9", "sample10")
#   }
#   colors = colorRampPalette(c(col[1], col[2]))(10)
#   plot(smooth.spline(dataToPlot[[1]],spar=0,x=xRange),col=colors[[1]],type='l',ylim=yLims,xlab="Distance from TSS (bp)",ylab="",cex.lab=fontSizes$axisNames,cex.axis=fontSizes$axisGrads)
#   abline(v=0,lty=5,col=col[[length(col)]])
#   lines(smooth.spline(dataToPlot[[2]],spar=0,x=xRange),col=colors[[2]],type='l')
#   lines(smooth.spline(dataToPlot[[3]],spar=0,x=xRange),col=colors[[3]],type='l')
#   lines(smooth.spline(dataToPlot[[4]],spar=0,x=xRange),col=colors[[4]],type='l')
#   lines(smooth.spline(dataToPlot[[5]],spar=0,x=xRange),col=colors[[5]],type='l')
#   lines(smooth.spline(dataToPlot[[6]],spar=0,x=xRange),col=colors[[6]],type='l')
#   lines(smooth.spline(dataToPlot[[7]],spar=0,x=xRange),col=colors[[7]],type='l')
#   lines(smooth.spline(dataToPlot[[8]],spar=0,x=xRange),col=colors[[8]],type='l')
#   lines(smooth.spline(dataToPlot[[9]],spar=0,x=xRange),col=colors[[9]],type='l')
#   lines(smooth.spline(dataToPlot[[10]],spar=0,x=xRange),col=colors[[10]],type='l')
#   legend('topright',legend=dataNames,lty=1,col=colors[1:(length(colors))],text.col=colors[1:(length(colors))],cex=fontSizes$legend)
#   if (!is.null(outputFile)) {
#     dev.off()
#   }
# }
#



print("USAGE :  avg_prof_plot_chipseq(LIST_of_AVG_prof, smoothed = 0.6, interval = c(1,10001), intervalx = c(-5000,5000), main_title = 'average plot of RNApol2 TSS(-5000+5000) genes profiles',0,22,'plot_RNApol2ALL_avg_prof_TSS_-5000+5000') ")


avg_prof_plot_chipseq=function(data_to_plot,smoothed,interval, intervalx,main_title, mlim,plim,output){
col_vect=c('black','red','red','orange','purple','blue','turquoise4')
pdf(paste(dirprof,output,'.pdf',sep=''))
##########################################################################################################
# WITH LEGEND
if (smoothed>0){
	plot(smooth.spline(data_to_plot[[1]][interval[1]:interval[2]],spar=smoothed,x=seq(intervalx[1],intervalx[2],1)),main=main_title,type='l',col=col_vect[1],ylim=c(mlim,plim),xlab='',ylab='')
    if(length(data_to_plot)>1){
      for (i in 2:length(data_to_plot)){
        lines(smooth.spline(data_to_plot[[i]][interval[1]:interval[2]],spar=smoothed,x=seq(intervalx[1],intervalx[2],1)),col=col_vect[i])
      }
      legend('topright',legend=names(data_to_plot),text.col=c('black','red','red','orange','purple','blue','turquoise4'))
    }
	}else{
	plot(seq(intervalx[1],intervalx[2],1),data_to_plot[[1]][interval[1]:interval[2]],main=main_title,type='l',col=col_vect[1],ylim=c(mlim,plim),xlab='',ylab='')
  if(length(data_to_plot)>1){
    for (i in 2:length(data_to_plot)){
      lines(seq(intervalx[1],intervalx[2],1),data_to_plot[[i]][interval[1]:interval[2]],col=col_vect[i])
		}
  }
	legend('topright',legend=names(data_to_plot),text.col=c('black','red','red','orange','purple','blue','turquoise4'))
	}
##########################################################################################################
# WITHOUT  LEGEND
if (smoothed>0){
	plot(smooth.spline(data_to_plot[[1]][interval[1]:interval[2]],spar=smoothed,x=seq(intervalx[1],intervalx[2],1)),main=main_title,type='l',col=col_vect[1],ylim=c(mlim,plim),xlab='',ylab='')
    if(length(data_to_plot)>1){
      for (i in 2:length(data_to_plot)){
        lines(smooth.spline(data_to_plot[[i]][interval[1]:interval[2]],spar=smoothed,x=seq(intervalx[1],intervalx[2],1)),col=col_vect[i])
      }
    }
	}else{
	plot(seq(intervalx[1],intervalx[2],1),data_to_plot[[1]][interval[1]:interval[2]],main=main_title,type='l',col=col_vect[1],ylim=c(mlim,plim),xlab='',ylab='')
  if(length(data_to_plot)>1){
    for (i in 2:length(data_to_plot)){
      lines(seq(intervalx[1],intervalx[2],1),data_to_plot[[i]][interval[1]:interval[2]],col=col_vect[i])
		}
  }
	}
dev.off()
}
