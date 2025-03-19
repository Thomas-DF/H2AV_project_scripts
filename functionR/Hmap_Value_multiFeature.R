
# Single Cell RNA seq project
# Cuvier''s Team
# Schaak - Heurteau - Depierre*
# 2017

# modified function from pascal/gael
# fisher_namesList :
# fisher function adapted to test intersection genes list against genes list
# list 1 and 2 in input are list of char vector, total correspond to the total number of genes, statistical "univers"

# Hmap_pval :
# plot pval from fisher_namesList
# input :
# mat = output return from fisher_namesList
# graph = "/path/to/output/mygraph.pdf"
# bool = T to exprt pdf grap according to graph input, if False, plot directed in X11
# threshold=1e-8 correspond to pval threshold for the beginning of the blue colorRamp

Get_mtxMeanValue_MultiFeature=function(mtx_Value,list2)
 {
 sublist = list2[[1]]
 require(GenomicRanges)
 m_MeanValue=matrix(nrow=length(list2),ncol=length(sublist))
 m_VARValue=matrix(nrow=length(list2),ncol=length(sublist))
    for (j in 1:length(list2)){
      sublist2=list2[[j]] #sublist2=list2[[2]]
      # liste2=list2[[j]]
      for(k in 1:length(sublist2)){
              liste2 = sublist2[[k]]
              m_MeanValue[j,k] = mean(mtx_Value[rownames(mtx_Value) %in% liste2,1])
              m_VARValue[j,k] = var(mtx_Value[rownames(mtx_Value) %in% liste2,1])
      }
    }
 colnames(m_MeanValue)=names(sublist)
 rownames(m_MeanValue)=names(list2)
 colnames(m_VARValue)=names(sublist)
 rownames(m_VARValue)=names(list2)
 l=list()
 l[[1]]=m_MeanValue
 l[[2]]=m_VARValue
 names(l)=c("mean_value", "VAR_value")
 return(l)
 }



 Hmap_Value = function(mat_ToPlot, title, graph,bool=T){
   if (bool==T){postscript(graph)}
   mat_ToPlot_mean = mat_ToPlot$mean_value
   # mat_or = mat_test$mean_value
   if(min(mat_ToPlot_mean)<0){
     mat_cut <- cut(as.vector(mat_ToPlot_mean), breaks = seq(min(as.vector(mat_ToPlot_mean)), max(as.vector(mat_ToPlot_mean)), len = 1000),include.lowest = TRUE)
     # mat_cut <- cut(as.vector(mat_ToPlot), breaks = seq(min(as.vector(mat_ToPlot)), max(as.vector(mat_ToPlot)), len = 200),include.lowest = TRUE)
     vec_ToPlot_ramp = as.vector(rep(NA, length(mat_cut)))
     for(i in 1:length(vec_ToPlot_ramp)){vec_ToPlot_ramp[i] = which(levels(mat_cut) %in% mat_cut[i])}
     mat_ToPlot_ramp = matrix(vec_ToPlot_ramp, ncol=10)

     ramp <- colorRamp(c("darkblue", "white", "darkred"))
     palette(rgb(ramp((1000:0)/1000),max=255))
     # palette(rgb(ramp((200:0)/200),max=255))
     plot.new()
     plot.window(c(-0.4,1.2),c(0,1))
     long=ncol(mat_ToPlot_mean)
     long2=nrow(mat_ToPlot_mean)
     for (j in 1:ncol(mat_ToPlot_mean)){
       for (k in 1:nrow(mat_ToPlot_mean)){
         rect((j-1)/long,(k-1)/long2,j/long,k/long2,col=((mat_ToPlot_ramp[k,j])))
         text(-0.13,(k-1/2)/long2,rownames(mat_ToPlot_mean)[k],col='black',cex=0.7)
       }
       text((j-1/2)/long,1.13,colnames(mat_ToPlot_mean)[j],col='black',srt=45,cex=0.7)
     }
     ramp <- colorRamp(c("darkblue", "white", "darkred"))
     palette(rgb(ramp((200:0)/200),max=255))
     plot.new()
     plot.window(c(-0.4,1.2),c(0,1))
     title(main=title, " min: ", min(mat_ToPlot_mean), " max: ", max(mat_ToPlot_mean))
     for (i in 0:200){
       rect((200-i)/200,0,(199-i)/200,1,col=(i),border=i)
     }
   }else{
     mat_cut <- cut(as.vector(mat_ToPlot_mean), breaks = seq(min(as.vector(mat_ToPlot_mean)), max(as.vector(mat_ToPlot_mean)), len = 1000),include.lowest = TRUE)
     # mat_cut <- cut(as.vector(mat_ToPlot), breaks = seq(min(as.vector(mat_ToPlot)), max(as.vector(mat_ToPlot)), len = 200),include.lowest = TRUE)
     vec_ToPlot_ramp = as.vector(rep(NA, length(mat_cut)))
     for(i in 1:length(vec_ToPlot_ramp)){vec_ToPlot_ramp[i] = which(levels(mat_cut) %in% mat_cut[i])}
     mat_ToPlot_ramp = matrix(vec_ToPlot_ramp, ncol=10)

     ramp <- colorRamp(c("darkblue", "white"))
     palette(rgb(ramp((1000:0)/1000),max=255))
     # palette(rgb(ramp((200:0)/200),max=255))
     plot.new()
     plot.window(c(-0.4,1.2),c(0,1))
     long=ncol(mat_ToPlot_mean)
     long2=nrow(mat_ToPlot_mean)
     for (j in 1:ncol(mat_ToPlot_mean)){
       for (k in 1:nrow(mat_ToPlot_mean)){
         rect((j-1)/long,(k-1)/long2,j/long,k/long2,col=((mat_ToPlot_ramp[k,j])))
         text(-0.13,(k-1/2)/long2,rownames(mat_ToPlot_mean)[k],col='black',cex=0.7)
       }
       text((j-1/2)/long,1.13,colnames(mat_ToPlot_mean)[j],col='black',srt=45,cex=0.7)
     }
     ramp <- colorRamp(c("darkblue", "white"))
     palette(rgb(ramp((200:0)/200),max=255))
     plot.new()
     plot.window(c(-0.4,1.2),c(0,1))
     title(main=title, " min: ", min(mat_ToPlot_mean), " max: ", max(mat_ToPlot_mean))
     for (i in 0:200){
       rect((200-i)/200,0,(199-i)/200,1,col=(i),border=i)
     }
   }
   mat_ToPlot_VAR = mat_ToPlot$VAR_value
   # mat_or = mat_test$mean_value
   if(min(mat_ToPlot_VAR)<0){
     mat_cut <- cut(as.vector(mat_ToPlot_VAR), breaks = seq(min(as.vector(mat_ToPlot_VAR)), max(as.vector(mat_ToPlot_VAR)), len = 1000),include.lowest = TRUE)
     # mat_cut <- cut(as.vector(mat_ToPlot), breaks = seq(min(as.vector(mat_ToPlot)), max(as.vector(mat_ToPlot)), len = 200),include.lowest = TRUE)
     vec_ToPlot_ramp = as.vector(rep(NA, length(mat_cut)))
     for(i in 1:length(vec_ToPlot_ramp)){vec_ToPlot_ramp[i] = which(levels(mat_cut) %in% mat_cut[i])}
     mat_ToPlot_ramp = matrix(vec_ToPlot_ramp, ncol=10)

     ramp <- colorRamp(c("darkblue", "white", "darkred"))
     palette(rgb(ramp((1000:0)/1000),max=255))
     # palette(rgb(ramp((200:0)/200),max=255))
     plot.new()
     plot.window(c(-0.4,1.2),c(0,1))
     long=ncol(mat_ToPlot_VAR)
     long2=nrow(mat_ToPlot_VAR)
     for (j in 1:ncol(mat_ToPlot_VAR)){
       for (k in 1:nrow(mat_ToPlot_VAR)){
         rect((j-1)/long,(k-1)/long2,j/long,k/long2,col=((mat_ToPlot_ramp[k,j])))
         text(-0.13,(k-1/2)/long2,rownames(mat_ToPlot_VAR)[k],col='black',cex=0.7)
       }
       text((j-1/2)/long,1.13,colnames(mat_ToPlot_VAR)[j],col='black',srt=45,cex=0.7)
     }
     ramp <- colorRamp(c("darkblue", "white", "darkred"))
     palette(rgb(ramp((200:0)/200),max=255))
     plot.new()
     plot.window(c(-0.4,1.2),c(0,1))
     title(main=title, " min: ", min(mat_ToPlot_VAR), " max: ", max(mat_ToPlot_VAR))
     for (i in 0:200){
       rect((200-i)/200,0,(199-i)/200,1,col=(i),border=i)
     }
   }else{
     mat_cut <- cut(as.vector(mat_ToPlot_VAR), breaks = seq(min(as.vector(mat_ToPlot_VAR)), max(as.vector(mat_ToPlot_VAR)), len = 1000),include.lowest = TRUE)
     # mat_cut <- cut(as.vector(mat_ToPlot), breaks = seq(min(as.vector(mat_ToPlot)), max(as.vector(mat_ToPlot)), len = 200),include.lowest = TRUE)
     vec_ToPlot_ramp = as.vector(rep(NA, length(mat_cut)))
     for(i in 1:length(vec_ToPlot_ramp)){vec_ToPlot_ramp[i] = which(levels(mat_cut) %in% mat_cut[i])}
     mat_ToPlot_ramp = matrix(vec_ToPlot_ramp, ncol=10)

     ramp <- colorRamp(c("black", "white"))
     palette(rgb(ramp((1000:0)/1000),max=255))
     # palette(rgb(ramp((200:0)/200),max=255))
     plot.new()
     plot.window(c(-0.4,1.2),c(0,1))
     long=ncol(mat_ToPlot_VAR)
     long2=nrow(mat_ToPlot_VAR)
     for (j in 1:ncol(mat_ToPlot_VAR)){
       for (k in 1:nrow(mat_ToPlot_VAR)){
         rect((j-1)/long,(k-1)/long2,j/long,k/long2,col=((mat_ToPlot_ramp[k,j])))
         text(-0.13,(k-1/2)/long2,rownames(mat_ToPlot_VAR)[k],col='black',cex=0.7)
       }
       text((j-1/2)/long,1.13,colnames(mat_ToPlot_VAR)[j],col='black',srt=45,cex=0.7)
     }
     ramp <- colorRamp(c("black", "white"))
     palette(rgb(ramp((200:0)/200),max=255))
     plot.new()
     plot.window(c(-0.4,1.2),c(0,1))
     title(main=title, " min: ", min(mat_ToPlot_VAR), " max: ", max(mat_ToPlot_VAR))
     for (i in 0:200){
       rect((200-i)/200,0,(199-i)/200,1,col=(i),border=i)
     }
   }


   if (bool==T){dev.off()}
 }









 # boxplot(values_ordonnée, main=deparse(substitute(values_ordonnée)), outline=F)
 # barplot(unlist(values_ordonnée), main=deparse(substitute(values_ordonnée)))
 # plot(density(unlist(values_ordonnée)), main=deparse(substitute(values_ordonnée)))
 # hist(unlist(values_ordonnée), 50, main=deparse(substitute(values_ordonnée)))
 # boxplot(values_abscisse, main=deparse(substitute(values_abscisse)), outline=F)
 # barplot(unlist(values_abscisse), main=deparse(substitute(values_abscisse)))
 # plot(density(unlist(values_abscisse)), main=deparse(substitute(values_abscisse)))
 # hist(unlist(values_abscisse), 50, main=deparse(substitute(values_abscisse)))



# end
