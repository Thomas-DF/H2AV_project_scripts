
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


fisher_namesList=function(list1, list2, total, effMax = F, test.side="greater")
 {
 require(GenomicRanges)
 if(effMax == T){
   effmax_l1 = min(unlist(lapply(list1, length)))
   list1 = lapply(list1, function(vec){
   	sample(vec, effmax_l1)
   })
   effmax_l2 = min(unlist(lapply(list2, length)))
   list2 = lapply(list2, function(vec){
   	sample(vec, effmax_l2)
   })
 }
 m_pval=matrix(nrow=length(list1),ncol=length(list2))
 m_or=matrix(nrow=length(list1),ncol=length(list2))
 for (i in 1:length(list1)){
   liste1=list1[[i]]
    for (j in 1:length(list2)){
      liste2=list2[[j]]
      inter=length(Reduce(intersect, list(liste1,liste2)))
      m_pval[i,j]=signif(fisher.test(matrix(c(inter,length(liste1)-inter,length(liste2)-inter,max(total-length(liste1)-length(liste2)+inter,0)),ncol=2),alternative=test.side)$p.value,2)
      m_or[i,j]=signif(fisher.test(matrix(c(inter,length(liste1)-inter,length(liste2)-inter,max(total-length(liste1)-length(liste2)+inter,0)),ncol=2),alternative=test.side)$estimate,2)
    }
 }
 colnames(m_pval)=names(list2)
 rownames(m_pval)=names(list1)
 colnames(m_or)=names(list2)
 rownames(m_or)=names(list1)
 l=list()
 l[[1]]=m_pval
 l[[2]]=m_or
 names(l)=c("p.value", "odds.ratio")
 return(l)
 }



 Hmap_pval=function(mat, values_ordonnée, values_abscisse, graph,bool, threshold=0.05, col_max =1e-8 ){
   mat_pv=mat$p.value
   mat_pv=mat_pv+10^-300
   mat_pv=log10(mat_pv)
   mat_or = mat$odds.ratio
   if (threshold>1 | threshold<=0) (stop("threshold must be in ]0,1]"))
   logt=log10(threshold*col_max)
   if (bool==T){postscript(graph)}
   # mat$p.value
   ramp <- colorRamp(c("blue", "white"))
   palette(rgb(ramp((200:0)/200),max=255))
   plot.new()
   plot.window(c(-0.4,1.2),c(-0.2,1.2))
   long=ncol(mat_pv)
   long2=nrow(mat_pv)
   for (j in 1:ncol(mat_pv)){
     for (k in 1:nrow(mat_pv)){
      if(mat$p.value[k,j]<=threshold){
        rect((j-1)/long,(k-1)/long2,j/long,k/long2,col=((mat_pv[k,j])*(mat_pv[k,j]>logt)*(200/logt))+200*(mat_pv[k,j]<logt))
        
      }
      if(mat$p.value[k,j]>threshold){
        rect((j-1)/long,(k-1)/long2,j/long,k/long2,col="white")


      }
        text(-0.13,(k-1/2)/long2,rownames(mat_pv)[k],col='black',cex=0.7)
      
     }
     text((j-1/2)/long,1.13,colnames(mat_pv)[j],col='black',srt=45,cex=0.7)
   }
   palette(rgb(ramp((200:0)/200),max=255))
   plot.new()
   plot.window(c(-0.4,1.2),c(-0.2,1.2))
   for (i in 0:200){
     rect((200-i)/200,0,(199-i)/200,1,col=(i),border=i)
   }
   par(mfrow=c(2,1))
   pval_info = matrix(c(1, threshold*col_max))
   rownames(pval_info) = c("white", "blue")
   textplot(pval_info)
   textplot(mat$p.value)
   par(mfrow=c(1,1))
   # mat$odds.ratio
   # ThOdds = 1
   mat_or = (log(mat_or)+1)*100
   mat_or[mat_or > 200] = 200
   mat_or[mat_or < 0] = 0
   ramp <- colorRamp(c("darkblue", "white", "yellow"))
   palette(rgb(ramp((200:0)/200),max=255))
   plot.new()
   plot.window(c(-0.4,1.2),c(-0.2,1.2))
   long=ncol(mat_or)
   long2=nrow(mat_or)
   for (j in 1:ncol(mat_or)){
     for (k in 1:nrow(mat_or)){
       rect((j-1)/long,(k-1)/long2,j/long,k/long2,col=((mat_or[k,j])))
       # rect((j-1)/long,(k-1)/long2,j/long,k/long2,col=((mat_or[k,j])*(mat_or[k,j]>ThOdds)*(200/ThOdds))+200*(mat_or[k,j]<ThOdds))
       text(-0.13,(k-1/2)/long2,rownames(mat_or)[k],col='black',cex=0.7)
     }
     text((j-1/2)/long,1.13,colnames(mat_or)[j],col='black',srt=45,cex=0.7)
   }
   palette(rgb(ramp((200:0)/200),max=255))
   plot.new()
   plot.window(c(-0.4,1.2),c(-0.2,1.2))
   for (i in 0:200){
     rect((200-i)/200,0,(199-i)/200,1,col=(i),border=i)
    }
    textplot(mat$odds.ratio)
    # textplot(lapply(values_ordonnée, length))
    # textplot(lapply(values_abscisse, length))
    if(is.character(values_ordonnée[[1]])){ #if the split of genes nmaes is not done on quantitatif value, plot effectif of each feature
      barplot(unlist(lapply(values_ordonnée, length)), cex.names=0.5)
    }else{
      boxplot(values_ordonnée, main=deparse(substitute(values_ordonnée)), outline=F)
      barplot(unlist(values_ordonnée), main=deparse(substitute(values_ordonnée)))
      plot(density(unlist(values_ordonnée)), main=deparse(substitute(values_ordonnée)))
      hist(unlist(values_ordonnée), 50, main=deparse(substitute(values_ordonnée)))
    }
    if(is.character(values_abscisse[[1]])){ #if the split of genes nmaes is not done on quantitatif value, plot effectif of each feature
      barplot(unlist(lapply(values_abscisse, length)), cex.names=0.5)
    }else{
      boxplot(values_abscisse, main=deparse(substitute(values_abscisse)), outline=F)
      barplot(unlist(values_abscisse), main=deparse(substitute(values_abscisse)))
      plot(density(unlist(values_abscisse)), main=deparse(substitute(values_abscisse)))
      hist(unlist(values_abscisse), 50, main=deparse(substitute(values_abscisse)))
    }
   if (bool==T){
     dev.off()
   }
 }






# end
