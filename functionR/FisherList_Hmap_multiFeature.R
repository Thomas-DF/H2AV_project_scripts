
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

# vec_gn = topUpK27 ; list2 = List_MultiFeature_splitted_gn ; total = comRownames_5994

# fisher_namesList_MultiFeature=function(vec_gn,list2,total)
#  {
#    if(is.numeric(total)%in%FALSE){total = length(total)}
#    list1 = rep(list(vec_gn), length(list2))
#    sublist2 = list2[[1]]
#  require(GenomicRanges)
#  m_pval=matrix(nrow=length(list1),ncol=length(sublist2))
#  m_or=matrix(nrow=length(list1),ncol=length(sublist2))
#  for (i in 1:length(list1)){
#    liste1=list1[[i]]
#     for (j in 1:length(list2)){
#       sublist2=list2[[j]]
#       # liste2=list2[[j]]
#       for(k in 1:length(sublist2)){
#               inter=length(Reduce(intersect, list(liste1,sublist2)))
#               m_pval[i,k]=signif(fisher.test(matrix(c(inter,length(liste1)-inter,length(sublist2)-inter,max(total-length(liste1)-length(sublist2)+inter,0)),ncol=2),alternative="greater")$p.value,2)
#               m_or[i,k]=signif(fisher.test(matrix(c(inter,length(liste1)-inter,length(sublist2)-inter,max(total-length(liste1)-length(sublist2)+inter,0)),ncol=2),alternative="greater")$estimate,2)
#       }
#     }
#  }
#  colnames(m_pval)=names(sublist2)
#  rownames(m_pval)=names(list1)
#  colnames(m_or)=names(sublist2)
#  rownames(m_or)=names(list1)
#  l=list()
#  l[[1]]=m_pval
#  l[[2]]=m_or
#  names(l)=c("p.value", "odds.ratio")
#  return(l)
#  }


fisher_namesList_MultiFeature=function(vec_gn,list2,total)
 {
   if(is.numeric(total)%in%FALSE){total = length(total)}
   liste1 = vec_gn
   sublist = list2[[1]]
 require(GenomicRanges)
 m_pval=matrix(nrow=length(list2),ncol=length(sublist))
 m_or=matrix(nrow=length(list2),ncol=length(sublist))
    for (j in 1:length(list2)){
      sublist2=list2[[j]] #sublist2=list2[[2]]
      # liste2=list2[[j]]
      for(k in 1:length(sublist2)){
              liste2 = sublist2[[k]]
              inter=length(Reduce(intersect, list(liste1,liste2)))
              m_pval[j,k]=signif(fisher.test(matrix(c(inter,length(liste1)-inter,length(liste2)-inter,max(total-length(liste1)-length(liste2)+inter,0)),ncol=2),alternative="greater")$p.value,2)
              m_or[j,k]=signif(fisher.test(matrix(c(inter,length(liste1)-inter,length(liste2)-inter,max(total-length(liste1)-length(liste2)+inter,0)),ncol=2),alternative="greater")$estimate,2)
      }
    }
 colnames(m_pval)=names(sublist)
 rownames(m_pval)=names(list2)
 colnames(m_or)=names(sublist)
 rownames(m_or)=names(list2)
 l=list()
 l[[1]]=m_pval
 l[[2]]=m_or
 names(l)=c("p.value", "odds.ratio")
 return(l)
 }



 Hmap_pval=function(mat, values_ordonnée=NULL, values_abscisse=NULL, graph,bool, threshold=1e-8){
   mat_pv=mat$p.value
   mat_pv=mat_pv+10^-300
   mat_pv=log10(mat_pv)
   mat_or = mat$odds.ratio
   if (threshold>1 | threshold<=0) (stop("threshold must be in ]0,1]"))
   logt=log10(threshold)
   if (bool==T){postscript(graph)}
   # mat$p.value
   ramp <- colorRamp(c("blue", "white"))
   palette(rgb(ramp((200:0)/200),max=255))
   plot.new()
   plot.window(c(-0.4,1.2),c(0,1))
   long=ncol(mat_pv)
   long2=nrow(mat_pv)
   for (j in 1:ncol(mat_pv)){
     for (k in 1:nrow(mat_pv)){
       rect((j-1)/long,(k-1)/long2,j/long,k/long2,col=((mat_pv[k,j])*(mat_pv[k,j]>logt)*(200/logt))+200*(mat_pv[k,j]<logt))
       text(-0.13,(k-1/2)/long2,rownames(mat_pv)[k],col='black',cex=0.7)
     }
     text((j-1/2)/long,1.13,colnames(mat_pv)[j],col='black',srt=45,cex=0.7)
   }
   # mat$odds.ratio
   # ThOdds = 1
   mat_or = (log(mat_or)+1)*100
   mat_or[mat_or > 200] = 200
   mat_or[mat_or < 0] = 0
   ramp <- colorRamp(c("darkblue", "white", "yellow"))
   palette(rgb(ramp((200:0)/200),max=255))
   plot.new()
   plot.window(c(-0.4,1.2),c(0,1))
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
   ramp <- colorRamp(c("blue", "white"))
   palette(rgb(ramp((200:0)/200),max=255))
   plot.new()
   plot.window(c(-0.4,1.2),c(0,1))
   for (i in 0:200){
     rect((200-i)/200,0,(199-i)/200,1,col=(i),border=i)
   }
   ramp <- colorRamp(c("darkblue", "white", "yellow"))
   palette(rgb(ramp((200:0)/200),max=255))
   plot.new()
   plot.window(c(-0.4,1.2),c(0,1))
   for (i in 0:200){
     rect((200-i)/200,0,(199-i)/200,1,col=(i),border=i)
    }
    # boxplot(values_ordonnée, main=deparse(substitute(values_ordonnée)), outline=F)
    # barplot(unlist(values_ordonnée), main=deparse(substitute(values_ordonnée)))
    # plot(density(unlist(values_ordonnée)), main=deparse(substitute(values_ordonnée)))
    # hist(unlist(values_ordonnée), 50, main=deparse(substitute(values_ordonnée)))
    # boxplot(values_abscisse, main=deparse(substitute(values_abscisse)), outline=F)
    # barplot(unlist(values_abscisse), main=deparse(substitute(values_abscisse)))
    # plot(density(unlist(values_abscisse)), main=deparse(substitute(values_abscisse)))
    # hist(unlist(values_abscisse), 50, main=deparse(substitute(values_abscisse)))
   if (bool==T){
     dev.off()
   }
 }






# end
