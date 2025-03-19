
#~ Single Cell RNA seq project
#~ Cuvier s Team
#~ Schaak - Heurteau - Depierre* 
#~ 2017


## plot fisher function
# dataToPlot.DF must have feat, name, odds and pval columns as:
#~  matFisher
#~       feat feat_nbr    name name_nbr inter  odds   pval P.Value
#~ 1  mes4_UR      566 mes4_UR      566   566   Inf    Inf    0.00
#~ 2  mes4_UR      566 mes4_DR      566     0  0.00   0.00    1.00
#~ 3  mes4_UR      566 hypb_UR      586   215 21.38 157.01    0.00
#~ 4  mes4_UR      566 hypb_DR      586    37  1.63   2.26    0.01
#~ 5  mes4_DR      566 mes4_UR      566     0  0.00   0.00    1.00
#~ 6  mes4_DR      566 mes4_DR      566   566   Inf    Inf    0.00
#~ 7  mes4_DR      566 hypb_UR      586    41  1.83   3.29    0.00
#~ 8  mes4_DR      566 hypb_DR      586   168 13.03  99.46    0.00
#~ 9  hypb_UR      586 mes4_UR      566   215 21.38 157.01    0.00
#~ 10 hypb_UR      586 mes4_DR      566    41  1.83   3.29    0.00
#~ 11 hypb_UR      586 hypb_UR      586   586   Inf    Inf    0.00
#~ 12 hypb_UR      586 hypb_DR      586     0  0.00   0.00    1.00
#~ 13 hypb_DR      586 mes4_UR      566    37  1.63   2.26    0.01
#~ 14 hypb_DR      586 mes4_DR      566   168 13.03  99.46    0.00
#~ 15 hypb_DR      586 hypb_UR      586     0  0.00   0.00    1.00
#~ 16 hypb_DR      586 hypb_DR      586   586   Inf    Inf    0.00

require(scales)
require(ggplot2)
require(gplots)
plotFisher =function(dataToPlot.DF, o_name){
	# plot Pvalue fisher
  p <- ggplot(data =  dataToPlot.DF, aes(x = name, y = feat)) +
    geom_tile(aes(fill=pval,width=0.97, height=0.97)) +
    geom_text(aes(fill=pval,label = as.character(round(pval,2))), size = 2) +
    with(dataToPlot.DF,scale_fill_gradientn(limits = c(0,max(pval[which(pval<Inf)])),
     values = rescale(c(0,-log10(0.05),max(pval[which(pval<Inf)]))),
    colours=c("grey", "white", "red"),
    breaks=c(0,-log10(0.05),max(pval[which(pval<Inf)])), labels=format(c(0,-log10(0.05),max(pval[which(pval<Inf)])))))+
    coord_fixed(ratio=0.7) +
    theme(
      # AXIS MODIFICATION
        axis.text.x = element_text(vjust=0.5,size=4,angle = 0,colour="black"),
        axis.text.y = element_text(size=4,colour="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
     # PANEL MODIFICATION
        panel.background = element_rect(fill = "white",colour="black",size=0.7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
     # LEGEND MODIFICATION
     legend.title = element_blank(),
     legend.key = element_rect(fill="white"), legend.background = element_rect(fill=NA))

	# plot odds ratio fisher
  q <- ggplot(data =  dataToPlot.DF, aes(x = name, y = feat)) +
    geom_tile(aes(fill=odds,width=0.97, height=0.97)) +
    geom_text(aes(fill=odds,label = as.character(round(odds,2))), size = 2) +
    with(dataToPlot.DF,scale_fill_gradientn(limits = c(0,max(odds[which(odds<Inf)])),
     values = rescale(c(0,1,max(odds[which(odds<Inf)]))),
    colours=c("#4A9060", "white", "red"),
    breaks=c(0,1,max(odds[which(odds<Inf)])), labels=format(c(0,1,max(odds[which(odds<Inf)])))))+
    coord_fixed(ratio=0.7) +
    theme(
      # AXIS MODIFICATION
        axis.text.x = element_text(vjust=0.5,size=4,angle = 0,colour="black"),
        axis.text.y = element_text(size=4,colour="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
     # PANEL MODIFICATION
        panel.background = element_rect(fill = "white",colour="black",size=0.7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
     # LEGEND MODIFICATION
     legend.title = element_blank(),
     legend.key = element_rect(fill="white"), legend.background = element_rect(fill=NA))

pdf(o_name)
print(p)
print(q)
dev.off()
}


