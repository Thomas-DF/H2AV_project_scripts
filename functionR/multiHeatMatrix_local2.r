
multiHeatMatrix_local2<-function(sml,grid=TRUE,col=NULL,xcoords=NULL,
                          group=NULL,group.col=NULL,KmeansUsed=FALSE,
                          order=FALSE,user.order=FALSE,
                          winsorize=c(0,100),
                          clustfun=FALSE,RangeToCluster = NULL,
                          column.scale=TRUE,rangesList=NULL,
                          matrix.main=NULL,
                          common.scale=FALSE,legend=TRUE,
                          legend.name=NULL,cex.legend=0.8,
                          xlab=NULL,
                          cex.lab=1, cex.main=1,cex.axis=0.8,newpage=TRUE)
{
.jets<-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                          "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))  
#   if( class(sml) !="ScoreMatrixList" ){stop("'sml' is not a ScoreMatrix object\n")}
  
  # check if col and xcoords are lists
  # if so their number should be equal to the number of matrices
  if(is.list(col) & length(col) != length(sml)){
    stop("'col' is a list and its length does not match the length of ScoreMatrixList\n")
    col=NULL
  }
  
  if( is.list(xcoords) & length(xcoords) != length(sml) ){
    stop("'xcoords' is a list and its length does not match the length of ScoreMatrixList\n")
    xcoords=NULL
  }
  
  # check legend.name argument if , it is a vector
  if(!is.null(legend.name) & length(legend.name)>1 & length(legend.name) != length(sml)){
    stop("'legend.name' should match the length of the 'sml' ",
         "if it is not a length 1 vector\n")
  }else if (length(legend.name)==1){
    legend.name=rep(legend.name,length(sml))
  }
  
  # check xlab argument and set it
  if(!is.null(xlab) & length(xlab)>1 & length(xlab) != length(sml)){
    stop("'xlab' should match the length of the 'sml' ",
         "if it is not a length 1 vector\n")
  }else if (length(xlab)==1){
    xlab=rep(xlab,length(sml))
  } 
  
  #check number of rows if they match or not
  # if not, warn and run unionScoreMatrixList
  if( length(unique(sapply(sml,nrow)))>1  ){
    warning("\nThe row numbers are different\n",
            "attempting to get common rows and to reorder rows\n",
            "using 'intersectScoreMatrixList()'\n")
    sml=intersectScoreMatrixList(sml,reorder=TRUE)
  }
  
  
  
  mat.list=lapply(sml,function(x) x) # get the matrix class, some operations are not transitive
  
  
  # if this is changed, a rowSide color map will be drawn
  # setting clustfun or group.list will populate this vector
  # and this function will check on its value later on
  # to decide to plot rowside colors or not
  group.vector=NULL
  group.names=NULL # will be plotted as annotation if filled in later
  
  
  # this can set extreme values to given percentile
  # better for visualization
  # alternative is to take log or sth
  # but that should be done prior to heatMatrix
  if(winsorize[2]<100 | winsorize[1]>0){
    mat.list=lapply(mat.list,function(x) .winsorize(x,winsorize) )
  }
  
  
  # if order is true | clustfun is provided
  # make a one large matrix by cbind
  if(!identical(clustfun, FALSE) | order){
    mat2=do.call("cbind",mat.list)
    if(column.scale){
      mat2=scale(mat2)
      mat2[is.nan(mat2)]=0
    }
  }
  
  # do clustfun if requested
  if(!identical(clustfun, FALSE)){
    if(!is.null(RangeToCluster)){
        mat_forClustering.list=lapply(sml,function(x) x[,RangeToCluster])
        if(winsorize[2]<100 | winsorize[1]>0){
            mat_forClustering.list=lapply(mat_forClustering.list,function(x) .winsorize(x,winsorize) )
        }
        if(!identical(clustfun, FALSE) | order){
            mat2_clu=do.call("cbind",mat_forClustering.list)
            if(column.scale){
                mat2_clu=scale(mat2_clu)
                mat2_clu[is.nan(mat2_clu)]=0
            }
        }
        if(any(is.na(mat2_clu)) ) {
            mat3_clu=impute.knn(mat2_clu ,k = 10, 
                            rowmax = 0.5, colmax = 0.8, 
                            maxp = 1500)$data
            if(KmeansUsed){
                clusterToReturn = clustfun(mat3_clu)
                clu = clusterToReturn$cluster
            }else {
               clu=clustfun(mat3_clu)
            }
        }
        else{
            if(KmeansUsed){
                clusterToReturn = clustfun(mat2_clu)
                clu = clusterToReturn$cluster
            }else {
               clu=clustfun(mat2_clu)
            }
            # cluster 
            # clu=clustfun(mat2_clu)
        }
    }else{
        # impute values if there are NAs
        if(any(is.na(mat2)) ) {
            mat3=impute.knn(mat2 ,k = 10, 
                            rowmax = 0.5, colmax = 0.8, 
                            maxp = 1500)$data
            if(KmeansUsed){
                clusterToReturn = clustfun(mat3)
                clu = clusterToReturn$cluster
            }else {
               clu=clustfun(mat3)
            }
        }
        else{
            if(KmeansUsed){
                clusterToReturn = clustfun(mat2)
                clu = clusterToReturn$cluster
            }else {
               clu=clustfun(mat2)
            }
        }
    }
    
    
    # get group.vector centers, will be used at ordering later
    group.vector=clu
    
    # order things by clusters only
    mat.list=lapply(mat.list,function(x) x[order(group.vector),])
    group.vector=group.vector[order(group.vector)]
    
    # if user wants to order
    if(order){
      
      # replicate center value for each cluster
      g.factor=factor(group.vector,levels=unique(group.vector))
      
      # do ordering within clusters
      my.order=order(group.vector,-rowSums(mat2,na.rm=TRUE))
      
      # commence the new order: Novus Ordo Seclorum
      mat.list=lapply(mat.list,function(x) x[my.order,])
      group.vector=group.vector[my.order]
    }
    
    group.names=unique(group.vector) # to be used for the rowSide colors
  }
  
  
  # check conditions of group.list
  # group must not have duplicated numbers
  # warn if total number group elements is below nrow(mat)
  if(!is.null(group)  & identical(clustfun, FALSE)){
    
    # if group is a list of rowids, row numbers for original windows argument
    if(is.list(group)){
      # group must not have duplicated numbers
      win.numbs=(lapply(group, function(x) unique(x) ))
      win.vec=unlist(win.numbs)
      if(any(table(unlist(win.numbs))>1)){
        stop("'group' containing a list must not have duplicated numbers\n")
      }
      row.ids=rownames(mat.list[[1]]) # get row ids from the matrix
      group.vector=rep(0,length(row.ids)) # make a starting vector of zeros
      for(i in 1:length(win.numbs)){
        group.vector[row.ids %in% win.numbs[[i]] ]=i # make group vector
      }
      
      # if list has names, take it as group.names
      # group.names will be plotted with rowSide colors
      if(!is.null(names(group))){
        group.names=names(group)
      }
      
      # stop if total number group elements is below nrow(mat)
      if(all(group.vector==0)){
        stop("None of the elements in 'group' are a part of rownames(mat) \n")
      }
      
      # warn if total number group elements is below nrow(mat)
      if(any(group.vector==0)){
        warning("Number of elements in 'group' argument is less then nrow(mat) \n",
                "Dropping rows from 'mat' that are not contained in 'group'\n")
        
        mat.list=lapply(mat.list,function(x) x[group.vector>0,])
        group.vector=group.vector[group.vector>0]
      }       
      
      
    }
    else if(is.factor(group)){
      
      # check if the length is same as nrow(mat)
      if(length(group) != nrow(mat.list[[1]])){
        stop("'group' is a factor, and its length should be equal to nrow(mat)\n")
      }
      # arrange factor levels by the order of appeareance
      #group=factor(as.character(group),levels=as.character(unique(group)))
      group.names=levels(group)
      levels(group)=1:length(levels(group))
      group.vector=as.numeric(group)
    }else{
      stop("'group' must be a factor or a list\n")
    }
    
    
    
    # order things by clusters only
    mat.list=lapply(mat.list,function(x) x[order(group.vector),])
    group.vector=group.vector[order(group.vector)]
    
    
    if(order){
      
      # get cbound matrix for ordering
      mat2=do.call("cbind",mat.list)
      if(column.scale){
        mat2=scale(mat2)
        mat2[is.nan(mat2)]=0
      }
      my.order=order(group.vector,-rowSums(mat2,na.rm=TRUE))
      
      # commence the new order: Novus Ordo Seclorum
      mat.list=lapply(mat.list,function(x) x[my.order,])
      group.vector=group.vector[my.order]
    }
    
  }else if(order & identical(clustfun, FALSE) ){ # if only ordering is needed no group or clustering
    
    # get cbound matrix for ordering
    mat2=do.call("cbind",mat.list)
    if(column.scale){
      mat2=scale(mat2)
      mat2[is.nan(mat2)]=0
    }
    order.vector       =rep(1,nrow(mat2))
    names(order.vector)=rownames(mat2)
    my.order=order(-rowSums(mat2,na.rm=TRUE))
    mat.list=lapply(mat.list,function(x) x[my.order,])
    order.vector       = order.vector[my.order]
    
  }
  
  #if user.order is provided
  if(!identical(user.order, FALSE) & !is.null(group.names)){
    if(length(user.order)!=length(group.names)){
      warning(paste0("length of 'user.order' vector (", length(user.order) ,") should be the the same as number of clusters (",length(group.names),"). Skipping it..\n"))
    }else{
        print(group.vector)
        gv.fac.ord <- order(factor(group.vector, levels = user.order))
    #   gv.fac.ord <- sort(group.vector)
      group.vector <- group.vector[gv.fac.ord] #convert factor to vector
      print(group.vector)
      group.names <- user.order
      print(group.names)
    }
  }else if(!identical(user.order, FALSE) & is.null(group.names)){
    warning("There are no groups or clusters to order. Skipping it..")
  }
  
  
  
  # THE PLOTTING STARTS HERE with ordered mat.list
  if(!grid){
    plot.new()
    vps <- baseViewports()
    pushViewport(vps$figure) # get the plot viewport from base graphics
    
  }else{
    if(newpage)
      grid.newpage()
  }
  
  # try to get matrix names from names of ScoreMatrixList
  if(is.null(matrix.main) & !is.null(names(sml)) & class(sml)=="ScoreMatrixList" )
  {
    matrix.main=  names(sml)
  }else if(!is.null(matrix.main) & length(matrix.main) != length(sml)){
    warning("'matrix.main' length does not match to the 'sml' length\n",
            "setting it to NULL")
    matrix.main=NULL
    
  }
  
  # calculate the width of heatmaps and everything else relative to that
  l.sml=length(mat.list) # get length of matrix
  # get the npc width of the heatmap
  # this will be our rudder when plotting
  # every coordinate and dimension will be relative to this value
  hw=40*(1-(convertX( unit(4,"lines"),"npc",valueOnly=TRUE)) )/(l.sml*44+7)
  hh=0.7 # this heatmap height in NPC
  hy=0.6 # heatmap y coordinate on the plot
  
  # if group.vector is not NULL, plot the rowside colors
  # make side colors
  if(!is.null(group.vector)){
    sideVp <- viewport(width=unit(hw/8, "npc"), height=unit(hh, "npc"),
                       x = convertX( unit(4,"lines"),"npc"), 
                       y = unit(hy, "npc"),just="left")
    pushViewport(sideVp) # push the viewport
    
    grid.rect()
    
    # make the rowside colors and labels
    .rowSideCol(group.vector,group.names=group.names,
                group.col=group.col,cex.lab=cex.lab)
    
    popViewport()
  }
  
  # heatmap start coordinates
  heat.startCoord=convertX( unit(4,"lines"),"npc",valueOnly=TRUE) + (hw/8) + (hw/20)
  
  # if the same scale to be used for all plots
  common.range=NULL
  if(common.scale){
    common.range=range(mat.list,na.rm=TRUE)
  }else{
    if(!is.null(rangesList)&length(rangesList)!=length(mat.list)){warning("length of ranges values list should be equal to the length of matrices to plot! Going back to auto scale or common.scale mode!")}
  }
  
  
  for(i in 1:length(mat.list)){
    
    
    # cxcoords stands for current xcoords
    if(is.list(xcoords)){
      cxcoords=xcoords[[i]]
    }else if( is.vector(xcoords) | is.null(xcoords) ){
      cxcoords=xcoords   
    }
    else if(!is.null(xcoords)){
      warning("xcoords should be a vector or a list or NULL,nothing else!!\n",
              "setting xcoords to NULL")
      cxcoords=NULL
    }
    
    # if it is a two element vector of start and end coordinates 
    # of the first and the last column of the matrix
    if(length(cxcoords)==2){
      cxcoords=seq(cxcoords[1],cxcoords[2],length.out=ncol(mat.list[[i]]))
    }
    
    # ccol: stands for current.color
    if(is.list(col)){
      ccol=col[[i]]
    }else if(is.vector(col) | is.null(col) ){
      ccol=col
    }else if(!is.null(col)){
      warning("col should be a vector or a list or NULL,nothing else!!\n",
              "setting colors to default")
      ccol=NULL
    }
    
    
    # get the default/given xcoordinates to plot
    if(!is.null(cxcoords)){
      if(length(cxcoords) != ncol( mat.list[[i]] ) ) 
        stop("xcoords has wrong length: ",length(cxcoords)," \n",
             " it should be equal to the number of columns of ScoreMatrix\n",
             " which is",ncol(mat.list[[i]]),"\n")    
    }else{
      cxcoords=1:ncol( mat.list[[i]] )
    }
    
    # get heatcolor scale
    if(is.null(ccol)) 
      ccol=.jets(100)
    
    
    # make heatmap viewport
    heatVp <- viewport(width=unit(hw, "npc"), height=unit(hh, "npc"),
                       x =unit(heat.startCoord,"npc"), 
                       y = unit(hy, "npc"),just="left")
    
    
    pushViewport(heatVp) # push the heatmap VP
    grid.rect()
    if(!is.null(rangesList) & length(rangesList)==length(mat.list)){
        mat_temp = mat.list[[i]]
        mat_temp[mat_temp > rangesList[[i]][2]] = rangesList[[i]][2]
        mat_temp[mat_temp < rangesList[[i]][1]] = rangesList[[i]][1]
        .gridHeat(mat_temp,ccol,xcoords=cxcoords,xlab[i],cex.lab,cex.axis,
              angle=60,hjust=0.6,vjust=1,rng=rangesList[[i]]) # make the heat  
    }else {
       .gridHeat(mat.list[[i]],ccol,xcoords=cxcoords,xlab[i],cex.lab,cex.axis,
              angle=60,hjust=0.6,vjust=1,rng=common.range) # make the heat  
    }
    
    #upViewport(1) # up one level current.vpTree()
    popViewport()
    
    if(legend){
      # make the legend viewport
      legendVp <- viewport(width=unit(hw*0.7, "npc"), 
                           height=unit(0.5, "lines"),
                           x =unit(heat.startCoord+hw*0.15,"npc"), 
                           y = unit(0.1, "npc"),
                           just="left")
      pushViewport(legendVp)
      grid.rect()
      
      if(common.scale){
        rng=common.range
      }else{
        if(!is.null(rangesList) & length(rangesList)==length(mat.list)){
            rng = rangesList[[i]]
        }else {
           rng=range(mat.list[[i]],na.rm=TRUE)
        }        
      }
      
      # make the legend 
      .heatLegendX(min=rng[1],max=rng[2],
                   ccol,legend.name[i],main=TRUE,cex.legend,
                   cex.lab,vjust=-0.5,hjust=0.5)
      popViewport()
    }
    
    if(!is.null(matrix.main)){
      title.y=unit(0.96, "npc")
      grid.text(matrix.main[i], y=title.y, 
                x = unit(heat.startCoord+hw*0.5,"npc"),
                gp=gpar(cex=cex.main),just="bottom")
    }
    
    # update start coord for next one
    heat.startCoord=heat.startCoord+ hw+ (hw/20)
    
  }
  
  if(!grid){
    popViewport()
  }
  
  if(!identical(clustfun, FALSE) | !is.null(group.vector)){
    if(KmeansUsed){
        return(invisible(clusterToReturn)) 
    }else {
       return(invisible(group.vector)) 
    }    
  }else if(order & is.null(group.vector) ){
    if(KmeansUsed){
        return(invisible(clusterToReturn)) 
    }else {
       return(invisible(order.vector)) 
    }
  }
  
}