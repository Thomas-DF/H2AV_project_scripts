
# Single Cell RNA seq project
# Cuvier''s Team
# Schaak - Heurteau* - Depierre
# 2017

#~ Function to bin profil matrix
#~ USAGE :

#~ matrix.binned = t(apply(matrix,1,binMean,xTend=ncol(matrix),binSize=10))

#~ Note that ncol(matrix) has to be multiple of binSize




binMean <- function(vec,xTend,binSize){
  cut <- xTend/binSize
  vec <- tapply(vec, createBins(vec,cut), mean)
  return(vec)
}


binMedian <- function(vec,xTend,binSize){
  cut <- xTend/binSize
  vec <- tapply(vec, createBins(vec,cut), median)
  return(vec)
}

createBins <- function(x, bin.count) {
  bin.size <- rep(length(x) %/% bin.count, bin.count)
  bin.size <- bin.size + ifelse(1:bin.count <= length(x) %% bin.count, 1, 0)
  bin <- rep(1:bin.count, bin.size)
  return(factor(bin, levels=1:bin.count))
}
