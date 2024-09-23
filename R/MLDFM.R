# Multi-Level Dynamic Factor Model

MLDFM <- function(data, outlier = TRUE, r = c(1), blocks = 1, block_ind = NULL, tol = 0.000001, max_iter = 1000) {
  
  ##### OUTLIERS #####
  outliers <- NULL
  # Handle outliers
  if (outlier) {
    OutlierResult <- HandleOutliers(data,r)
    data <- OutlierResult$data
    outliers <- OutlierResult$outliers
  }
  
  
  ##### FACTOR EXTRACTION #####
  
  if (blocks==1){
    result <- SingleBlock(data,r=r)
  }else if(blocks>1){
    #result <- MultipleBlocks(data, r=r,block_ind = block_ind, tol = tol, max_iter = max_iter)
    result <- MultipleBlocks3(data, r=r,block_ind = block_ind, tol = tol, max_iter = max_iter)
    
  }else{
    print('Error - Invalid number of block')
  }
  
  return(list(Data = data, Outliers = outliers, Factors = result$Factors, Loadings = result$Loadings ))
}


# Compute and correct outliers
HandleOutliers <- function(data, r ) {
  
  # Scale data
  X<-scale(data,TRUE,TRUE)
  
  # Compute dimensions
  t<-nrow(X)
  N<-ncol(X)
  
  # Number of obs to compute the median value
  N_median = 5
  
  # Compute first r eigenvalues and eigenvectors
  eR<-eigen(X%*%t(X))
  values<-eR$values[c(1:r)]
  vectors<-matrix(eR$vectors[,c(1:r)],t,r)
  
  # Compute PC and loadings
  F_hat<-sqrt(t)*vectors
  P_hat<-(1/t)*t(F_hat)%*%X
  
  # Compute idiosyncratic components
  Idio_comp<-X-F_hat%*%P_hat
  
  # Search and correct outliers
  outliers <- matrix(0,t,N)
  
  for (i in 1:N) {
    max <- quantile(Idio_comp[,i],0.75) + (IQR(Idio_comp[,i]) * 6)
    min <- quantile(Idio_comp[,i],0.25) - (IQR(Idio_comp[,i]) * 6)
    for (ii in (N_median+1):t) {
      if (Idio_comp[ii,i] < min | Idio_comp[ii,i] > max) {
        outliers[ii,i]<-1 
        data[ii,i]<-median(data[c((ii-N_median):(ii-1)),i])} 
    }
  }
  
  return(list(data = data, outliers = outliers))
}
  
