# Compute and correct outliers
HandleOutliers <- function(data, r ) {
  
  X<-scale(data,TRUE,TRUE)
  
  
  t<-nrow(X)
  N<-ncol(X)
 
  N_median = 5
  
  eR<-eigen(X%*%t(X))
  values<-eR$values[c(1:r)]
  vectors<-matrix(eR$vectors[,c(1:r)],t,r)
  
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

