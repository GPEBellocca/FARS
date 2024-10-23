# Multi-Level Dynamic Factor Model - Single block

SingleBlock <- function(data, r) {
  
  ##### PRINCIPAL COMPONENTS #####
  
  # Scale data
  X<-scale(data,TRUE,TRUE)
  
  # Compute dimensions
  t<-nrow(X)
  N<-ncol(X)
  
  # Compute first r eigenvalues and eigenvectors
  eR<-eigen(X%*%t(X))
  values<-eR$values[c(1:r)]
  vectors<-matrix(eR$vectors[,c(1:r)],t,r)
  
  # Compute PC and loadings
  F_hat<-sqrt(t)*vectors
  P_hat<-(1/t)*t(F_hat)%*%X
  
  # compute residuals
  Residuals <- X - F_hat %*% P_hat
  
  return(list(Factors = F_hat, Lamba = P_hat, Loadings = P_hat, Residuals = Residuals))
}



