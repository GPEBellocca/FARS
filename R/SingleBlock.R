# Multi-Level Dynamic Factor Model - Single block

SingleBlock <- function(data, r) {
  
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
  F_tilde<-sqrt(t)*vectors
  P_tilde<-(1/t)*t(F_tilde)%*%X
  
  # Compute PC hat
  H_tilde<-(1/N)*P_tilde%*%t(P_tilde)
  F_hat<-(1/N)*X%*%t(P_tilde)
  
  
  # compute residuals
  Residuals <- X - F_tilde %*% P_tilde
  
  Final_list <- list()
  Final_list['1'] <- r
  
  
  return(list(Factors = F_tilde, Lambda = P_tilde, Residuals = Residuals , Factors_hat = F_hat, Factors_list= Final_list))
}



