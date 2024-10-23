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
  F_tilde<-sqrt(t)*vectors
  P_tilde<-(1/t)*t(F_tilde)%*%X
  
  # Compute PC hat
  H_tilde<-(1/N)*P_tilde%*%t(P_tilde)
  F_hat<-(1/N)*X%*%t(P_tilde)
  
  
  # compute residuals
  Residuals <- X - F_tilde %*% P_tilde
  
  Factors_list <- list()
  Factors_hat_list <- list()
  Loadings_list <- list()
  Residuals_list <- list()
  
  Factors_list[[1]] <- F_tilde
  Factors_hat_list[[1]] <- F_hat
  Loadings_list[[1]] <- P_tilde
  Residuals_list[[1]] <- Residuals
  
  
  return(list(Factors = Factors_list, Factors_hat = Factors_hat_list,  
              Lambda = P_tilde, Loadings = Loadings_list, Residuals = Residuals_list))
}



