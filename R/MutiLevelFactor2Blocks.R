MutiLevelFactors2Blocks <- function(data, r, block_ind = block_ind) {
  # Based on Djogbenou, A. A. (2023)
  # Last modified: January 2024
  # It gives us loadings and factors estimated from a 2-level factor panel
  # model using the sequential principal component
  
  # Outputs:
  # lambda: estimated loadings
  # factors: estimated factors
  # Sigma_ee: estimated variance for idiosyncratic errors
  # tol: precision
  # count: number of iterations
  
  # Inputs:
  # X: t*n matrix,
  # n1: Number of developed economies real activity variables
  # r0: Number of global factor,
  # r01,r02: Number of specific factors.
  
  
  r0 <- r[1]
  r01 <- r[2]
  r02 <- r[3]
  
  
  X <- as.matrix(data)
  
  t <- nrow(X)
  #n <- ncol(X)
  n <- block_ind[2] # total number of var
  
  n1 <-  block_ind[1] # number of var in block 1
  n2 <- block_ind[2] - block_ind[1] # number of var in block 2
 
  X <- X[,1:block_ind[2]] # global block
  X1 <- X[, 1:block_ind[1]] # block 1
  X2 <- X[, (block_ind[1] + 1):block_ind[2]] # block 2
  
  
  
  # First step Factor estimation
  # First step: Initial estimates
  #Lambda_0<-eigen(X %*% t(X) / (n * t))$values # lambda
  #Factor_0<-eigen(X %*% t(X) / (n * t))$vectors # factor
  #round(Factor_0[,1]-svd_result$u[,1],5)
  #round(Lambda_0-svd_result$d,5)
  svd_result <- svd(X %*% t(X) / (n * t))
  vectors <- svd_result$u
  values <- svd_result$d
  v <- svd_result$v
  Ghat <- sqrt(t) * vectors[, 1:r0]
  

  
  # Initialization
  count <- 0
  tol <- 1000
  critere <- 1000
  Lhat <- t(X) %*% Ghat %*% solve(t(Ghat) %*% Ghat) # BETA OLS
  if(r0>1) Y <- X - Ghat[, 1:r0] %*% t(Lhat) # ERROR
  if(r0==1) Y <- X - Ghat %*% t(Lhat) # ERROR
  
  while (tol > 0.000001) {
    y1 <- Y[, 1:n1]
    y2 <- Y[, (1 + n1):n]
    
    
    # Step 2
    svd_result1 <- svd(y1 %*% t(y1) / (n1 * t))
    Vectors1 <- svd_result1$u
    Values1 <- svd_result1$d
    V1 <- svd_result1$v
    
    Ghat1 <- sqrt(t) * Vectors1[, 1:r01]
    Lhat1 <- t(y1) %*% Ghat1 %*% solve(t(Ghat1) %*% Ghat1) #t(y1) %*% Ghat1 / t
    if(r01>1) Ehat1 <- X1 - Ghat1[,1:r01] %*% t(Lhat1)
    if(r01==1) Ehat1 <- X1 - Ghat1 %*% t(Lhat1) #Ghat1 %*% Lhat1
    
    
    svd_result2 <- svd(y2 %*% t(y2) / (n2 * t))
    Vectors2 <- svd_result2$u
    Values2 <- svd_result2$d
    V2 <- svd_result2$v
    
    Ghat2 <- sqrt(t) * Vectors2[, 1:r02]
    Lhat2 <- t(y2) %*% Ghat2 / t
    if(r02>1) Ehat2 <- X2 - Ghat2[,1:r02] %*% t(Lhat2)
    if(r02==1) Ehat2 <- X2 - Ghat2 %*% t(Lhat2) #X2 - Ghat2 %*% Lhat2
    
    
    # Step 3
    Ehat <- cbind(Ehat1, Ehat2)
    svd_result3 <- svd(Ehat %*% t(Ehat) / (n * t))
    Vectors <- svd_result3$u
    Values <- svd_result3$d
    V <- svd_result3$v
    
    Ghat <- sqrt(t) * Vectors[, 1:r0] #Vectors
    Lhat <- t(Ehat) %*% Ghat / t
    Y <- X - Ghat %*% t(Lhat) #X - Ghat[, 1:r0] %*% t(Lhat)
    
    
    critere0 <- mean((Ehat - Ghat %*% t(Lhat))^2)
    count <- count + 1
    tol <- abs(critere - critere0)
    critere <- critere0
    cat("Iter:", count,"\n")
    Sigma_ee <- diag(mean((Ehat - Ghat %*% t(Lhat))^2),1)
    if(r02==2) lambda <- cbind(Lhat, rbind(Lhat1, matrix(0, n2, 1)),rbind(matrix(0, n1, 2), as.matrix(Lhat2)))
    if(r02==1) lambda <- cbind(Lhat, rbind(Lhat1, matrix(0, n2, 1)), rbind(matrix(0, n1, 1), Lhat2))
    factors <- cbind(Ghat, cbind(Ghat1, Ghat2))
  }
  
  return(list(Factors = factors, Loadings = NULL))
  
  return(list(lambda = lambda, factors = factors, Sigma_ee = Sigma_ee, tol = tol, count = count))
}