# Multi-Level Dynamic Factor Model - 3 blocks full


# beta_ols<-function(X,Y){
#   t(Y) %*% X %*% solve(t(X) %*% X)
# }


beta_ols <- function(X, Y) {
  solve(t(X) %*% X) %*% t(X) %*% Y
}

MutiLevelFactor3Blocks<-function(Yorig,r,block_ind,tol,max_iter){

  # Standardize the original data
  Yorig <- scale(Yorig)
  
  # Initialize variables
  results <- list()  # List to save results
  num_blocks <- 3 # Number of blocks
  num_vars <- ncol(Yorig) # Total number of variables in Yorig
  num_obs <- nrow(Yorig) # Total number of observations
  num_factors <- 2^num_blocks-1 # Total number of factors
  
  # Define block ranges based on input indices
  range1 <- 1:block_ind[1]
  range2 <- (block_ind[1] + 1):block_ind[2]
  range3 <- (block_ind[2] + 1):block_ind[3]
  
  
  # Get number of variables in each block
  num_vars_block1 <- length(range1)
  num_vars_block2 <- length(range2)
  num_vars_block3 <- length(range3)
  
  
  ### STEP 1 ###
  
  # a)
  # Canonical Correlation Analysis (CCA)
  GlobalFactors <- blockfact0(Yorig, c(num_vars_block1, num_vars_block2, num_vars_block3), 1, c(1, 1, 1))
  GlobalFactors <- GlobalFactors / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(GlobalFactors) %*% GlobalFactors))))
  
  # b)
  # Remove global component to get pairwise residuals
  
  Block12 <- cbind(Yorig[,range1],Yorig[,range2])
  Block13 <- cbind(Yorig[,range1],Yorig[,range3])
  Block23 <- cbind(Yorig[,range2],Yorig[,range3])
  
  Beta12 <-  beta_ols(GlobalFactors, Block12)
  Beta13 <-  beta_ols(GlobalFactors, Block13)
  Beta23 <-  beta_ols(GlobalFactors, Block23)
  
  Residuals12 <- Block12 - GlobalFactors %*% Beta12
  Residuals13 <- Block13 - GlobalFactors %*% Beta13
  Residuals23 <- Block23 - GlobalFactors %*% Beta23
  
  # c)
  # Apply CCA for lower-level block factors
  Factors12 <- blockfact0(Residuals12, c(num_vars_block1, num_vars_block2), 1, c(1, 1))
  Factors13 <- blockfact0(Residuals13, c(num_vars_block1, num_vars_block3), 1, c(1, 1))
  Factors23 <- blockfact0(Residuals23, c(num_vars_block2, num_vars_block3), 1, c(1, 1))
  
  # d)
  # Calculate residuals for each block using the lower-level factors
  
  Block1 <- Yorig[, range1]
  Block2 <- Yorig[, range2]
  Block3 <- Yorig[, range3]
  
  Beta1 <- beta_ols(cbind(Factors12, Factors13), Block1)
  Beta2 <- beta_ols(cbind(Factors12, Factors23), Block2)
  Beta3 <- beta_ols(cbind(Factors13, Factors23), Block3)
  
  Residuals1 <- Block1 - cbind(Factors12, Factors13) %*% Beta1
  Residuals2 <- Block2 - cbind(Factors12, Factors23) %*% Beta2
  Residuals3 <- Block3 - cbind(Factors13, Factors23) %*% Beta3
 
  # f)
  # Estimate initial factors for each single block using eigenvalue decomposition (PCA)
  number_of_factor = 1
  
  eig_vectors1 <- eigrs2(t(Residuals1) %*% Residuals1)$evec
  Factors1 <-  Residuals1 %*% eig_vectors1[, (ncol(eig_vectors1) - number_of_factor + 1):ncol(eig_vectors1)]
  Factors1 <- Factors1 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors1) %*% Factors1))))
  
  
  eig_vectors2 <- eigrs2(t(Residuals2) %*% Residuals2)$evec
  Factors2 <-  Residuals2 %*% eig_vectors2[, (ncol(eig_vectors2) - number_of_factor + 1):ncol(eig_vectors2)]
  Factors2 <- Factors2 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors2) %*% Factors2))))
  
  eig_vectors3 <- eigrs2(t(Residuals3) %*% Residuals3)$evec
  Factors3 <-  Residuals3 %*% eig_vectors3[, (ncol(eig_vectors3) - number_of_factor + 1):ncol(eig_vectors3)]
  Factors3 <- Factors3 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors3) %*% Factors3))))
  
  # Store all initial factors
  Factors <- cbind(GlobalFactors, Factors12, Factors13, Factors23, Factors1, Factors2, Factors3)
  InitialFactors <- Factors
  
  # g)
  # Compute initial loadings
  
  #Lambda <- beta_ols(InitialFactors, Yorig)
  
  GlobalLoadings <- beta_ols(GlobalFactors, Yorig)
  Loadings12 <- beta_ols(Factors12, Block12)
  Loadings13 <- beta_ols(Factors13, Block13)
  Loadings23 <- beta_ols(Factors23, Block23)
  Loadings1 <- beta_ols(Factors1, Block1)
  Loadings2 <- beta_ols(Factors2, Block2)
  Loadings3 <- beta_ols(Factors3, Block3)
  
  Lambda <- matrix(0, nrow = num_factors, ncol = num_vars)
  Lambda[1,] <- GlobalLoadings
  Lambda[2,c(range1,range2)] <- Loadings12
  Lambda[3,c(range1,range3)] <- Loadings13
  Lambda[4,c(range2,range3)] <- Loadings23
  Lambda[5,range1] <- Loadings1
  Lambda[6,range2] <- Loadings2
  Lambda[7,range3] <- Loadings3
  
  InitialLoadings <- Lambda


  
  ### STEP 2 ###
  
  # Initialize residual sum of squares (RSS) for convergence
  iteration <- 0
  Residuals <- Yorig - Factors %*% Lambda
  RSS_previous <- sum(diag(t(Residuals) %*% Residuals))
  

  
  # Iterative procedure for convergence
  while (iteration < max_iter) {
    
    
    
    iteration <- iteration + 1
    
    #cat("Iteration:", iteration, "\n")
    
    # Run least squares to estimate global factors 
    GlobalFactors <- t(beta_ols(t(GlobalLoadings), t(Yorig)))
    
    
    # Regress to filter out the global component
    Beta12 <-  beta_ols(GlobalFactors, Block12)
    Beta13 <-  beta_ols(GlobalFactors, Block13)
    Beta23 <-  beta_ols(GlobalFactors, Block23)
    
    Residuals12 <- Block12 - GlobalFactors %*% Beta12
    Residuals13 <- Block13 - GlobalFactors %*% Beta13
    Residuals23 <- Block23 - GlobalFactors %*% Beta23
    
    # Run least squares to estimate pairwise factors
    Factors12 <- t(beta_ols(t(Loadings12), t(Residuals12)))
    Factors13 <- t(beta_ols(t(Loadings13), t(Residuals13)))
    Factors23 <- t(beta_ols(t(Loadings23), t(Residuals23)))
    
    # Regress to filter out the pairwise component
    Beta1 <-  beta_ols(cbind(Factors12,Factors13), Block1)
    Beta2 <-  beta_ols(cbind(Factors12,Factors23), Block2)
    Beta3 <-  beta_ols(cbind(Factors13,Factors23), Block3)

    Residuals1 <- Block1 - cbind(Factors12, Factors13) %*% Beta1
    Residuals2 <- Block2 - cbind(Factors12, Factors23) %*% Beta2
    Residuals3 <- Block3 - cbind(Factors13, Factors23) %*% Beta3
    
    
    # Run least squares to estimate single block factors
    Factors1 <- t(beta_ols(t(Loadings1), t(Residuals1)))
    Factors2 <- t(beta_ols(t(Loadings2), t(Residuals2)))
    Factors3 <- t(beta_ols(t(Loadings3), t(Residuals3)))
    
    # Update and normalize factors 
    GlobalFactors <- GlobalFactors / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(GlobalFactors) %*% GlobalFactors))))

    Factors12 <- Factors12 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors12) %*% Factors12))))
    Factors13 <- Factors13 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors13) %*% Factors13))))
    Factors23 <- Factors23 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors23) %*% Factors23))))

    Factors1 <- Factors1 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors1) %*% Factors1))))
    Factors2 <- Factors2 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors2) %*% Factors2))))
    Factors3 <- Factors3 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors3) %*% Factors3))))
    
    Factors <- cbind(GlobalFactors, Factors12, Factors13, Factors23, Factors1, Factors2, Factors3)
    
    # Compute and update loadings
    #Lambda <- beta_ols(Factors, Yorig)
    
    GlobalLoadings <- beta_ols(GlobalFactors, Yorig)
    Loadings12 <- beta_ols(Factors12, Block12)
    Loadings13 <- beta_ols(Factors13, Block13)
    Loadings23 <- beta_ols(Factors23, Block23)
    Loadings1 <- beta_ols(Factors1, Block1)
    Loadings2 <- beta_ols(Factors2, Block2)
    Loadings3 <- beta_ols(Factors3, Block3)
    
    Lambda[1,] <- GlobalLoadings
    Lambda[2,c(range1,range2)] <- Loadings12
    Lambda[3,c(range1,range3)] <- Loadings13
    Lambda[4,c(range2,range3)] <- Loadings23
    Lambda[5,range1] <- Loadings1
    Lambda[6,range2] <- Loadings2
    Lambda[7,range3] <- Loadings3
    

    # Check RSS
    Residuals <- Yorig - Factors %*% Lambda
    RSS_new <- sum(diag(t(Residuals) %*% Residuals))
    
  
    #print(RSS_previous)
    #print(RSS_new)
    

    if (abs(RSS_previous - RSS_new) < tol) {
      break  # Converged
    }
    
    RSS_previous <- RSS_new
    
   
   
  }
  
  
  
  results[["Factors"]] <- Factors
  results[["Loadings"]] <- Lambda
  
  
  return(results)
  
  
}


