# Multi-Level Dynamic Factor Model - 3 blocks full


# beta_ols<-function(X,Y){
#   t(Y) %*% X %*% solve(t(X) %*% X)
# }


beta_ols <- function(X, Y) {
  solve(t(X) %*% X) %*% t(X) %*% Y
}

MutiLevelFactor4Blocks<-function(Yorig,r,block_ind,tol,max_iter){

  # Standardize the original data
  Yorig <- scale(Yorig)
  
  # Initialize variables
  results <- list()  # List to save results
  num_blocks <- 4 # Number of blocks
  num_vars <- ncol(Yorig) # Total number of variables in Yorig
  num_obs <- nrow(Yorig) # Total number of observations
  num_factors <- 2^num_blocks-1 # Total number of factors
  
  # Define block ranges based on input indices
  range1 <- 1:block_ind[1]
  range2 <- (block_ind[1] + 1):block_ind[2]
  range3 <- (block_ind[2] + 1):block_ind[3]
  range4 <- (block_ind[3] + 1):block_ind[4]
  
  
  # Get number of variables in each block
  num_vars_block1 <- length(range1)
  num_vars_block2 <- length(range2)
  num_vars_block3 <- length(range3)
  num_vars_block4 <- length(range4)
  
  
  ### STEP 1 ###
  
  
  # Canonical Correlation Analysis (CCA)
  GlobalFactors <- blockfact0(Yorig, c(num_vars_block1, num_vars_block2, num_vars_block3,num_vars_block4), 1, c(1, 1, 1,1))
  GlobalFactors <- GlobalFactors / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(GlobalFactors) %*% GlobalFactors))))
  
  # Remove global component to get triplewise residuals
  
  Block123 <- cbind(Yorig[,range1],Yorig[,range2],Yorig[,range3])
  Block124 <- cbind(Yorig[,range1],Yorig[,range2],Yorig[,range4])
  Block134 <- cbind(Yorig[,range1],Yorig[,range3],Yorig[,range4])
  Block234 <- cbind(Yorig[,range2],Yorig[,range3],Yorig[,range4])
  
  Beta123 <-  beta_ols(GlobalFactors, Block123)
  Beta124 <-  beta_ols(GlobalFactors, Block124)
  Beta134 <-  beta_ols(GlobalFactors, Block134)
  Beta234 <-  beta_ols(GlobalFactors, Block234)
  
  Residuals123 <- Block123 - GlobalFactors %*% Beta123
  Residuals124 <- Block124 - GlobalFactors %*% Beta124
  Residuals134 <- Block134 - GlobalFactors %*% Beta134
  Residuals234 <- Block234 - GlobalFactors %*% Beta234
  
  # Apply CCA for lower-level block factors
  
  Factors123 <- blockfact0(Residuals123, c(num_vars_block1, num_vars_block2, num_vars_block3 ), 1, c(1, 1, 1))
  Factors124 <- blockfact0(Residuals124, c(num_vars_block1, num_vars_block2, num_vars_block4 ), 1, c(1, 1, 1))
  Factors134 <- blockfact0(Residuals134, c(num_vars_block1, num_vars_block3, num_vars_block4 ), 1, c(1, 1, 1))
  Factors234 <- blockfact0(Residuals234, c(num_vars_block2, num_vars_block3, num_vars_block4 ), 1, c(1, 1, 1))
  
  
  # Remove triplewise component to get pairwise residuals
  
  Block12 <- cbind(Yorig[,range1],Yorig[,range2])
  Block13 <- cbind(Yorig[,range1],Yorig[,range3])
  Block14 <- cbind(Yorig[,range1],Yorig[,range4])
  Block23 <- cbind(Yorig[,range2],Yorig[,range3])
  Block24 <- cbind(Yorig[,range2],Yorig[,range4])
  Block34 <- cbind(Yorig[,range3],Yorig[,range4])
  
  Beta12 <-  beta_ols(cbind(Factors123, Factors124), Block12)
  Beta13 <-  beta_ols(cbind(Factors123, Factors134), Block13)
  Beta14 <-  beta_ols(cbind(Factors124, Factors134), Block14)
  Beta23 <-  beta_ols(cbind(Factors123, Factors234), Block23)
  Beta24 <-  beta_ols(cbind(Factors124, Factors234), Block24)
  Beta34 <-  beta_ols(cbind(Factors134, Factors234), Block34)
  
  
  Residuals12 <- Block12 - cbind(Factors123, Factors124) %*% Beta12
  Residuals13 <- Block13 - cbind(Factors123, Factors134) %*% Beta13
  Residuals14 <- Block14 - cbind(Factors124, Factors134) %*% Beta14
  Residuals23 <- Block23 - cbind(Factors123, Factors234) %*% Beta23
  Residuals24 <- Block24 - cbind(Factors124, Factors234) %*% Beta24
  Residuals34 <- Block34 - cbind(Factors134, Factors234) %*% Beta34
  
  
  # Apply CCA for lower-level block factors
  Factors12 <- blockfact0(Residuals12, c(num_vars_block1, num_vars_block2), 1, c(1, 1))
  Factors13 <- blockfact0(Residuals13, c(num_vars_block1, num_vars_block3), 1, c(1, 1))
  Factors14 <- blockfact0(Residuals14, c(num_vars_block1, num_vars_block4), 1, c(1, 1))
  Factors23 <- blockfact0(Residuals23, c(num_vars_block2, num_vars_block3), 1, c(1, 1))
  Factors24 <- blockfact0(Residuals24, c(num_vars_block2, num_vars_block4), 1, c(1, 1))
  Factors34 <- blockfact0(Residuals34, c(num_vars_block3, num_vars_block4), 1, c(1, 1))
  
  # Calculate residuals for each block using the lower-level factors
  
  Block1 <- Yorig[, range1]
  Block2 <- Yorig[, range2]
  Block3 <- Yorig[, range3]
  Block4 <- Yorig[, range4]
  
  Beta1 <- beta_ols(cbind(Factors12, Factors13, Factors14), Block1)
  Beta2 <- beta_ols(cbind(Factors12, Factors23, Factors24), Block2)
  Beta3 <- beta_ols(cbind(Factors13, Factors23, Factors34), Block3)
  Beta4 <- beta_ols(cbind(Factors14, Factors24, Factors34), Block4)
  
  Residuals1 <- Block1 - cbind(Factors12, Factors13, Factors14) %*% Beta1
  Residuals2 <- Block2 - cbind(Factors12, Factors23, Factors24) %*% Beta2
  Residuals3 <- Block3 - cbind(Factors13, Factors23, Factors34) %*% Beta3
  Residuals4 <- Block4 - cbind(Factors14, Factors24, Factors34) %*% Beta4
 
  
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
  
  eig_vectors4 <- eigrs2(t(Residuals4) %*% Residuals4)$evec
  Factors4 <-  Residuals4 %*% eig_vectors4[, (ncol(eig_vectors4) - number_of_factor + 1):ncol(eig_vectors4)]
  Factors4 <- Factors4 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors4) %*% Factors4))))
  
  # Store all initial factors
  Factors <- cbind(GlobalFactors, Factors123, Factors124, Factors134, Factors234,
                   Factors12, Factors13, Factors14, Factors23, Factors24, Factors34, 
                   Factors1, Factors2, Factors3, Factors4)
  InitialFactors <- Factors
  
  
  # Compute initial loadings
  
  GlobalLoadings <- beta_ols(GlobalFactors, Yorig)
  Loadings123 <- beta_ols(Factors123, Block123)
  Loadings124 <- beta_ols(Factors124, Block124)
  Loadings134 <- beta_ols(Factors134, Block134)
  Loadings234 <- beta_ols(Factors234, Block234)
  Loadings12 <- beta_ols(Factors12, Block12)
  Loadings13 <- beta_ols(Factors13, Block13)
  Loadings14 <- beta_ols(Factors14, Block14)
  Loadings23 <- beta_ols(Factors23, Block23)
  Loadings24 <- beta_ols(Factors24, Block24)
  Loadings34 <- beta_ols(Factors34, Block34)
  Loadings1 <- beta_ols(Factors1, Block1)
  Loadings2 <- beta_ols(Factors2, Block2)
  Loadings3 <- beta_ols(Factors3, Block3)
  Loadings4 <- beta_ols(Factors4, Block4)
  
  Lambda <- matrix(0, nrow = num_factors, ncol = num_vars)
  Lambda[1,] <- GlobalLoadings
  Lambda[2,c(range1,range2,range3)] <- Loadings123
  Lambda[3,c(range1,range2,range4)] <- Loadings124
  Lambda[4,c(range1,range3,range4)] <- Loadings134
  Lambda[5,c(range2,range3,range4)] <- Loadings234
  Lambda[6,c(range1,range2)] <- Loadings12
  Lambda[7,c(range1,range3)] <- Loadings13
  Lambda[8,c(range1,range4)] <- Loadings14
  Lambda[9,c(range2,range3)] <- Loadings23
  Lambda[10,c(range2,range4)] <- Loadings24
  Lambda[11,c(range3,range4)] <- Loadings34
  
  Lambda[12,range1] <- Loadings1
  Lambda[13,range2] <- Loadings2
  Lambda[14,range3] <- Loadings3
  Lambda[15,range4] <- Loadings4
  
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
    Beta123 <-  beta_ols(GlobalFactors, Block123)
    Beta124 <-  beta_ols(GlobalFactors, Block124)
    Beta134 <-  beta_ols(GlobalFactors, Block134)
    Beta234 <-  beta_ols(GlobalFactors, Block234)
    
    Residuals123 <- Block123 - GlobalFactors %*% Beta123
    Residuals124 <- Block124 - GlobalFactors %*% Beta124
    Residuals134 <- Block134 - GlobalFactors %*% Beta134
    Residuals234 <- Block234 - GlobalFactors %*% Beta234
    
    # Run least squares to estimate triplewise factors
    Factors123 <- t(beta_ols(t(Loadings123), t(Residuals123)))
    Factors124 <- t(beta_ols(t(Loadings124), t(Residuals124)))
    Factors134 <- t(beta_ols(t(Loadings134), t(Residuals134)))
    Factors234 <- t(beta_ols(t(Loadings234), t(Residuals234)))
    
    
    # Regress to filter out the triplewise component
    
    Beta12 <-  beta_ols(cbind(Factors123, Factors124), Block12)
    Beta13 <-  beta_ols(cbind(Factors123, Factors134), Block13)
    Beta14 <-  beta_ols(cbind(Factors124, Factors134), Block14)
    Beta23 <-  beta_ols(cbind(Factors123, Factors234), Block23)
    Beta24 <-  beta_ols(cbind(Factors124, Factors234), Block24)
    Beta34 <-  beta_ols(cbind(Factors134, Factors234), Block34)
    
    
    Residuals12 <- Block12 - cbind(Factors123, Factors124) %*% Beta12
    Residuals13 <- Block13 - cbind(Factors123, Factors134) %*% Beta13
    Residuals14 <- Block14 - cbind(Factors124, Factors134) %*% Beta14
    Residuals23 <- Block23 - cbind(Factors123, Factors234) %*% Beta23
    Residuals24 <- Block24 - cbind(Factors124, Factors234) %*% Beta24
    Residuals34 <- Block34 - cbind(Factors134, Factors234) %*% Beta34
    
    # Run least squares to estimate pairwise factors
    Factors12 <- t(beta_ols(t(Loadings12), t(Residuals12)))
    Factors13 <- t(beta_ols(t(Loadings13), t(Residuals13)))
    Factors14 <- t(beta_ols(t(Loadings14), t(Residuals14)))
    Factors23 <- t(beta_ols(t(Loadings23), t(Residuals23)))
    Factors24 <- t(beta_ols(t(Loadings24), t(Residuals24)))
    Factors34 <- t(beta_ols(t(Loadings34), t(Residuals34)))
    
    # Regress to filter out the pairwise component
    Beta1 <- beta_ols(cbind(Factors12, Factors13, Factors14), Block1)
    Beta2 <- beta_ols(cbind(Factors12, Factors23, Factors24), Block2)
    Beta3 <- beta_ols(cbind(Factors13, Factors23, Factors34), Block3)
    Beta4 <- beta_ols(cbind(Factors14, Factors24, Factors34), Block4)
    
    Residuals1 <- Block1 - cbind(Factors12, Factors13, Factors14) %*% Beta1
    Residuals2 <- Block2 - cbind(Factors12, Factors23, Factors24) %*% Beta2
    Residuals3 <- Block3 - cbind(Factors13, Factors23, Factors34) %*% Beta3
    Residuals4 <- Block4 - cbind(Factors14, Factors24, Factors34) %*% Beta4
    
    
    # Run least squares to estimate single block factors
    Factors1 <- t(beta_ols(t(Loadings1), t(Residuals1)))
    Factors2 <- t(beta_ols(t(Loadings2), t(Residuals2)))
    Factors3 <- t(beta_ols(t(Loadings3), t(Residuals3)))
    Factors4 <- t(beta_ols(t(Loadings4), t(Residuals4)))
    
    
    # Update and normalize factors 
    GlobalFactors <- GlobalFactors / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(GlobalFactors) %*% GlobalFactors))))

    Factors123 <- Factors123 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors123) %*% Factors123))))
    Factors124 <- Factors124 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors124) %*% Factors124))))
    Factors134 <- Factors134 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors134) %*% Factors134))))
    Factors234 <- Factors234 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors234) %*% Factors234))))
    
    
    Factors12 <- Factors12 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors12) %*% Factors12))))
    Factors13 <- Factors13 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors13) %*% Factors13))))
    Factors14 <- Factors14 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors14) %*% Factors14))))
    Factors23 <- Factors23 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors23) %*% Factors23))))
    Factors24 <- Factors24 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors24) %*% Factors24))))
    Factors34 <- Factors34 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors34) %*% Factors34))))
    

    Factors1 <- Factors1 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors1) %*% Factors1))))
    Factors2 <- Factors2 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors2) %*% Factors2))))
    Factors3 <- Factors3 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors3) %*% Factors3))))
    Factors4 <- Factors4 / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors4) %*% Factors4))))
    
    Factors <- cbind(GlobalFactors, Factors123, Factors124, Factors134, Factors234,
                     Factors12, Factors13, Factors14, Factors23, Factors24, Factors34, 
                     Factors1, Factors2, Factors3, Factors4)
    
    # Compute and update loadings
    #Lambda <- beta_ols(Factors, Yorig)
    
    GlobalLoadings <- beta_ols(GlobalFactors, Yorig)
    Loadings123 <- beta_ols(Factors123, Block123)
    Loadings124 <- beta_ols(Factors124, Block124)
    Loadings134 <- beta_ols(Factors134, Block134)
    Loadings234 <- beta_ols(Factors234, Block234)
    Loadings12 <- beta_ols(Factors12, Block12)
    Loadings13 <- beta_ols(Factors13, Block13)
    Loadings14 <- beta_ols(Factors14, Block14)
    Loadings23 <- beta_ols(Factors23, Block23)
    Loadings24 <- beta_ols(Factors24, Block24)
    Loadings34 <- beta_ols(Factors34, Block34)
    Loadings1 <- beta_ols(Factors1, Block1)
    Loadings2 <- beta_ols(Factors2, Block2)
    Loadings3 <- beta_ols(Factors3, Block3)
    Loadings4 <- beta_ols(Factors4, Block4)
    
    
    Lambda[1,] <- GlobalLoadings
    Lambda[2,c(range1,range2,range3)] <- Loadings123
    Lambda[3,c(range1,range2,range4)] <- Loadings124
    Lambda[4,c(range1,range3,range4)] <- Loadings134
    Lambda[5,c(range2,range3,range4)] <- Loadings234
    Lambda[6,c(range1,range2)] <- Loadings12
    Lambda[7,c(range1,range3)] <- Loadings13
    Lambda[8,c(range1,range4)] <- Loadings14
    Lambda[9,c(range2,range3)] <- Loadings23
    Lambda[10,c(range2,range4)] <- Loadings24
    Lambda[11,c(range3,range4)] <- Loadings34
    
    Lambda[12,range1] <- Loadings1
    Lambda[13,range2] <- Loadings2
    Lambda[14,range3] <- Loadings3
    Lambda[15,range4] <- Loadings4
    

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


