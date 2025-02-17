# Multi-Level Dynamic Factor Model - Multiple blocks

library(MASS)

beta_ols <- function(X, Y) {
  solve(t(X) %*% X) %*% t(X) %*% Y
}



MultipleBlocks4<-function(Yorig,r,block_ind,tol,max_iter,method){

  ### PRE-PREOCESSING ###
  
  # Standardize the original data
  Yorig <- scale(Yorig,TRUE,TRUE)
 
  # Initialize variables
  results <- list()  # List to save results
  num_blocks <- length(block_ind) # Number of blocks
  num_obs <- nrow(Yorig) # Total number of observations
  num_factors <- sum(r) # Total number of factors
  
  
  # Define block ranges and count the number of variables in each range
  ranges <- list()
  num_vars <- numeric(length(block_ind))  
 
  for (i in 1:length(block_ind)) {
    if (i == 1) {
      ranges[[i]] <- 1:block_ind[i]
    } else {
      ranges[[i]] <- (block_ind[i - 1] + 1):block_ind[i]
    }
    
    num_vars[i] <- length(ranges[[i]]) 
  }
  
  
  ### STEP 1 - COMPUTE INITIAL FACTORS ###
  res <- Compute_Initial_Factors(Yorig, num_vars, num_obs, num_blocks, ranges, num_factors, r, method)
  InitialFactors <- res$InitialFactors
  Factor_list <- res$Factor_list
  
  
  ### STEP 2 - ITERATIVE PROCEDURE ###
  RSS_previous <- 1000000000000
  rss_values <- c()
  iteration <- 0


  # Iterative procedure for convergence
  while (iteration < max_iter) {
    iteration <- iteration + 1
    print(iteration)

    
    # Compute Lambda
    Lambda <- Compute_Lambda(Yorig,num_blocks,ranges,num_factors,r,Factor_list)
    
    # Compute new factors
    FinalFactors <- t(solve(Lambda %*% t(Lambda)) %*% Lambda %*% t(Yorig))
    
    # Update factor list
    Factor_list <- update_factor_list(Factor_list, FinalFactors, r)
      

    # Check RSS
    FinalResiduals <- Yorig - FinalFactors %*% Lambda
    RSS_new <- sum(FinalResiduals^2)
    print(RSS_new)
    rss_values <- c(rss_values,RSS_new)
    
    if ((log(RSS_previous) - log(RSS_new)) < tol) {
      break  # Converged
    }

    RSS_previous <- RSS_new

  }

  
  
  
  
  
 
  #print(RSS_new)
  
  
  #orthogonal_InitialFactors <- orthogonalize_factors(InitialFactors) # orthogonalization
  orthogonal_InitialFactors <- InitialFactors
  #orthogonal_FinalFactors <- orthogonalize_factors(FinalFactors) # orthogonalization
  orthogonal_FinalFactors <- FinalFactors
  
  
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(-1,-1,1)) #hcpi
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(-1,-1,-1)) #ccpi
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(-1,-1,1)) #ecpi
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(-1,-1,-1,-1)) #ecpi
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(1,1,1,1)) #ecpi
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(-1,-1,-1)) #fcpi
  
  
  # Scale factors
  orthogonal_InitialFactors <- scale(orthogonal_InitialFactors,TRUE,TRUE)
  orthogonal_FinalFactors <- scale(orthogonal_FinalFactors,TRUE,TRUE)
  
 
  # Factor list
  Final_list <- list()
  for (key in names(Factor_list)) {
    
    factors <- Factor_list[[key]]
    n_factors <-  ncol(factors)
    
    Final_list[[key]] <- n_factors
  }
  
  
  
  # Compute factor hat 
  Factors_hat <- matrix(nrow = num_obs, ncol = 0) 
  
  
  factor_index <- 1
 
  for (key in names(Final_list)){
    
    # extract combination
    combination <- as.numeric(unlist(strsplit(key, "-")))
    
    # Extract block data 
    Block <- do.call(cbind, lapply(combination, function(idx) Yorig[, ranges[[idx]]]))
   
    # Extract corresponding factors
    n_factors <- Final_list[[key]]
    Facts <- orthogonal_FinalFactors[,factor_index:(factor_index+n_factors-1)]
    
    # Compute Loadings 
    Loads <- beta_ols(Facts, Block)
    Lambda[factor_index:(factor_index+n_factors-1), unlist(ranges[combination])] <- Loads
   
    # Compute Residuals
    Resid <- Block - Facts %*% Loads
    
    # Compute Factors Hat
    N <- ncol(Loads)
    Facts_hat<-(1/N)*Block%*%t(Loads)
    
    
    Factors_hat <- cbind(Factors_hat, Facts_hat)
    
    factor_index <- factor_index + n_factors
    
  }
 
  # Compute final residuals
  Residuals <- Yorig - orthogonal_FinalFactors %*% Lambda
  
  
  # Store results
  results[["Initial_Factors"]] <- orthogonal_InitialFactors
  results[["Factors"]] <- orthogonal_FinalFactors
  results[["Factors_hat"]] <- Factors_hat
  results[['Lambda']] <- t(Lambda)
  results[['Residuals']] <- Residuals
  results[['Factors_list']] <- Final_list
  results[['RSS_list']] <- rss_values
  
 
  return(results)
}


