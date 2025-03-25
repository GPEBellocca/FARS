# Multi-Level Dynamic Factor Model - Multiple blocks

#library(MASS)

beta_ols <- function(X, Y) {
  solve(t(X) %*% X) %*% t(X) %*% Y
}



MultipleBlocks<-function(Yorig,r,block_ind,tol,max_iter,method){

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
  Loadings_list <- list()
  while (iteration < max_iter) {
    iteration <- iteration + 1
    #print(iteration)

    # Compute Lambda
    L_res <- Compute_Lambda(Yorig,num_blocks,ranges,num_factors,r,Factor_list,Loadings_list)
    Lambda <- L_res$Lambda
    Loadings_list <- L_res$Loadings_list
    
    # Compute new factors
    FinalFactors <- t(solve(Lambda %*% t(Lambda)) %*% Lambda %*% t(Yorig))
    
    # Update factor list
    Factor_list <- update_factor_list(Factor_list, FinalFactors, r)
      

    # Check RSS
    FinalResiduals <- Yorig - FinalFactors %*% Lambda
    RSS_new <- sum(FinalResiduals^2)
    #print(RSS_new)
    rss_values <- c(rss_values,RSS_new)
    
    if ((log(RSS_previous) - log(RSS_new)) < tol) {
      break  # Converged
    }

    RSS_previous <- RSS_new

  }

 
  #print(RSS_new)
  
  
  # Impose identification on final factors and lambda
  Id_result <- PC_identifications(Yorig, num_blocks, ranges, num_factors, r, FinalFactors, Factor_list, Loadings_list)
  Factor_list <- Id_result$Factor_list
  orthogonal_FinalFactors <- Id_result$FinalFactors
  Lambda <- Id_result$Lambda
  
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(-1,-1,1)) #hcpi
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(-1,-1,-1)) #ccpi
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(-1,-1,1)) #ecpi
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(-1,-1,-1,-1)) #ecpi
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(1,1,1,1)) #ecpi
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(-1,-1,-1)) #fcpi
  
  
  # Compute final residuals
  Residuals <- Yorig - orthogonal_FinalFactors %*% Lambda
  
  
  # Factor list
  Final_list <- list()
  for (key in names(Factor_list)) {
    factors <- Factor_list[[key]]
    n_factors <-  ncol(factors)
    Final_list[[key]] <- n_factors
  }
  
  
  # Factors hat
  Factors_hat <- Compute_factors_hat(Yorig, ranges, Final_list, Factor_list,Loadings_list)
    
  

  # # Compute factor hat
  # Factors_hat <- matrix(nrow = num_obs, ncol = 0)
  # factor_index <- 1
  # 
  # for (key in names(Final_list)){
  # 
  #   # extract combination
  #   combination <- as.numeric(unlist(strsplit(key, "-")))
  # 
  #   # Extract block data
  #   Block <- do.call(cbind, lapply(combination, function(idx) Yorig[, ranges[[idx]]]))
  # 
  #   # Extract corresponding factors
  #   n_factors <- Final_list[[key]]
  #   Facts <- orthogonal_FinalFactors[,factor_index:(factor_index+n_factors-1)]
  # 
  #   # Compute Loadings
  #   Loads <- beta_ols(Facts, Block)
  #   Lambda[factor_index:(factor_index+n_factors-1), unlist(ranges[combination])] <- Loads
  # 
  #   # Compute Residuals
  #   Resid <- Block - Facts %*% Loads
  # 
  #   # Compute Factors Hat
  #   N <- ncol(Loads)
  #   Facts_hat<-(1/N)*Block%*%t(Loads)
  # 
  # 
  #   Factors_hat <- cbind(Factors_hat, Facts_hat)
  #   factor_index <- factor_index + n_factors
  # 
  # }
  
 

  
  # Store results
  results[["Initial_Factors"]] <- InitialFactors
  results[["Factors"]] <- orthogonal_FinalFactors
  results[["Factors_hat"]] <- Factors_hat
  results[['Lambda']] <- t(Lambda)
  results[['Residuals']] <- Residuals
  results[['Factors_list']] <- Final_list
  results[['RSS_list']] <- rss_values
  
 
  return(results)
}


