# Multi-Level Dynamic Factor Model - Multiple blocks

library(MASS)


beta_ols <- function(X, Y) {
  solve(t(X) %*% X) %*% t(X) %*% Y
}

# Orthogonalize factors
orthogonalize_factors <- function(X) {
  n_factors <- ncol(X)  
  
  # Iterative orthogonalization
  for (i in 1:n_factors) {
    # Regress the i-th factor on all other factors (excluding itself)
    other_factors <- X[, -i, drop = FALSE]  
    model <- lm(X[, i] ~ other_factors)
    
    # Update the i-th factor with the residuals 
    X[, i] <- residuals(model)
  }
  
  return(X)
}


# Extract factors from a given level
get_Factors <- function(Factor_list, combination, level) {
 
  matching_values <- list()
  
  
  for (key in names(Factor_list)) {
    key_elements <- as.integer(strsplit(key, "-")[[1]])
    
    if (length(key_elements) == (level) && all(combination %in% key_elements)) {
      matching_values <- cbind(matching_values, Factor_list[[key]])
    }
  }
  
  matching_values_numeric <- as.numeric(matching_values)
  original_dims <- dim(matching_values)
  
  if (is.null(original_dims)){
    matching_values <- NULL
  }else{
    matching_values <- matrix(matching_values_numeric, nrow = original_dims[1], ncol = original_dims[2])
    
  }

  return(matching_values)
}




MultipleBlocks_3<-function(Yorig,r,block_ind,tol,max_iter,method){

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
  
  
  # Define Factors data structures
  Factor_list <- list()
  Loading_list <- list()
  InitialFactors <- matrix(nrow = num_obs, ncol = 0)  
  
  # Loadings 
  Lambda <- matrix(0, nrow = num_factors, ncol =  ncol(Yorig))
  r_index <- 1 
  counter <- 1
  
  ### STEP 1 ###
  
  
  # Compute Global factors
  if(r[r_index]>0){
    
    number_of_factor <- r[r_index] # number of factor to be extracted with PCA
    if (method == 0){
      # CCA
      GlobalFactors <- blockfact0(Yorig, num_vars, number_of_factor, rep(1, length(num_vars)))
    }else{
      # PCA 
      pca_result <- prcomp(Yorig, scale. = FALSE)
      GlobalFactors <- pca_result$x[, 1:number_of_factor]
      GlobalFactors <- GlobalFactors / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(GlobalFactors) %*% GlobalFactors))))
    }
    
    
    # Store Global factors
    key <- paste(seq(1, num_blocks), collapse = "-")  
    Factor_list[[key]] <- GlobalFactors  
    InitialFactors <- cbind(InitialFactors, GlobalFactors)
    
    
    # Compute Global Loadings
    GlobalLoadings <- beta_ols(GlobalFactors, Yorig)
    Loading_list[[key]] <- GlobalLoadings  
    
    # Update Lambda
    combination <- seq(1, num_blocks)
    Lambda[counter:(counter+r[r_index]-1), unlist(ranges[combination])] <- GlobalLoadings
    counter <- counter + r[r_index]
  }
  
  

  # Loop on lower levels to compute Factors 
  for (i in 1:(num_blocks-1)) {
    k <-  num_blocks - i
    combinations_matrix <- t(combn(num_blocks,k))
    for (j in 1:nrow(combinations_matrix)) {
      combination <- combinations_matrix[j, ]
      
      
      r_index <- r_index + 1
      
      # Skip blocks where Factors are not needed
      if (r[r_index] == 0){
        next
      }
        
     
      #Extract Residuals filtering out upper levels factors (start with global factors)
      level <- num_blocks

      Residuals <- do.call(cbind, lapply(combination, function(idx) Yorig[, ranges[[idx]]])) # initialize with block data

      while (level > length(combination)) {
        Factors <- get_Factors(Factor_list, combination, level)

        # filter out
        if(!is.null(Factors)){
          ols_result <- beta_ols(Factors, Residuals)
          Residuals <- Residuals - Factors %*% ols_result
        }

        level <- level - 1

      }
          
      number_of_factor <- r[r_index] # number of factor to be extracted with PCA
      
      # Compute factors
      if (i < num_blocks - 1 && method == 0) {
        # Use CCA for middle level
        Factors <- blockfact0(Residuals, num_vars[combination], number_of_factor, rep(1, length(combination)))
      }else{
        # Use PCA
        pca_result <- prcomp(Residuals, scale. = FALSE)
        Factors <- pca_result$x[, 1:number_of_factor]
        Factors <- Factors / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors) %*% Factors))))
        
      }
      
      # Store Factors
      key <- paste(combination, collapse = "-")
      Factor_list[[key]] <- Factors  
      InitialFactors <- cbind(InitialFactors, Factors)
      
      # Compute Loadings
      Loadings <- beta_ols(Factors, Residuals)
      Loading_list[[key]] <- Loadings  
      
      
      # Update Lambda
      Lambda[counter:(counter+r[r_index]-1), unlist(ranges[combination])] <- Loadings
      counter <- counter + r[r_index]

    }
    
  }
  
  
  
 
  
  
  
  
  
 

  
  

  RSS_previous <- 1000000000000
  rss_values <- c()
  iteration <- 0


  # Iterative procedure for convergence
  while (iteration < max_iter) {
    iteration <- iteration + 1

    FinalFactors <- matrix(nrow = num_obs, ncol = 0)  
    
    r_index <- 1 
    counter <- 1
    
    
    # Extract Global Loadings
    key <- paste(seq(1, num_blocks), collapse = "-")
    GlobalLoadings <- Loading_list[[key]]
    
    # Compute Global Factors
    GlobalFactors <- t(beta_ols(t(GlobalLoadings), t(Yorig)))
    #GlobalFactors <- Yorig %*% t(GlobalLoadings) %*% solve(GlobalLoadings %*% t(GlobalLoadings))
    #GlobalFactors <- t(solve(GlobalLoadings %*% t(GlobalLoadings)) %*% GlobalLoadings %*% t(Yorig))
    
    # Store Global factors
    key <- paste(seq(1, num_blocks), collapse = "-")  
    Factor_list[[key]] <- GlobalFactors  
    FinalFactors <- cbind(FinalFactors, GlobalFactors)
    
    
    # Compute Global Loadings
    GlobalLoadings <- beta_ols(GlobalFactors, Yorig)
    Loading_list[[key]] <- GlobalLoadings  
   
    
    # Update Lambda
    combination <- seq(1, num_blocks)
    Lambda[counter:(counter+r[r_index]-1), unlist(ranges[combination])] <- GlobalLoadings
    counter <- counter + r[r_index]
    
    


    # Loop on lower levels
    for (i in 1:(num_blocks-1)) {
      k <-  num_blocks - i
      combinations_matrix <- t(combn(num_blocks,k))
      for (j in 1:nrow(combinations_matrix)) {
        combination <- combinations_matrix[j,]

        r_index <- r_index + 1

        # Skip blocks where Factors are not needed
        if (r[r_index] == 0){
          next
        }


        # Extract Residuals filtering out upper levels factors
        level <- num_blocks

        Residuals <- do.call(cbind, lapply(combination, function(idx) Yorig[, ranges[[idx]]]))


        while (level > length(combination)) {
          Factors <- get_Factors(Factor_list, combination, level)

          # filter out
          if(!is.null(Factors)){
            ols_result <- beta_ols(Factors, Residuals)
            Residuals <- Residuals - Factors %*% ols_result
          }

          level <- level - 1

        }

        

        
        # Extract current block Loadings
        key <- paste(combination, collapse = "-")
        Loadings <- Loading_list[[key]]
        
        
        # Compute new Factors
        Factors <- t(beta_ols(t(Loadings), t(Residuals)))
        
        # Store Factors
        Factor_list[[key]] <- Factors  
        FinalFactors <- cbind(FinalFactors, Factors)
        
        # Compute Loadings
        Loadings <- beta_ols(Factors, Residuals)
        Loading_list[[key]] <- Loadings  
        
        
        # Update Lambda
        Lambda[counter:(counter+r[r_index]-1), unlist(ranges[combination])] <- Loadings
        counter <- counter + r[r_index]

      }
    }

    
    
  


    # Check RSS
    FinalResiduals <- Yorig - FinalFactors %*% Lambda
    #RSS_new <- sum(diag(t(FinalResiduals) %*% FinalResiduals))
    RSS_new <- sum(FinalResiduals^2)
    #print(RSS_new)
    
    rss_values <- c(rss_values,RSS_new)
    
    
    #if (abs(RSS_previous - RSS_new) / RSS_previous < tol) {
    if ((log(RSS_previous) - log(RSS_new)) < tol) {
      break  # Converged
    }

    # Update RSS
    RSS_previous <- RSS_new


  }

  #Final RSS
  #print('Final RSS')
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


