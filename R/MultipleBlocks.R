# Multi-Level Dynamic Factor Model - Multiple blocks

beta_ols <- function(X, Y) {
  solve(t(X) %*% X) %*% t(X) %*% Y
}

# Function to orthogonalize the factors 
orthogonalize_factors <- function(X) {
  n_factors <- ncol(X)  # Number of factors 
  n_obs <- nrow(X)      # Number of observations 
  
  # Initialize a matrix to store the orthogonalized factors
  orthogonal_factors <- matrix(0, nrow = n_obs, ncol = n_factors)
  
  # first factor remains 
  orthogonal_factors[, 1] <- X[, 1]
  
  # Iterative orthogonalization 
  for (i in 2:n_factors) {
    
    # Regress factor i on the previous orthogonalized factors
    model <- lm(X[, i] ~ orthogonal_factors[, 1:(i - 1)])
    
    # Compute residuals to get new factor
    orthogonal_factors[, i] <- residuals(model)
  }
  
  # Return
  return(orthogonal_factors)
}

# Provide previous level factors 
get_Factors <- function(Factor_list, combination, len) {
 
  matching_values <- list()
  #comb_str <- as.character(combination)
  
  for (key in names(Factor_list)) {
    key_elements <- as.integer(strsplit(key, "-")[[1]])
    
    if (length(key_elements) == (len + 1) && all(combination %in% key_elements)) {
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


MultipleBlocks<-function(Yorig,r,block_ind,tol,max_iter){

  # Standardize the original data
  Yorig <- scale(Yorig)
  
  # Initialize variables
  results <- list()  # List to save results
  num_blocks <- length(block_ind) # Number of blocks
  num_obs <- nrow(Yorig) # Total number of observations
  num_factors <- 2^num_blocks-1 # Total number of factors
  num_factors <- sum(r)
  
  
  # Define block ranges and count the number of var in each range
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
  
  # Define data structure
  Factor_list <- list()
  InitialFactors <- matrix(nrow = num_obs, ncol = 0)  # start with zero columns
  
  
  Loadings_list <- list()
  Lambda <- matrix(0, nrow = num_factors, ncol =  ncol(Yorig))
  counter <- 1 # Lambda matrix index
  
  
  ### STEP 1 ###
  factor_index <- 1
  
  # Compute Global factors
  if(r[factor_index]>0){
    GlobalFactors <- blockfact0(Yorig, num_vars, r[factor_index], rep(1, length(num_vars)))
    GlobalFactors <- GlobalFactors / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(GlobalFactors) %*% GlobalFactors))))
    
    # Store Global factors
    key <- paste(seq(1, num_blocks), collapse = "-")  
    Factor_list[[key]] <- GlobalFactors  
    InitialFactors <- cbind(InitialFactors, GlobalFactors)
    
    
    # Compute Global Loadings
    Loadings <- beta_ols(GlobalFactors, Yorig)
    
    
    
    # Store Global Loadings
    Loadings_list[[key]] <- Loadings  
    combination <- seq(1, num_blocks)
    Lambda[counter:(counter+r[factor_index]-1), unlist(ranges[combination])] <- Loadings
    counter <- counter + r[factor_index]
  }
  
  
  factor_index <- factor_index + 1

  # Loop on lower levels to compute Factors and Loadings
  for (i in 1:(num_blocks-1)) {
    k <-  num_blocks - i
    combinations_matrix <- t(combn(num_blocks,k))
    for (j in 1:nrow(combinations_matrix)) {
      combination <- combinations_matrix[j, ]
      
     
      
      if(r[factor_index]>0){
        
        # Extract block
        Block <- do.call(cbind, lapply(combination, function(idx) Yorig[, ranges[[idx]]]))
        
        # Get previous level factors
        
        len <- length(combination)
        Factors <- NULL
        while (is.null(Factors) && len < num_blocks) {
          Factors <- get_Factors(Factor_list, combination, len)
          len <- len + 1
        }
        
        if (is.null(Factors)) {
          message("Factors were not found before reaching the maximum len value.")
          return(1)
        }
        
        # Compute residuals
        Beta <-  beta_ols(Factors, Block)
        Residuals <- Block - Factors %*% Beta
          
        # Compute factors
        if (i < num_blocks -1){
          # Middle levels
          Factors <- blockfact0(Residuals, num_vars[combination], r[factor_index], rep(1, length(combination)))
        }else{
          # Last level - Apply simple PCA
          number_of_factor <- r[factor_index] # number of factor to be extracted with PCA
          eig_vectors <- eigen_sorted(t(Residuals) %*% Residuals)$eigenvectors
          Factors <-  Residuals %*% eig_vectors[, (ncol(eig_vectors) - number_of_factor + 1):ncol(eig_vectors)]
          Factors <- Factors / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors) %*% Factors))))
        }
        
        
        # Store Factors
        key <- paste(combination, collapse = "-")
        Factor_list[[key]] <- Factors  
        InitialFactors <- cbind(InitialFactors, Factors)
        
        
        # Compute Loadings
        
        Loadings <- beta_ols(Factors, Block)
        
        # Store Loadings
        Loadings_list[[key]] <- Loadings  
        Lambda[counter:(counter+r[factor_index]-1), unlist(ranges[combination])] <- Loadings
        counter <- counter + r[factor_index]
      
      }
     
    
      factor_index = factor_index + 1
     
    }
    
  }
  
  InitialLoadings <- Lambda
  

  
  
  ### STEP 2 ###
  
  # Initialize residual sum of squares (RSS) for convergence
  iteration <- 0
  Residuals <- Yorig - InitialFactors %*% Lambda
  RSS_previous <- sum(diag(t(Residuals) %*% Residuals))
  
 
  
  # Iterative procedure for convergence
  while (iteration < max_iter) {
    
    
    
    factor_index <- 1
    counter <- 1
    
    iteration <- iteration + 1
    
    
    FinalFactors <- matrix(nrow = num_obs, ncol = 0)  
    
    
    if(r[factor_index]>0){
      
      # Extract Global Loadings
      key <- paste(seq(1, num_blocks), collapse = "-")
      Loadings <- Loadings_list[[key]]
      
      Factors <- t(beta_ols(t(Loadings), t(Yorig)))
      Factors <- Factors / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors) %*% Factors))))
    
      # Store new Global Factors
      Factor_list[[key]] <- Factors 
      FinalFactors <- cbind(FinalFactors, Factors)
      
      Loadings <- beta_ols(Factors, Yorig)
      
      # Store new Global Loadings
      Loadings_list[[key]] <- Loadings
      Lambda[counter:(counter+r[factor_index]-1),] <- Loadings
      counter <- counter + r[factor_index]
      
    }
    
    
    factor_index = factor_index + 1
    
    
    
   
    
    # Lower levels
   
    for (i in 1:(num_blocks-1)) {
      k <-  num_blocks - i
      combinations_matrix <- t(combn(num_blocks,k))
      for (j in 1:nrow(combinations_matrix)) {
        combination <- combinations_matrix[j,]
        
        
        
        if(r[factor_index]>0){
          
          # Extract block
          Block <- do.call(cbind, lapply(combination, function(idx) Yorig[, ranges[[idx]]]))
          
          # Get previous level Factors 

          len <- length(combination)
          Factors <- NULL
          while (is.null(Factors) && len < num_blocks) {
            Factors <- get_Factors(Factor_list, combination, len)
            len <- len + 1
          }
          
          if (is.null(Factors)) {
            message("Factors were not found before reaching the maximum len value.")
            return(1)
          }
          
          # Compute residuals
          Beta <-  beta_ols(Factors, Block)
          Residuals <- Block - Factors %*% Beta
          
          
          # Extract Loadings
          key <- paste(combination, collapse = "-")
          Loadings <- Loadings_list[[key]]
          

          # Compute new Factors
          Factors <- t(beta_ols(t(Loadings), t(Residuals)))
          Factors <- Factors / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors) %*% Factors))))
        
          
          # Store new Factors
          key <- paste(combination, collapse = "-")
          Factor_list[[key]] <- Factors  
          FinalFactors <- cbind(FinalFactors, Factors)
          
          # Compute new Loadings
          Loadings <- beta_ols(Factors, Block)

          
          # Store new Loadings
          Loadings_list[[key]] <- Loadings
          Lambda[counter:(counter+r[factor_index]-1), unlist(ranges[combination])] <- Loadings
          
          
          counter <- counter + r[factor_index]
          
          
        }
    
        factor_index = factor_index + 1
        
      }
    }
    
    
    # Check RSS
    Residuals <- Yorig - FinalFactors %*% Lambda
    RSS_new <- sum(diag(t(Residuals) %*% Residuals))
    
    if (abs(RSS_previous - RSS_new) < tol) {
      break  # Converged
    }
    
    # Update RSS
    RSS_previous <- RSS_new
  }
  
  #Final RSS
  print(RSS_previous)
  
  # Factor Orthogonalization
  orthogonal_FinalFactors <- orthogonalize_factors(FinalFactors)
  
  
  # Store results
  results[["Factors"]] <- orthogonal_FinalFactors
  results[["Loadings"]] <- Lambda
  
  # Return
  return(results)
}


