# Multi-Level Dynamic Factor Model - Multiple blocks

library(MASS)

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


orthogonalize_factors2 <- function(X) {
  n_factors <- ncol(X)  # Number of factors (columns)
  
  # Iterative orthogonalization
  for (i in 1:n_factors) {
    # Regress the i-th factor on all other factors (excluding itself)
    other_factors <- X[, -i, drop = FALSE]  # Drop the ith column
    model <- lm(X[, i] ~ other_factors)
    
    # Update the i-th factor with the residuals (orthogonalized)
    X[, i] <- residuals(model)
  }
  
  # Return the matrix with orthogonalized factors
  return(X)
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


MultipleBlocks4<-function(Yorig,r,block_ind,tol,max_iter){

  # Standardize the original data
  Yorig <- scale(Yorig,TRUE,TRUE)
 
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
  
  
  ### STEP 1 ###
  factor_index <- 1
  
  # Compute Global factors
  if(r[factor_index]>0){
    GlobalFactors <- blockfact0(Yorig, num_vars, r[factor_index], rep(1, length(num_vars)))
    
    # Store Global factors
    key <- paste(seq(1, num_blocks), collapse = "-")  
    Factor_list[[key]] <- GlobalFactors  
    InitialFactors <- cbind(InitialFactors, GlobalFactors)
  }
  
  

  # Loop on lower levels to compute Factors 
  for (i in 1:(num_blocks-1)) {
    k <-  num_blocks - i
    combinations_matrix <- t(combn(num_blocks,k))
    for (j in 1:nrow(combinations_matrix)) {
      combination <- combinations_matrix[j, ]
      
      factor_index = factor_index + 1
      
      # Skip blocks where Factors are not needed
      if (r[factor_index] == 0){
        next
      }
        
      #----------------------------------------------------
      # # Extract block
      # Block <- do.call(cbind, lapply(combination, function(idx) Yorig[, ranges[[idx]]]))
      # 
      # # Get previous level factors
      # len <- length(combination)
      # Factors <- get_Factors(Factor_list, combination, len)
      # 
      # # Compute residuals
      # if(is.null(Factors)){
      #   Residuals <- Block
      # }else{
      #   Beta <-  beta_ols(Factors, Block)
      #   Residuals <- Block - Factors %*% Beta
      # }

      # ----------------------------------------------------
      
      #Extract Residuals filtering out upper levels factors
      len <- num_blocks

      Residuals <- do.call(cbind, lapply(combination, function(idx) Yorig[, ranges[[idx]]])) # initialize with block data

      while (len > length(combination)) {
        Factors <- get_Factors(Factor_list, combination, len-1)

        # filter out
        if(!is.null(Factors)){
          ols_result <- beta_ols(Factors, Residuals)
          Residuals <- Residuals - Factors %*% ols_result
        }

        len <- len - 1

      }
          
      # ----------------------------------------------------
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

    }
    
  }
  
 
  
  
  RSS_previous <- 1000000000000
 
  iteration <- 0
  
  
  # Iterative procedure for convergence
  while (iteration < max_iter) {
    iteration <- iteration + 1
    #print(iteration)
    
    r_index <- 1 
    counter <- 1
    
    # Initialize Lambda
    Lambda <- matrix(0, nrow = num_factors, ncol =  ncol(Yorig))
    
    # Extract Global Factors
    key <- paste(seq(1, num_blocks), collapse = "-")
    GlobalFactors <- Factor_list[[key]]
    
    # Compute Global Loadings
    GlobalLoadings <- beta_ols(GlobalFactors, Yorig)
    
    # Update Lambda
    combination <- seq(1, num_blocks)
    Lambda[counter:(counter+r[r_index]-1), unlist(ranges[combination])] <- GlobalLoadings
    counter <- counter + r[r_index]
    
    
    
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
        len <- num_blocks
       
        Residuals <- do.call(cbind, lapply(combination, function(idx) Yorig[, ranges[[idx]]]))
        
        
        while (len > length(combination)) {
          Factors <- get_Factors(Factor_list, combination, len-1)
          
          # filter out
          if(!is.null(Factors)){
            ols_result <- beta_ols(Factors, Residuals)
            Residuals <- Residuals - Factors %*% ols_result
          }
          
          len <- len - 1
          
        }
        
        # Extract current block factors
        key <- paste(combination, collapse = "-")
        Factors <- Factor_list[[key]]
        
        # Compute Loadings
        Loadings <- beta_ols(Factors, Residuals)
        
        
        # Update Lambda
        Lambda[counter:(counter+r[r_index]-1), unlist(ranges[combination])] <- Loadings
        counter <- counter + r[r_index]
        
      }
    }
    
   
    
    # Compute new factors
    FinalFactors <- t(solve(Lambda %*% t(Lambda)) %*% Lambda %*% t(Yorig))
    
    
    
    
    r_index <- 1 
    counter <- 1
    
    filtered_r <- r[r != 0]
    # Update  Factors list
    for (key in names(Factor_list)) {
     
      select_factors <- FinalFactors[,counter:(counter+filtered_r[r_index]-1)]
      
      Factor_list[[key]] <- matrix(select_factors, ncol = filtered_r[r_index])
      
      counter <- counter + filtered_r[r_index]
      r_index <- r_index + 1
    }
    
   
    
    # Check RSS
    FinalResiduals <- Yorig - FinalFactors %*% Lambda
    RSS_new <- sum(diag(t(FinalResiduals) %*% FinalResiduals))
    
   
    
    
    if ((log(RSS_previous) - log(RSS_new)) < tol) {
      break  # Converged
    }
    
    # Update RSS
    RSS_previous <- RSS_new
    
   
  }
  
  #Final RSS
  print(RSS_new)
  
 

 
  # Factor Orthogonalization
  #orthogonal_FinalFactors <- orthogonalize_factors(FinalFactors)
  orthogonal_FinalFactors <- orthogonalize_factors2(FinalFactors)
  
 
  
  
  # 
  Final_list <- list()
 
  for (key in names(Factor_list)) {
    
    
    factors <- Factor_list[[key]]
    n_factors <-  ncol(factors)
    
    Final_list[[key]] <- n_factors
    
  }
  
  
  
  
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(-1,-1,-1)) #hcpi
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(-1,-1,-1)) #ccpi
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(-1,-1,1)) #ecpi
  #orthogonal_FinalFactors <- orthogonal_FinalFactors %*% diag(c(-1,-1,-1)) #fcpi
  
  orthogonal_FinalFactors <- scale(orthogonal_FinalFactors,TRUE,TRUE)
  
  
  # Compute factor hat 
  Factors_hat <- matrix(nrow = num_obs, ncol = 0)  # start with zero columns
  
  
  factor_index <- 1
 
  for (key in names(Final_list)){
    
    # extract combination
    combination <- as.numeric(unlist(strsplit(key, "-")))
    #Block <- do.call(cbind, lapply(combination, function(idx) data[, ranges[[idx]]]))
    Block <- do.call(cbind, lapply(combination, function(idx) Yorig[, ranges[[idx]]]))
    n_factors <- Final_list[[key]]
    
    
    # Extract corresponding factors
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
 
  Residuals <- Yorig - orthogonal_FinalFactors %*% Lambda
  
  
  # Store results
  results[["Factors"]] <- orthogonal_FinalFactors
  results[["Factors_hat"]] <- Factors_hat
  results[['Lambda']] <- t(Lambda)
  results[['Residuals']] <- Residuals
  results[['Factors_list']] <- Final_list
  
  
  
  #results[["Residuals"]] <- Yorig - orthogonal_FinalFactors %*% Lambda
  
 
  return(results)
}


