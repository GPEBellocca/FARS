normalize_factors <- function(Factors) {
  # Ensure Factors is a matrix
  Factors <- as.matrix(Factors)  
  k <- ncol(Factors)  # Number of factors
  
  # Compute standard deviation of each factor
  factor_sd <- sqrt(diag(t(Factors) %*% Factors)) 
  
  if (k == 1) {
    # Special case: Only one factor (avoid using sweep)
    Factors <- Factors / factor_sd  # Direct scalar division
  } else {
    # Multiple factors case: Use sweep for column-wise division
    Factors <- sweep(Factors, 2, factor_sd, "/")
  }
  
  return(Factors)
}





# Compute Initial Factors
Compute_Initial_Factors <- function(Yorig, num_vars, num_obs, num_blocks, ranges, num_factors, r, method) {
  
  # Define Factors data structures
  Factor_list <- list()
  InitialFactors <- matrix(nrow = num_obs, ncol = 0)  
  
  
  
  # Compute Global factors
  r_index <- 1 
  number_of_factor <- r[r_index] # number of factor to be extracted with PCA
  
  if (method == 0){
    # CCA
    GlobalFactors <- blockfact0(Yorig, num_vars, number_of_factor, rep(1, num_blocks))
    #GlobalFactors <- blockfact_cca(Yorig, num_vars, number_of_factor, rep(1, num_blocks))
    
    
    
  }else{
    # PCA 
    pca_result <- prcomp(Yorig, scale. = FALSE)
    GlobalFactors <- pca_result$x[, 1:number_of_factor]
    #GlobalFactors <- GlobalFactors / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(GlobalFactors) %*% GlobalFactors))))
  }
  
  
  
  GlobalFactors <- normalize_factors(GlobalFactors)
  #check_orthonormality(GlobalFactors)
  GlobalLoadings <- beta_ols(GlobalFactors, Yorig)
  #check_loadings(t(GlobalLoadings))
  
 
  
  
  # Store Global factors
  key <- paste(seq(1, num_blocks), collapse = "-")  
  Factor_list[[key]] <- GlobalFactors  
  InitialFactors <- cbind(InitialFactors, GlobalFactors)
    
  
  
  
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
      
      
      # Compute factors
      number_of_factor <- r[r_index] # number of factor to be extracted with PCA
      if (i < num_blocks - 1 && method == 0) {
        # Use CCA for middle level
        Factors <- blockfact0(Residuals, num_vars[combination], number_of_factor, rep(1, num_blocks))
        #Factors <- blockfact_cca(Residuals, num_vars[combination], number_of_factor, rep(1, num_blocks))
       
      }else{
        # Use PCA
        pca_result <- prcomp(Residuals, scale. = FALSE)
        Factors <- pca_result$x[, 1:number_of_factor]
      }
      
      
      
      Factors <- normalize_factors(Factors)
      #check_orthonormality(Factors)
      Loadings <- beta_ols(Factors, Residuals)
      #check_loadings(t(Loadings))
      
      # Store Factors
      key <- paste(combination, collapse = "-")
      Factor_list[[key]] <- Factors  
      InitialFactors <- cbind(InitialFactors, Factors)
    }
    
  }
  
  
  
  results <- list()
  results[["InitialFactors"]] <- InitialFactors
  results[["Factor_list"]] <- Factor_list
  
  return(results)
}

