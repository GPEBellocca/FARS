

orthonormalize_factors <- function(factors) {
  factors <- as.matrix(factors)
  qr_decomp <- qr(factors)
  Q <- qr.Q(qr_decomp)
  Q <- Q * sqrt(nrow(factors))  # Scale to ensure F'F/T = I
  return(Q)
}


check_orthonormality <- function(factors) {
 
  # Compute orthonormality: (F'F) / T
  #orthonormality_matrix <- t(factors) %*% factors / nrow(factors)
  cov_matrix <- t(factors) %*% factors
  
  # Check if matrix is identity
  identity_matrix <- diag(ncol(factors))
  if (all(round(cov_matrix, 4) == identity_matrix)) {
    cat("\n✅ Factors are orthonormal (F'F/T = I)\n")
  } else {
    cat("\n❌ Factors are NOT orthonormal.\n")
  }
  
  # Return the matrix
  return(cov_matrix)
}

# Function to check if LAMBDA'LAMBDA is diagonal
check_loadings_diagonal <- function(loadings) {
  cov_matrix <- t(loadings) %*% loadings  # Compute covariance
  
  # Check if off-diagonal elements are effectively zero
  off_diag <- cov_matrix[upper.tri(cov_matrix) | lower.tri(cov_matrix)]
  is_diagonal <- all(abs(off_diag) < 1e-6)
  
  # Print results
  cat("\nLAMBDA'LAMBDA Covariance Matrix:\n")
  
  
  if (is_diagonal) {
    cat("✅ LAMBDA'LAMBDA is diagonal.\n")
  } else {
    cat("❌ LAMBDA'LAMBDA is NOT diagonal.\n")
  }
  
  return(is_diagonal)
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
    GlobalFactors <- blockfact0(Yorig, num_vars, number_of_factor, rep(1, length(num_vars)))
  }else{
    # PCA 
    pca_result <- prcomp(Yorig, scale. = FALSE)
    GlobalFactors <- pca_result$x[, 1:number_of_factor]
    GlobalFactors <- GlobalFactors / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(GlobalFactors) %*% GlobalFactors))))
  }
  #GlobalFactors <- orthonormalize_factors(GlobalFactors)
  
  check_orthonormality(GlobalFactors)
  
  
  
  GlobalLoadings <- beta_ols(GlobalFactors, Yorig)
  #GlobalLoadings <- pca_result$rotation[, 1:number_of_factor]
  
  
  check_loadings_diagonal(t(GlobalLoadings))
  
  
  
  
 
  #print(check_factor_orthonormality(GlobalFactors))
  
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
        Factors <- blockfact0(Residuals, num_vars[combination], number_of_factor, rep(1, length(combination)))
      }else{
        # Use PCA
        pca_result <- prcomp(Residuals, scale. = FALSE)
        Factors <- pca_result$x[, 1:number_of_factor]
        Factors <- Factors / kronecker(matrix(1, nrow = num_obs, ncol = 1), t(sqrt(diag(t(Factors) %*% Factors))))
        
      }
      
      #Factors <- orthonormalize_factors(Factors)
      check_orthonormality(Factors)
      Loadings <- beta_ols(Factors, Residuals)
      check_loadings_diagonal(t(Loadings))
      
      # Store Factors
      key <- paste(combination, collapse = "-")
      Factor_list[[key]] <- Factors  
      InitialFactors <- cbind(InitialFactors, Factors)
    }
    
  }
  
  return(1)
  
  results <- list()
  results[["InitialFactors"]] <- InitialFactors
  results[["Factor_list"]] <- Factor_list
  
  return(results)
}

