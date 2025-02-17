# check_factor_orthonormality <- function(factors) {
#   T <- nrow(factors)
#   cov_matrix <- t(factors) %*% factors / T
#   cat("\nFactor Covariance (F'F/T):\n")
#   print(round(cov_matrix, 4))
#   
#   diag_ok <- all(abs(diag(cov_matrix) - 1) < 1e-6)
#   offdiag_ok <- all(abs(cov_matrix[upper.tri(cov_matrix)]) < 1e-6)
#   
#   if (diag_ok && offdiag_ok) {
#     cat("✅ Factors are properly orthonormalized (F'F/T = I).\n")
#   } else {
#     cat("❌ Factors are NOT orthonormalized.\n")
#   }
#   
#   return(diag_ok && offdiag_ok)
# }


# Function to check F'F/T = I
check_factor_orthonormality <- function(factors) {
  T <- nrow(factors)
  cov_matrix <- t(factors) %*% factors / T
  #print("Factor Covariance (F'F/T):")
  #print(round(cov_matrix, 4))
  return(all(abs(diag(cov_matrix) - 1) < 1e-6) && all(abs(cov_matrix[upper.tri(cov_matrix)]) < 1e-6))
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



# Compute Lambda
Compute_Lambda <- function(Yorig, num_blocks, ranges, num_factors, r, Factor_list) {
  
  # Initialize Lambda
  Lambda <- matrix(0, nrow = num_factors, ncol =  ncol(Yorig))
  
  r_index <- 1
  counter <- 1
  
  # Extract Global Factors
  key <- paste(seq(1, num_blocks), collapse = "-")
  GlobalFactors <- Factor_list[[key]]
  
  
  
  # Compute Global Loadings
  GlobalLoadings <- beta_ols(GlobalFactors, Yorig)
  #GlobalLoadings <- t(Yorig) %*% GlobalFactors / nrow(GlobalFactors)
  
  
  
  check_loadings_diagonal(t(GlobalLoadings))
  
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
      
      # Extract current block factors
      key <- paste(combination, collapse = "-")
      Factors <- Factor_list[[key]]
      
      
      
      # Compute Loadings
      Loadings <- beta_ols(Factors, Residuals)
      #Loadings <- t(Residuals) %*% Factors / nrow(Factors)
      
      
      check_loadings_diagonal(t(Loadings))
      
      
      # Update Lambda
      Lambda[counter:(counter+r[r_index]-1), unlist(ranges[combination])] <- Loadings
      counter <- counter + r[r_index]
      
    }
  }
  return(Lambda)
}

