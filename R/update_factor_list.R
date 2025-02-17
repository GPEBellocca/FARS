orthonormalize_factors <- function(factors) {
  factors <- as.matrix(factors)
  
  qr_decomp <- qr(factors)
  Q <- qr.Q(qr_decomp)
  #Q <- Q * sqrt(nrow(factors))  # Scale to ensure F'F/T = I
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

# Update factor list
update_factor_list <- function(Factor_list, FinalFactors, r) {
  
  r_index <- 1
  counter <- 1
  
  filtered_r <- r[r != 0]
  # Update  Factors list
  for (key in names(Factor_list)) {
    
    select_factors <- FinalFactors[,counter:(counter+filtered_r[r_index]-1)]
    
    factor_matrix <- matrix(select_factors, ncol = filtered_r[r_index])
    #factor_matrix <- factor_matrix / kronecker(matrix(1, nrow = nrow(factor_matrix), ncol = 1), t(sqrt(diag(t(factor_matrix) %*% factor_matrix))))
    
    factor_matrix <- orthonormalize_factors(factor_matrix)
    check_orthonormality(factor_matrix)
    
    #Factor_list[[key]] <- matrix(select_factors, ncol = filtered_r[r_index])
    Factor_list[[key]] <- factor_matrix
    
    counter <- counter + filtered_r[r_index]
    r_index <- r_index + 1
  }
  
  
  return(Factor_list)
}

