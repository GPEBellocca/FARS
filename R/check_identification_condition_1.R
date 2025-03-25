# check_identification_condition_1 (factors)

check_identification_condition_1 <- function(Factors, tol = 1e-8) {
 
  cov_matrix <- cov(Factors)
  
  I <- diag(ncol(Factors))
  
  # Compute max absolute deviation from identity
  deviation <- cov_matrix - I
  max_dev <- max(abs(deviation))
  
  # Result
  is_orthonormal <- max_dev < tol

  
  if (is_orthonormal) {
    cat("✅ Identification (i) satisfied? YES \n")
  } else {
    cat("❌ Identification (i) satisfied? NO \n")
  }
  
  
  return(is_orthonormal)
}
