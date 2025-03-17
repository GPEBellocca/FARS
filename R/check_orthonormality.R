check_orthonormality <- function(factors) {

  # Compute orthonormality: (F'F) / T
  cov_matrix <- (t(factors) %*% factors) 

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

# check_orthonormality <- function(factors, tol = 1e-6) {
#   
#   T <- nrow(factors)  # Number of observations (rows)
#   
#   # Compute (F'F) / T
#   cov_matrix <- (t(factors) %*% factors) / T
#   
#   print(cov_matrix)
#   
#   # Identity matrix for comparison
#   identity_matrix <- diag(ncol(factors))
#   
#   # Check if matrix is numerically close to identity
#   if (max(abs(cov_matrix - identity_matrix)) < tol) {
#     cat("\n✅ Factors are orthonormal (F'F/T ≈ I)\n")
#   } else {
#     cat("\n❌ Factors are NOT orthonormal.\n")
#   }
#   
#   # Return the computed matrix for inspection
#   return(cov_matrix)
# }
