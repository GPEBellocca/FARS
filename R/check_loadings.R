# Function to check if LAMBDA'LAMBDA is diagonal
check_loadings <- function(loadings) {
  cov_matrix <- t(loadings) %*% loadings  # Compute covariance
  
  # Check if off-diagonal elements are effectively zero
  off_diag <- cov_matrix[upper.tri(cov_matrix) | lower.tri(cov_matrix)]
  is_diagonal <- all(abs(off_diag) < 1e-6)
  
  if (is_diagonal) {
    cat("✅ LAMBDA'LAMBDA is diagonal.\n")
  } else {
    cat("❌ LAMBDA'LAMBDA is NOT diagonal.\n")
  }
  
  return(is_diagonal)
}
