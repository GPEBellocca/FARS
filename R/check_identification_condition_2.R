# check_identification_condition_2 (loadings)
check_identification_condition_2 <- function(loadings) {
  cov_matrix <- t(loadings) %*% loadings  
  
  # Check if off-diagonal elements are effectively zero
  off_diag <- cov_matrix[upper.tri(cov_matrix) | lower.tri(cov_matrix)]
  is_diagonal <- all(abs(off_diag) < 1e-6)
  
  if (is_diagonal) {
    cat("✅ Identification (ii)satisfied?  YES \n")
  } else {
    cat("❌ Identification (ii) satisfied?  NO\n")
  }
  
  return(is_diagonal)
}
