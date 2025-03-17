prewhiten_factors <- function(Factors) {
  T <- nrow(Factors)  # Number of observations
  k <- ncol(Factors)  # Number of factors
  
  # Compute covariance matrix (divided by T for correct normalization)
  Sigma_F <- (t(Factors) %*% Factors) / T
  
  if (k == 1) {
    # Special case: only one factor (1x1 covariance matrix)
    Sigma_F_inv_sqrt <- 1 / sqrt(Sigma_F)  # Just invert the scalar
  } else {
    # Standard case: multiple factors
    eig_decomp <- eigen(Sigma_F)
    V <- eig_decomp$vectors  # Eigenvectors
    D <- diag(eig_decomp$values)  # Eigenvalues
    
    # Compute inverse square root of Sigma_F
    D_inv_sqrt <- diag(1 / sqrt(diag(D)))
    Sigma_F_inv_sqrt <- V %*% D_inv_sqrt %*% t(V)
  }
  
  # Transform the factors
  Adjusted_Factors <- Factors %*% Sigma_F_inv_sqrt
  
  # FINAL NORMALIZATION: Ensure unit variance
  Adjusted_Factors <- Adjusted_Factors / sqrt(T)
  
  return(Adjusted_Factors)
}
