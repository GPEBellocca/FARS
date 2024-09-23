


StressFactors <- function(Factors, alpha = 0.95) {
  # Calculate the covariance matrix of the factors
  cov_matrix <- cov(Factors)
  
  # Calculate the mean of the factors
  factor_mean <- colMeans(Factors)
  
  # Compute the quantile for the stress condition
  stress_threshold <- qnorm(alpha)
  
  # Apply the stress condition to each factor (shifting by alpha quantile)
  stress_factors <- sweep(Factors, 2, factor_mean) + stress_threshold * sqrt(diag(cov_matrix))
  
  return(stress_factors)
}



