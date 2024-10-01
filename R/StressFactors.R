StressFactors <- function(factors, quantile_models, tau_vec, alpha, target_tau) {
  # factors: Data frame or matrix of factors (n x p)
  # quantile_models: List of quantile regression models for each τ* in tau_vec
  # tau_vec: Vector of quantiles τ*
  # alpha: Level α defining the stress (e.g., 0.05)
  # target_tau: τ, the quantile at which to compute GiS (e.g., 0.05)
  
  # Load necessary libraries
  if (!requireNamespace("sn", quietly = TRUE)) {
    install.packages("sn")
  }
  library(sn)
  
  # Compute mean and covariance of factors
  mu_F <- colMeans(factors)
  Sigma_F <- cov(factors)
  p <- ncol(factors)
  
  # Compute c = sqrt(χ²ₚ(1 - α))
  c <- sqrt(qchisq(1 - alpha, df = p))
  
  # Initialize vector to store minimal growths
  minimal_growths <- numeric(length(tau_vec))
  
  for (i in seq_along(tau_vec)) {
    tau_star <- tau_vec[i]
    model <- quantile_models[[i]]
    
    # Extract coefficients
    beta <- as.numeric(coef(model))
    beta_0 <- beta[1]     # Intercept
    beta_factors <- beta[-1]  # Coefficients for factors
    
    # Compute sqrt(β' Σ_F β)
    beta_Sigma_beta <- t(beta_factors) %*% Sigma_F %*% beta_factors
    beta_Sigma_beta_sqrt <- sqrt(beta_Sigma_beta)
    
    # Compute minimal q_{τ*}
    q_tau_star <- beta_0 + sum(beta_factors * mu_F) + c * beta_Sigma_beta_sqrt
    
    minimal_growths[i] <- q_tau_star
  }
  
  # Fit a skewed-t distribution to the minimal growths
  fit_st <- selm(minimal_growths ~ 1, family = "ST")
  
  # Compute GiS at target_tau
  GiS <- qst(target_tau, dp = fit_st@param$dp)
  
  # Return GiS and additional results
  return(list(GiS = GiS, minimal_growths = minimal_growths, fit = fit_st))
}
