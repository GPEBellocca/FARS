# Density3 with exp and tanh for stability
library(sn)
library(nloptr)

Density3 <- function(All_q_matrix, edge = 0.05, est_points = 512, random_samples = 5000) {
  
  # Prepare quantiles
  quintiles <- c(0.00, 0.25, 0.50, 0.75, 1)
  quintiles[1] <- quintiles[1] + edge 
  quintiles[5] <- quintiles[5] - edge 
  
  # Extract number of obs
  n_obs <- nrow(All_q_matrix)
  
  # Initialize variables
  density <- c() # Density array
  density_matrix <- matrix(NA, nrow = n_obs, ncol = est_points) # Density matrix
  distribution <- matrix(0, n_obs, random_samples) # Skew-t distribution 
  
  for (tt in 1:n_obs) {
    
    # Initial values
    iqn <- qnorm(0.75) - qnorm(0.25) # Interquartile range of standard normal distribution
    l0 <- All_q_matrix[tt, 3]  # Location
    s0 <- log(max(1, (All_q_matrix[tt, 4] - All_q_matrix[tt, 2]) / iqn)) # Log of scale to use with exp
    sh0 <- 0 # Shape, no transformation needed here
    
    # Initial parameters for optimization
    x0 <- c(l0, s0, sh0)
    
    # Bounds
    LB <- c(l0 - 10, -Inf, -Inf) 
    UB <- c(l0 + 20, Inf, Inf)
    
    # Objective func
    objective_fn <- function(x) {
      transformed_values <- qst(quintiles, xi = x[1], omega = exp(x[2]), alpha = tanh(x[3]))
      return(sum((as.numeric(All_q_matrix[tt, ]) - transformed_values)^2))
    }
    
    # Non-linear optimization using nloptr
    result <- nloptr(
      x0 = x0,
      eval_f = objective_fn,
      lb = LB,
      ub = UB,
      opts = list(
        algorithm = "NLOPT_LN_SBPLX", 
        xtol_rel = 1.0e-8
      )
    )
    
    # Extract optimized parameters and apply transformations
    xi_opt <- result$solution[1]
    omega_opt <- exp(result$solution[2]) # Use exp to ensure omega is positive
    alpha_opt <- tanh(result$solution[3]) # Use tanh to keep alpha bounded
    
    # Generate n random samples from skew-t distribution 
    skt <- rst(n = random_samples, xi = xi_opt, omega = omega_opt, alpha = alpha_opt)
    
    # Store samples
    distribution[tt, ] <- skt
    
    # Compute density of generated samples
    fit <- dst(seq(-30, 10, length.out = est_points), xi = xi_opt, omega = omega_opt, alpha = alpha_opt)
    
    # Store density
    density <- c(density, fit)
    density_matrix[tt, ] <- fit
  }
  
  return(list(density = density, density_matrix = density_matrix, distribution = distribution))
}
