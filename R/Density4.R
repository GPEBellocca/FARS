# Density4 using nlm for optimization with parameter stabilization
library(sn)

Density4 <- function(All_q_matrix, edge = 0.05, est_points = 512, random_samples = 5000) {
  
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
    s0 <- max(1, (All_q_matrix[tt, 4] - All_q_matrix[tt, 2]) / iqn) # Scale, constrained to be >= 1
    sh0 <- 0 # Shape
    
    # Initial parameters for optimization
    x0 <- c(l0, log(s0), sh0) # Use log(s0) to work with exp()
    
    # Objective func
    objective_fn <- function(x) {
      xi <- x[1]
      omega <- exp(x[2]) # Ensure omega is positive
      alpha <- tanh(x[3]) # Constrain alpha to be within a reasonable range
      # Use solver="RFB" explicitly for qst()
      transformed_values <- qst(quintiles, xi = xi, omega = omega, alpha = alpha, solver = "RFB")
      return(sum((as.numeric(All_q_matrix[tt, ]) - transformed_values)^2))
    }
    
    # Non-linear optimization using nlm
    result <- nlm(objective_fn, p = x0)
    
    # Extract optimized parameters
    optimized_params <- result$estimate
    xi_opt <- optimized_params[1]
    omega_opt <- exp(optimized_params[2]) # Transform back
    alpha_opt <- tanh(optimized_params[3]) # Constrain alpha
    
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
