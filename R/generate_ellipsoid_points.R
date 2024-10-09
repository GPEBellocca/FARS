library(MASS) 

generate_ellipsoid_points <- function(center, sigma_inv, calpha, n_points = 300) {
 

  n_dim <- length(center)
  
  # Generate random points from a multivariate normal distribution
  random_points <- mvrnorm(n_points, mu = rep(0, n_dim), Sigma = diag(n_dim))
  
  # Scale the random points to lie on the surface of the ellipsoid
  # sqrt(calpha) scales the points based on the size of the ellipsoid
  scaling_factors <- sqrt(calpha) / sqrt(rowSums(random_points^2))
  scaled_points <- random_points * scaling_factors
  
  # Apply the transformation to match the hyperellipsoid defined by sigma_inv
  ellipsoid_points <- t(solve(chol(sigma_inv))) %*% t(scaled_points)
  
  # Shift the points to the center of the ellipsoid
  ellipsoid_points <- t(ellipsoid_points) + matrix(center, nrow = n_points, ncol = n_dim, byrow = TRUE)
  
  return(ellipsoid_points)
}

