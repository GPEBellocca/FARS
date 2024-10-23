CreateScenario <- function(MLDFM_result, Sub_sampling_result,n_samples,alpha=0.95) {
  
  
  Factors <- MLDFM_result$Factors
  Loadings <- MLDFM_result$Loadings
  Residuals <- MLDFM_result$Residuals
  
  n_obs <- nrow(Factors)
  n_factors <- ncol(Factors)
  n_var <- ncol(Loadings)
  
  # Compute ellipsoid center for each obs
  CenterHE_list <- list()
  
  
  for (obs in 1:n_obs) {
    CenterHE_list[[obs]] <- Factors[obs, 1:n_factors]
  }
  
  # Initialize sigma 
  Sigma_list <- list()
  
  # Precompute the inverse of the loading matrix
  inv_Loadings <- solve((Loadings %*% t(Loadings)) / n_var)
  
  # Compute sigma (covariance matrix) based on the subsampling 
  for (obs in 1:n_obs) {
    Gamma <- matrix(0, nrow = n_factors, ncol = n_factors)
    
    for(v in 1:n_var){
      term <- (Loadings[,v] %*% t(Loadings[,v]))*(Residuals[obs,v]^2)
      Gamma <- Gamma + term
    }
    
    Gamma <- Gamma / n_var
    
    term2 <- matrix(0, nrow = n_factors, ncol = n_factors)
    for(s in 1:n_samples){
      Factors_s <- Sub_sampling_result$Factors_s[[s]]
      Factors_s_obs <- Factors_s[obs,]
      Factors_obs <- Factors[obs,]
      
      diff <- Factors_s_obs - Factors_obs
      term2 <- term2 + (diff %*% t(diff)) 
      
    }
    
    Sigma <- inv_Loadings %*% ((term2/n_samples)+Gamma) %*% inv_Loadings
    
    # diag_matrix <- Sigma
    # diag_matrix[upper.tri(diag_matrix) | lower.tri(diag_matrix)] <- 0
    
    # Store the sigma matrix for this observation
    Sigma_list[[obs]] <- Sigma
  }
  
  
  # Hyperellipsoids 
  Hyperellipsoids <- list()
  
  # Loop over each observation and compute the hyperellipsoid
  for (obs in 1:n_obs) {
    center_obs <- CenterHE_list[[obs]]  # Center 
    sigma_obs <- Sigma_list[[obs]]      # Covariance matrix 
    
    calpha <- sizeparam_normal_distn(alpha, d=n_factors)  # Size parameter for the hyperellipsoid
    #calpha <- qchisq(alpha_stress, df = n_factors) #quantile of chi-square distribution
    
    if(n_factors > 2){
      hellip <- hyperellipsoid(center_obs, solve(sigma_obs), calpha)
      Hyperellipsoids[[obs]] <- hypercube_mesh(8,hellip,TRUE)
    }else{
      Hyperellipsoids[[obs]] <- ellipse(sigma_obs, centre=center_obs, level = alpha, npoints = 300)
      #Hyperellipsoids[[obs]] <- generate_ellipsoid_points(center_obs, solve(sigma_obs), calpha, n_points=300)
    }

    
    
  }
  
  return(Hyperellipsoids)
}
