CreateScenario <- function(Factors, MLDFM_SubSampling,n_samples,alpha=0.95) {
  
  n_obs <- nrow(Factors)
  n_factors <- ncol(Factors)
  
  
  # Compute ellipsoid center for each obs
  CenterHE_list <- list()
  
  for (obs in 1:n_obs) {
    CenterHE_list[[obs]] <- Factors[obs, 1:n_factors]
  }
  
  # Initialize sigma 
  Sigma_list <- list()
  
  # Compute sigma (covariance matrix) based on the subsampling matrices for each observation
  for (obs in 1:n_obs) {
    Sigma <- matrix(0, nrow = n_factors, ncol = n_factors)
    
    
    for (i in 1:n_factors) {
      
      subsample_i <- MLDFM_result_sub_sampling[[i]][obs, 1:n_samples]
      
      # Variance for the diagonal elements
      Sigma[i, i] <- var(subsample_i)
      
     
      if (i < n_factors) {  
        for (j in (i+1):n_factors) {
          subsample_j <- MLDFM_result_sub_sampling[[j]][obs, 1:n_samples]
          
          # Compute the covariance between the two factors for the off diagonal elements
          Sigma[i, j] <- cov(subsample_i, subsample_j)
          Sigma[j, i] <- Sigma[i, j]  
        }
      }
    }
    
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
      Hyperellipsoids[[obs]] <- ellipse(sigma_obs, centre=center_obs, level = alpha,npoints = 300)
      #Hyperellipsoids[[obs]] <- generate_ellipsoid_points(center_obs, solve(sigma_obs), calpha, n_points=300)
    }

    
    
  }
  
  return(Hyperellipsoids)
}
