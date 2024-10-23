CreateScenario <- function(MLDFM_result, Sub_sampling_result,n_factors, n_samples,alpha=0.95) {
  
  
  Factors <- MLDFM_result$Factors
  Factors_hat <- MLDFM_result$Factors_hat
  Loadings <- MLDFM_result$Loadings
  Residuals <- MLDFM_result$Residuals
  
  
  
  
  #first_key <- names(Factors)[1]
 
  #n_obs <- ifelse(!is.null(nrow(Factors[[first_key]])), nrow(Factors[[first_key]]), length(Factors[[first_key]]))
  n_obs <- nrow(Factors[[1]])
  n_keys <- length(Factors)
  
  
  
  # Compute ellipsoid center for each obs
  #CenterHE_list <- list()
  CenterHE_matrix <- matrix(NA, nrow = n_obs, ncol = n_factors)
  
  
  
  
  for (obs in 1:n_obs) {
    center <- c()
    
    for(k in 1:n_keys){
      factors_key <- Factors[[k]]
      if (!is.matrix(factors_key)) {
        factors_key <- matrix(unlist(factors_key), ncol = 1)  # Convert list (or vector) to a one-column matrix
      }
      center <- c(center,factors_key[obs,])
      
    }
    
    
    #CenterHE_list[[obs]] <- center
    CenterHE_matrix[obs, ] <- center
  }
  
 
 
  
  
  
  
  # Initialize sigma 
  Sigma_list <- list()
  
  for (obs in 1:n_obs) {
    
    # initialize Sigma diagonal for current obs
    #Sigma <- matrix(0, nrow = n_factors, ncol = n_factors)
    Sigma_diag <- c()
    # loop over each factor
    for (k in 1:n_keys){
      
      # Extract Factors, Loadings and Residuals
      Facts_hat <- Factors_hat[[k]]
      
      
      Loads <- Loadings[[k]]
      Resid <- Residuals[[k]]
      n_var <- ncol(Loads)
      
      
      
      # Precompute the inverse of the loading matrix
      inv_Loads <- solve((Loads %*% t(Loads)) / n_var)
      
      Gamma <- matrix(0, nrow = ncol(Facts_hat), ncol = ncol(Facts_hat))
      # Compute Gamma
      for(v in 1:n_var){
        term <- (Loads[,v] %*% t(Loads[,v]))*(Resid[obs,v]^2)
        Gamma <- Gamma + term
      }
      
      Gamma <- Gamma / n_var
      
      
      
     
      
      term2 <- matrix(0, nrow = ncol(Facts_hat), ncol = ncol(Facts_hat))
      for(s in 1:n_samples){
        Factors_hat_s <- Sub_sampling_result$Factors_hat_s[[s]]
        Facts_hat_s <- Factors_hat_s[[k]]
        Facts_hat_s_obs <- Facts_hat_s[obs,]
        Facts_hat_obs <- Facts_hat[obs,]
        
        diff <- Facts_hat_s_obs - Facts_hat_obs
        # print(Facts_hat_s_obs)
        # print(Facts_hat_obs)
        # print(diff)
        
       
        
        term2 <- term2 + (diff %*% t(diff)) 
        
      }
      
      
      
      Sig <- inv_Loads %*% ((term2/n_samples)+Gamma) %*% inv_Loads
      
      
      
      if(ncol(Sig)==1){
        Sigma_diag <- c(Sigma_diag,Sig)
      }else{
        Sigma_diag <- c(Sigma_diag,diag(Sig))
      }
      
    
    }
   
    Sigma <- diag(Sigma_diag)
    Sigma <- Sig
    
    
    # Store the sigma matrix for this observation
    Sigma_list[[obs]] <- Sigma
  
    
  }
  
  
  
  # Hyperellipsoids 
  Hyperellipsoids <- list()
  
  # Loop over each observation and compute the hyperellipsoid
  for (obs in 1:n_obs) {
    #center_obs <- CenterHE_list[[obs]]  # Center 
    center_obs <- CenterHE_matrix[obs,]
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
