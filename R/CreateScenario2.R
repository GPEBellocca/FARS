beta_ols <- function(X, Y) {
  solve(t(X) %*% X) %*% t(X) %*% Y
}


CreateScenario2 <- function(MLDFM_result, Sub_sampling_result,data,block_ind ,n_samples,alpha=0.95) {
  
  data <- scale(data,TRUE,TRUE)
  
  # Extract factors
  Factors <- MLDFM_result$Factors
  Factors_hat <- MLDFM_result$Factors_hat
  Factors_list <- MLDFM_result$Factors_list
  
  Factors_samples <- Sub_sampling_result$Factors_samples
  Factors_hat_samples <- Sub_sampling_result$Factors_hat_samples
  
  n_obs <- nrow(Factors)
  tot_n_factors <-  ncol(Factors)
  #keys <- names(Factors_list)
  #n_keys <- length(keys)
  
  # Define block ranges and count the number of var in each range
  ranges <- list()
  num_vars <- numeric(length(block_ind))  
  
  
  
  for (i in 1:length(block_ind)) {
    if (i == 1) {
      ranges[[i]] <- 1:block_ind[i]
    } else {
      ranges[[i]] <- (block_ind[i - 1] + 1):block_ind[i]
    }
    
    num_vars[i] <- length(ranges[[i]]) 
  }
  
  
  
  
  # Compute ellipsoid center for each obs
  CenterHE_matrix <- Factors
  
  
  
 
  
  
  # Initialize sigma 
  Sigma_list <- list()
  
  for (obs in 1:n_obs) {
    
    # initialize Sigma diagonal for current obs
    #Sigma <- matrix(0, nrow = n_factors, ncol = n_factors)
    Sigma_diag <- c()
    factor_index <- 1
    # loop over  factors
    for (key in names(Factors_list)){
      
      # extract combination
      combination <- as.numeric(unlist(strsplit(key, "-")))
      Block <- do.call(cbind, lapply(combination, function(idx) data[, ranges[[idx]]]))
      n_factors <- Factors_list[[key]]
    
      
      # Extract corresponding factors
      Facts <- Factors[,factor_index:(factor_index+n_factors-1)]
      
      
      # Compute Loadings and Residuals
      Loads <- beta_ols(Facts, Block)
      Resid <- Block - Facts %*% Loads
      
      # Compute Factors Hat
      N <- ncol(Loads)
      Facts_hat<-(1/N)*Block%*%t(Loads)
      
     
      
      # Compute number of variables
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
        # Extract sample's Factors
        Factors_hat_s <- Factors_hat_samples[[s]]
        Facts_hat_s <- Factors_hat_s[,factor_index:(factor_index+n_factors-1)]
        Facts_hat_s <- as.matrix(Facts_hat_s) # impose matrix structure 
        
        # Compute sample's Loadings 
        #Loads_s <- beta_ols(Facts_s, Block)
        
        # Compute sample's Factors Hat of the sample
        #N <- ncol(Loads_s)
        #Facts_hat_s <-(1/N)*Block%*%t(Loads_s)
        
        
        # Extract current obs
        Facts_hat_s_obs <- Facts_hat_s[obs,]
        Facts_hat_obs <- Facts_hat[obs,]
        
        # Compute diff
        diff <- Facts_hat_s_obs - Facts_hat_obs
        
        
        term2 <- term2 + (diff %*% t(diff)) 
        
      }
      
      
      
      Sig <- inv_Loads %*% ((term2/n_samples)+Gamma) %*% inv_Loads
     
      
      
      if(ncol(Sig)==1){
        Sigma_diag <- c(Sigma_diag,Sig)
      }else{
        Sigma_diag <- c(Sigma_diag,diag(Sig))
      }
      
      factor_index <- factor_index + n_factors
    
    }
   
    Sigma <- diag(Sigma_diag)
    #Sigma <- Sig
   
    
    
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
    
    if(tot_n_factors > 2){
      hellip <- hyperellipsoid(center_obs, solve(sigma_obs), calpha)
      Hyperellipsoids[[obs]] <- t(hypercube_mesh(8,hellip,TRUE))
    }else{
      Hyperellipsoids[[obs]] <- ellipse(sigma_obs, centre=center_obs, level = alpha, npoints = 300)
      #Hyperellipsoids[[obs]] <- generate_ellipsoid_points(center_obs, solve(sigma_obs), calpha, n_points=300)
    }

    
    
  }
  
  return(Hyperellipsoids)
}
