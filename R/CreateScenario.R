beta_ols <- function(X, Y) {
  solve(t(X) %*% X) %*% t(X) %*% Y
}


CreateScenario <- function(MLDFM_result, Sub_sampling_result,data,block_ind ,n_samples,alpha=0.95) {
  
  data <- scale(data,TRUE,TRUE)
  
  # Extract factors
  Factors <- MLDFM_result$Factors
  Factors_hat <- MLDFM_result$Factors_hat
  Factors_list <- MLDFM_result$Factors_list
  
  Factors_samples <- Sub_sampling_result$Factors_samples
  Factors_hat_samples <- Sub_sampling_result$Factors_hat_samples
  
  n_obs <- nrow(Factors)
  tot_n_factors <-  ncol(Factors)
  
  
  # Define block ranges and count the number of variables in each range
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
  
  
  
  
  # Set ellipsoid center for each obs
  CenterHE_matrix <- Factors

  
  # Initialize sigma 
  Sigma_list <- list()
  
  for (obs in 1:n_obs) {
    
    # initialize Sigma diagonal for current obs
    
    Sigma_diag <- c()
    factor_index <- 1
    
    Sigma <- matrix(0, nrow = tot_n_factors, ncol = tot_n_factors)
    
    # loop over  factors
    for (key in names(Factors_list)){
      
      # Extract data
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
      
      # Compute the inverse of the loading matrix
      inv_Loads <- solve((Loads %*% t(Loads)) / n_var)
      
      # Compute Gamma
      Gamma <- matrix(0, nrow = ncol(Facts_hat), ncol = ncol(Facts_hat))
      for(v in 1:n_var){
        term <- (Loads[,v] %*% t(Loads[,v]))*(Resid[obs,v]^2)
        Gamma <- Gamma + term
      }
      
      Gamma <- Gamma / n_var
      
      
      # Compute Sigma
      term2 <- matrix(0, nrow = ncol(Facts_hat), ncol = ncol(Facts_hat))
      for(s in 1:n_samples){
        # Extract sample's Factors
        Factors_hat_s <- Factors_hat_samples[[s]]
        Facts_hat_s <- Factors_hat_s[,factor_index:(factor_index+n_factors-1)]
        Facts_hat_s <- as.matrix(Facts_hat_s) # impose matrix structure 
        
        # Extract current obs
        Facts_hat_s_obs <- Facts_hat_s[obs,]
        Facts_hat_obs <- Facts_hat[obs,]
        
        # Compute diff
        diff <- Facts_hat_s_obs - Facts_hat_obs
        term2 <- term2 + (diff %*% t(diff)) 
        
      }
      
      Sig <- inv_Loads %*% ((term2/n_samples)+Gamma) %*% inv_Loads
      size <- ncol(Sig)
      
      Sigma[factor_index:(factor_index + size - 1), factor_index:(factor_index + size - 1)] <- Sig
      
      
      factor_index <- factor_index + n_factors
    
    }
   
    #Sigma <- diag(Sigma_diag)
    #Sigma <- Sig
   
    # Store the sigma matrix for this observation
    Sigma_list[[obs]] <- Sigma
  
   
  }
  
  # Hyperellipsoids 
  Hyperellipsoids <- list()
  
  # Loop over each observation and compute the hyperellipsoid
  for (obs in 1:n_obs) {
    
    center_obs <- CenterHE_matrix[obs,]
    sigma_obs <- Sigma_list[[obs]]     
    
    calpha <- sizeparam_normal_distn(alpha, d=n_factors)  # Size parameter 
   
    
    if(tot_n_factors > 2){
      # more than 2 dimensions
      hellip <- hyperellipsoid(center_obs, solve(sigma_obs), calpha)
      Hyperellipsoids[[obs]] <- t(hypercube_mesh(8,hellip,TRUE))
    }else{
      # 2D ellips
      Hyperellipsoids[[obs]] <- ellipse(sigma_obs, centre=center_obs, level = alpha, npoints = 300)
    }
    
  }
  
  return(Hyperellipsoids)
}
