


# Compute Factors hat
Compute_factors_hat<- function(Yorig, ranges, Final_list, Factor_list,Loadings_list) {
  
  num_obs <- nrow(Yorig)
  
  # Compute factor hat
  Factors_hat <- matrix(nrow = num_obs, ncol = 0)
  
  for (key in names(Final_list)){
    # extract combination
    combination <- as.numeric(unlist(strsplit(key, "-")))
    
    # Extract block data
    Block <- do.call(cbind, lapply(combination, function(idx) Yorig[, ranges[[idx]]]))
    
    # Extract Factors
    n_facts <- Final_list[[key]]
    Facts <- Factor_list[[key]]
    
    # Extract Loadings
    Loads <- Loadings_list[[key]]
    
    # Compute Factors Hat
    N <- ncol(Loads)
    Facts_hat<-(1/N)*Block%*%t(Loads)
    
    Factors_hat <- cbind(Factors_hat, Facts_hat)
    
  }
 
  return(Factors_hat)
}

