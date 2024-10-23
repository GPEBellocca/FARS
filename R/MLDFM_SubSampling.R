# Multi-Level Dynamic Factor Model with subsampling procedure

MLDFM_SubSampling <- function(data, r = c(1), blocks = 1, block_ind = NULL, tol = 0.000001, max_iter = 1000, samples = 10) {
  
  n_factors <- sum(r)
  n_obs <- nrow(data)
  
  # Create output factors structure 
  sub_sampling_factors <- vector("list", samples)
  
  # Create output factors hat structure 
  sub_sampling_factors_hat <- vector("list", samples)
  
  # Create output loadings structure 
  sub_sampling_lambda <- vector("list", samples)
  
  # Create output residuals structure 
  sub_sampling_residuals <- vector("list", samples)
  
 
  # Extract factors from each sample
  for (i in 1:samples) {
    
    
    sub_sample_result <- SubSample(data, block_ind = block_ind, n = n_blocks,r)
    
    MLDFM_result <- MLDFM(sub_sample_result$sample_data, r=r, blocks = n_blocks, block_ind = sub_sample_result$sample_block_ind, tol = 0.000001, max_iter = 1000) 
    
    sub_sampling_factors[[i]] <- MLDFM_result$Factors
    sub_sampling_factors_hat[[i]] <- MLDFM_result$Factors_hat
    sub_sampling_lambda[[i]] <- MLDFM_result$Lambda
    
    #sub_sampling_residuals[[i]] <- MLDFM_result$Residuals
  }
  
  
  return(list(Factors_s= sub_sampling_factors, Factors_hat_s= sub_sampling_factors_hat, Lambda_s = sub_sampling_lambda))
}
  
