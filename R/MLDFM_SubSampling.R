# Multi-Level Dynamic Factor Model with subsampling procedure


MLDFM_SubSampling <- function(data, r = c(1), blocks = 1, block_ind = NULL, tol = 0.000001, max_iter = 1000, n_sample = 10, sample_size = 1) {
  
  n_factors <- sum(r)
  n_obs <- nrow(data)
  
  # Create output  structures
  sub_sampling_factors <- vector("list", n_sample)
  sub_sampling_factors_hat <- vector("list", n_sample)
  sub_sampling_lambda <- vector("list", n_sample)
  sub_sampling_residuals <- vector("list", n_sample)
  
 
  for (i in 1:n_sample) {
    
    # compute sample 
    sub_sample_result <- SubSample(data, block_ind = block_ind, n = n_blocks,r,sample_size)
    
    # call MLDFM for each sample
    MLDFM_result <- MLDFM(sub_sample_result$sample_data, r=r, blocks = n_blocks, block_ind = sub_sample_result$sample_block_ind, tol = 0.000001, max_iter = 1000) 
    
    sub_sampling_factors[[i]] <- MLDFM_result$Factors
    sub_sampling_factors_hat[[i]] <- MLDFM_result$Factors_hat
    sub_sampling_lambda[[i]] <- MLDFM_result$Lambda
    sub_sampling_residuals[[i]] <- MLDFM_result$Residuals
  }
  
  
  return(list(Factors_samples = sub_sampling_factors, Factors_hat_samples = sub_sampling_factors_hat,
         Lamba_samples = sub_sampling_lambda, Residuals_samples = sub_sampling_residuals))
}
  
