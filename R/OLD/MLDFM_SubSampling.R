# Multi-Level Dynamic Factor Model with subsampling procedure

MLDFM_SubSampling <- function(data, outlier = TRUE, r = c(1), blocks = 1, block_ind = NULL, tol = 0.000001, max_iter = 1000, samples = 10, sample_size = 0.8) {
  
  n_factors <- sum(r)
  n_obs <- nrow(data)
  
  # Create output factors structure 
  sub_sampling_factors <- vector("list", samples)
  
  # Create output loadings structure 
  sub_sampling_loadings <- vector("list", samples)
  
  # Create output residuals structure 
  sub_sampling_residuals <- vector("list", samples)
  
 
  
  # Extract factors from each sample
  for (i in 1:samples) {
    sub_sample_result <- SubSample(data, block_ind = block_ind, n = n_blocks, sample_size = sample_size)
    MLDFM_result <- MLDFM(sub_sample_result$sample_data, outlier = FALSE, r=r, blocks = n_blocks, block_ind = sub_sample_result$sample_block_ind, tol = 0.000001, max_iter = 1000) 
   
    sub_sampling_factors[[i]] <- MLDFM_result$Factors
    sub_sampling_loadings[[i]] <- MLDFM_result$Loadings
    sub_sampling_residuals[[i]] <- MLDFM_result$Residuals
  }
  
  
  return(list(Factors_s= sub_sampling_factors, Loadings_s = sub_sampling_loadings, Residuals_s =  sub_sampling_residuals))
}
  
