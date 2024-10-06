# Multi-Level Dynamic Factor Model with subsampling procedure

MLDFM_SubSampling <- function(data, outlier = TRUE, r = c(1), blocks = 1, block_ind = NULL, tol = 0.000001, max_iter = 1000, samples = 10, sample_size = 0.8) {
  
  n_factors <- sum(r)
  n_obs <- nrow(data)
  
  # Create output structure 
  sub_sampling_factors <- vector("list", n_factors)
  
  for (i in 1:n_factors) {
    sub_sampling_factors[[i]] <- matrix(0, nrow = n_obs, ncol = samples)  
  }
  
  # Extract factors from each sample
  for (i in 1:samples) {
    sub_sample_result <- SubSample(data, block_ind = block_ind, n = n_blocks, sample_size = sample_size)
    MLDFM_result <- MLDFM(sub_sample_result$sample_data, outlier = FALSE, r=r, blocks = n_blocks, block_ind = sub_sample_result$sample_block_ind, tol = 0.000001, max_iter = 1000) 
    
    for (j in 1:n_factors) {
      samples <- sub_sampling_factors[[j]]
      samples[,i] <- MLDFM_result$Factors[,j]
      sub_sampling_factors[[j]] <- samples
      
    }
    
  }
  
  return(sub_sampling_factors)
}
  
