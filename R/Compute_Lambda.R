# Compute Lambda
Compute_Lambda <- function(Yorig, num_blocks, ranges, num_factors, r, Factor_list) {
  
  # Initialize Lambda
  Lambda <- matrix(0, nrow = num_factors, ncol =  ncol(Yorig))
  
  r_index <- 1
  counter <- 1
  
  # Extract Global Factors
  key <- paste(seq(1, num_blocks), collapse = "-")
  GlobalFactors <- Factor_list[[key]]
  
  
  # Compute Global Loadings
  GlobalLoadings <- beta_ols(GlobalFactors, Yorig)
  
  # Update Lambda
  combination <- seq(1, num_blocks)
  Lambda[counter:(counter+r[r_index]-1), unlist(ranges[combination])] <- GlobalLoadings
  counter <- counter + r[r_index]
  
  
  # Loop on lower levels
  for (i in 1:(num_blocks-1)) {
    k <-  num_blocks - i
    combinations_matrix <- t(combn(num_blocks,k))
    for (j in 1:nrow(combinations_matrix)) {
      combination <- combinations_matrix[j,]
      
      r_index <- r_index + 1
      
      # Skip blocks where Factors are not needed
      if (r[r_index] == 0){
        next
      }
      
      
      # Extract Residuals filtering out upper levels factors
      level <- num_blocks
      
      Residuals <- do.call(cbind, lapply(combination, function(idx) Yorig[, ranges[[idx]]]))
      
      
      while (level > length(combination)) {
        Factors <- get_Factors(Factor_list, combination, level)
        
        # filter out
        if(!is.null(Factors)){
          ols_result <- beta_ols(Factors, Residuals)
          Residuals <- Residuals - Factors %*% ols_result
        }
        
        level <- level - 1
        
      }
      
      # Extract current block factors
      key <- paste(combination, collapse = "-")
      Factors <- Factor_list[[key]]
      
      
      # Compute Loadings
      Loadings <- beta_ols(Factors, Residuals)
      
      
      # Update Lambda
      Lambda[counter:(counter+r[r_index]-1), unlist(ranges[combination])] <- Loadings
      counter <- counter + r[r_index]
      
    }
  }
  return(Lambda)
}

