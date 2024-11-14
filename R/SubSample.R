

SubSample <- function(data, block_ind, n,r,sample_size = 1) {
  new_data <- list()
  new_block_ind <- numeric(n)  #block indices for the sample data
  
  for (i in 1:n) {
    # Define the start and end indices for each block 
    start_idx <- ifelse(i == 1, 1, block_ind[i-1] + 1)
    end_idx <- block_ind[i]
    
    # Extract the current block 
    block_data <- data[, start_idx:end_idx]
    
    # Randomly select the columns from the current block
    n_vars <- ncol(block_data)
    n_obs <- nrow(block_data)
    if (sample_size>=1){
    portion=1-235/((n_vars+25)^2)-0.2*(sqrt(n_vars)/n_obs)-(1/(n_obs*n_vars))
    }else{
      portion <- sample_size
    }
    
    #set.seed(123)
    selected_cols <- sample(n_vars, size = round(n_vars*portion,0)
                            ,replace =FALSE ,prob = NULL)
    
    # Store the reduced block data
    new_data[[i]] <- block_data[, selected_cols]
    
    # Update the new block index
    if (i == 1) {
      new_block_ind[i] <- ncol(new_data[[i]])
    } else {
      new_block_ind[i] <- new_block_ind[i-1] + ncol(new_data[[i]])
    }
  }
  
  # Combine all the new block data into a single dataframe
  final_data <- do.call(cbind, new_data)
  
  return(list(sample_data = final_data, sample_block_ind = new_block_ind))
}