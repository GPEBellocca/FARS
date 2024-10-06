#sub_sample <- function(data, block_ind, n, sample_size) {
SubSample <- function(data, block_ind, n, sample_size) {
  new_data <- list()
  new_block_ind <- numeric(n)  # Will store the new block indices for the reduced data
  
  for (i in 1:n) {
    # Define the start and end indices for each block (columns)
    start_idx <- ifelse(i == 1, 1, block_ind[i-1] + 1)
    end_idx <- block_ind[i]
    
    # Extract the current block (subset of columns)
    block_data <- data[, start_idx:end_idx]
    
    # Randomly select 80% of the columns from the current block
    selected_cols <- sample(ncol(block_data), size = floor(0.8 * ncol(block_data)))
    
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