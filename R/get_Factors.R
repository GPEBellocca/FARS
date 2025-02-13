# Extract factors from a given level
get_Factors <- function(Factor_list, combination, level) {
  
  matching_values <- list()
  
  
  for (key in names(Factor_list)) {
    key_elements <- as.integer(strsplit(key, "-")[[1]])
    
    if (length(key_elements) == (level) && all(combination %in% key_elements)) {
      matching_values <- cbind(matching_values, Factor_list[[key]])
    }
  }
  
  matching_values_numeric <- as.numeric(matching_values)
  original_dims <- dim(matching_values)
  
  if (is.null(original_dims)){
    matching_values <- NULL
  }else{
    matching_values <- matrix(matching_values_numeric, nrow = original_dims[1], ncol = original_dims[2])
    
  }
  
  return(matching_values)
}

