#' Extract Factors at a Given Hierarchical Level
#'
#' @keywords internal
#' 

get_factors <- function(Factor_list, combination, level) {
  
  matching_values <- list()
  
  
  for (key in names(Factor_list)) {
    key_elements <- as.integer(strsplit(key, "-")[[1]])
    
    # Check if this key corresponds to the requested level and includes all in combination
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


# get_factors <- function(Factor_list, combination, level) {
#   
#   matched <- list()
#   
#   for (key in names(Factor_list)) {
#     key_elements <- as.integer(strsplit(key, "-")[[1]])
#     
#     # Check if this key corresponds to the requested level and includes all in combination
#     if (length(key_elements) == level && all(combination %in% key_elements)) {
#       matched[[length(matched) + 1]] <- Factor_list[[key]]
#     }
#   }
#   
#   # Combine matched factors if any found
#   if (length(matched) > 0) {
#     return(do.call(cbind, matched))
#   } else {
#     return(NULL)
#   }
# }
