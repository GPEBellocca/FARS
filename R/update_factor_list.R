# Update factor list
update_factor_list <- function(Factor_list, FinalFactors, r) {
  
  r_index <- 1
  counter <- 1
  
  filtered_r <- r[r != 0]
  # Update  Factors list
  for (key in names(Factor_list)) {
    
    select_factors <- FinalFactors[,counter:(counter+filtered_r[r_index]-1)]
    
    Factor_list[[key]] <- matrix(select_factors, ncol = filtered_r[r_index])
    
    counter <- counter + filtered_r[r_index]
    r_index <- r_index + 1
  }
  
  
  return(Factor_list)
}

