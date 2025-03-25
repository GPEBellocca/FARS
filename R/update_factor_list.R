# orthonormalize_factors <- function(factors) {
#   factors <- as.matrix(factors)
#   
#   qr_decomp <- qr(factors)
#   Q <- qr.Q(qr_decomp)
#   #Q <- Q * sqrt(nrow(factors))  # Scale to ensure F'F/T = I
#   return(Q)
# }
# 



# Update factor list
update_factor_list <- function(Factor_list, FinalFactors, r) {
  
  r_index <- 1
  counter <- 1
  
  filtered_r <- r[r != 0]
  # Update  Factors list
  for (key in names(Factor_list)) {
    
    
    select_factors <- FinalFactors[,counter:(counter+filtered_r[r_index]-1)]
    
    factor_matrix <- matrix(select_factors, ncol = filtered_r[r_index])
    #factor_matrix <- prewhiten_factors(factor_matrix)
    #check_orthonormality(factor_matrix)
   
    
    Factor_list[[key]] <- matrix(select_factors, ncol = filtered_r[r_index])
    Factor_list[[key]] <- factor_matrix
    
    counter <- counter + filtered_r[r_index]
    r_index <- r_index + 1
  }
  
  
  return(Factor_list)
}

