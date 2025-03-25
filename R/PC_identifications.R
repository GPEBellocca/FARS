#library(GPArotation)
#library(MASS)




#loadings_signs <- c("1-2" = 1, "1" = -1, "2" = -1) #hcpi
#loadings_signs <- c("1-2" = 1, "1" = 1, "2" = 1) #ccpi
#loadings_signs <- c("1-2" = -1, "1" = -1, "2" = -1) #fcpi
loadings_signs <- c("1-2" = -1, "1" = -1, "2" = 1) #ecpi

# Impose identifications
PC_identifications <- function(Yorig, num_blocks, ranges, num_factors, r, currentFactors, Factor_list, Loadings_list) {
  
  T_obs <- nrow(Yorig)
  FinalFactors <- matrix(nrow = T_obs, ncol = 0)  
  
  # orthogonalization of all factors (identification 3)
  orthogonal_Factors <- orthogonalize_factors(currentFactors) 
  
  
  Factor_list <- update_factor_list(Factor_list, orthogonal_Factors, r)
  L_res <- Compute_Lambda(Yorig,num_blocks,ranges,num_factors,r,Factor_list,Loadings_list)
  Loadings_list <- L_res$Loadings_list
  
  # Initialize Lambda
  Lambda <- matrix(0, nrow = num_factors, ncol =  ncol(Yorig))
  
  r_index <- 1
  counter <- 1
  
  # Extract Global Factors and Loadings
  key <- paste(seq(1, num_blocks), collapse = "-")
  GlobalFactors <- Factor_list[[key]]
  GlobalLoadings <- Loadings_list[[key]]
  
  

  CommonComponent <- GlobalFactors %*% GlobalLoadings
  pca_result <- prcomp(CommonComponent, center = FALSE, scale. = FALSE)
  GlobalFactors_new <- pca_result$x[, 1:ncol(GlobalFactors), drop = FALSE]  
  GlobalFactors_new<- scale(GlobalFactors_new,TRUE,TRUE)
  
  GlobalFactors_new <- loadings_signs["1-2"]*GlobalFactors_new #adjust sign
  
  GlobalLoadings_new <- qr.solve(GlobalFactors_new, CommonComponent)

  # Update lists
  FinalFactors <- cbind(FinalFactors, GlobalFactors_new)
  Factor_list[[key]] <- GlobalFactors_new
  Loadings_list[[key]] <- GlobalLoadings_new
  
  # Update Lambda
  combination <- seq(1, num_blocks)
  Lambda[counter:(counter+r[r_index]-1), unlist(ranges[combination])] <- GlobalLoadings_new
  counter <- counter + r[r_index]
  
  # Check identifications
  # print(key)
  # check_identification_condition_1(GlobalFactors_new)
  # check_identification_condition_2(t(GlobalLoadings_new))
  
  
  
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
      Loadings <- Loadings_list[[key]]
      
      
      CommonComponent <- Factors %*% Loadings
      pca_result <- prcomp(CommonComponent, center = FALSE, scale. = FALSE)
      Factors_new <- pca_result$x[, 1:ncol(Factors), drop = FALSE]  # 59 x r
      Factors_new<- scale(Factors_new,TRUE,TRUE)
      
      
      Factors_new <- loadings_signs[key]*Factors_new #adjust sign
      
      Loadings_new <- qr.solve(Factors_new, CommonComponent)

      # Update lists
      FinalFactors <- cbind(FinalFactors, Factors_new)
      Factor_list[[key]] <- Factors_new
      Loadings_list[[key]] <- Loadings_new
      
      # Update Lambda
      Lambda[counter:(counter+r[r_index]-1), unlist(ranges[combination])] <- Loadings_new
      counter <- counter + r[r_index]
      
      # Check identifications
      # print(key)
      # check_identification_condition_1(Factors_new)
      # check_identification_condition_2(t(Loadings_new))
      
    }
  }
  results <- list()  
  results[["Factor_list"]] <- Factor_list
  results[["FinalFactors"]] <- FinalFactors
  results[["Lambda"]] <- Lambda
  return(results)
}

