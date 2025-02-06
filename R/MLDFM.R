# Multi-Level Dynamic Factor Model





MLDFM <- function(data, r = c(1), blocks = 1, block_ind = NULL, tol = 0.000001, max_iter = 1000, method = 0) {
  
  
  ##### FACTOR EXTRACTION #####
  
  if (blocks==1){
    # one block
    result <- SingleBlock(data,r=r)
    
    return(list(Factors = result$Factors,
                Factors_hat = result$Factors_hat,
                Lambda = result$Lambda,
                Residuals = result$Residuals,
                Factors_list = result$Factors_list
                
    )) 
    
  }else if(blocks>1){
    # multiple blocks
    #result <- MultipleBlocks(data, r=r,block_ind = block_ind, tol = tol, max_iter = max_iter)
    result <- MultipleBlocks(data, r=r,block_ind = block_ind, tol = tol, 
                                  max_iter = max_iter, method = method)
    
    return(list(Initial_Factors = result$Initial_Factors,
                Factors = result$Factors,
                Factors_hat = result$Factors_hat,
                Lambda = result$Lambda,
                Residuals = result$Residuals,
                Factors_list = result$Factors_list,
                RSS_list = result$RSS_list
    )) 
    
  }else{
    # invalid
    print('Error - Invalid number of block')
    return()
  }
  
 
  # return(list(Initial_Factors = result$Initial_Factors,
  #             Factors = result$Factors,
  #             Factors_hat = result$Factors_hat,
  #             Lambda = result$Lambda,
  #             Residuals = result$Residuals,
  #             Factors_list = result$Factors_list,
  #             RSS_list = result$RSS_list
  #             )) 
}


