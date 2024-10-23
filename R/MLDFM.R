# Multi-Level Dynamic Factor Model

MLDFM <- function(data, r = c(1), blocks = 1, block_ind = NULL, tol = 0.000001, max_iter = 1000) {
  
  
  
  ##### FACTOR EXTRACTION #####
  
  if (blocks==1){
    result <- SingleBlock(data,r=r)
  }else if(blocks>1){
    #result <- MultipleBlocks(data, r=r,block_ind = block_ind, tol = tol, max_iter = max_iter)
    result <- MultipleBlocks(data, r=r,block_ind = block_ind, tol = tol, max_iter = max_iter)
    
  }else{
    print('Error - Invalid number of block')
  }
  
  
  return(list(Data = data, 
              Factors = result$Factors, 
              Factors_hat = result$Factors_hat, 
              Lambda = result$Lambda,  
              Loadings = result$Loadings, 
              Residuals = result$Residuals ))
}


