# GaR


GaR <- function(distribution, QTAU=0.05) {
 
  
 
  # compute extreme quantiles of the skew-t distribution
  #left <- quintiles[1]
  #right <- quintiles[5]
  
  #QGaR.left <- array(dim = n_obs)
  #QGaR.right <- array(dim = n_obs)
  
  #for (tt in 1:n_obs){
    #QGaR.left[tt]=quantile(distirbution[,tt],left)
    #QGaR.right[tt]=quantile(distirbution[,tt],right)
  #}
  
  # compute extreme quantiles of the skew-t distribution
  
  
  n_obs = nrow(distribution) 
  
  QGaR <- array(dim = n_obs)
  
  for (tt in 1:n_obs){
  QGaR[tt]=quantile(distribution[tt],QTAU)
  }
  
  return(QGaR)
  
}



