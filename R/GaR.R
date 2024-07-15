# GaR


GaR <- function(distribution, QTAU = 0.05) {
 
  
  n_obs = nrow(distribution) 
  QGaR <- array(dim = n_obs)
  
  for (tt in 1:n_obs){
    QGaR[tt]=quantile(distribution[tt,],QTAU)
  }
  
  return(QGaR)
  
}



