


ComputeQuantiles <- function(dep_variable, Factors, scenario = NULL, h = 1,   edge = 0.05, min = TRUE) {
 
  
  # prepare quantiles
  quintiles <- c(0.00, 0.25, 0.50, 0.75, 1)
  quintiles[1] <- quintiles[1]+edge # adjust left edge
  quintiles[5] <- quintiles[5]-edge # adjust right edge
  
  
  Quantiles <- matrix(nrow = length(dep_variable), ncol = length(quintiles))
  Scenario_Quantiles <- matrix(nrow = length(dep_variable), ncol = length(quintiles))
  
  # Loop through each quantile and compute Qreg
  for (i in seq_along(quintiles)) {
    q <- quintiles[i]
    QReg_result <- QReg(dep_variable, Factors, scenario, h=h, QTAU = q, min = min)
    Quantiles[, i] <- QReg_result$Pred_q  
    Scenario_Quantiles[, i] <- QReg_result$Scenario_Pred_q  
  }
  
  return(list(Quantiles = Quantiles, Scenario_Quantiles = Scenario_Quantiles))
 
}





