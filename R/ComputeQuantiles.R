

source("./FARS/R/QReg.R")
ComputeQuantiles <- function(dep_variable, factors, scenario, h = 1,   edge = 0.05, min = TRUE) {
 
  
  # Prepare quantiles
  quantiles <- c(0.00, 0.25, 0.50, 0.75, 1)
  quantiles[1] <- quantiles[1]+edge # adjust left edge
  quantiles[5] <- quantiles[5]-edge # adjust right edge
  
  # Output structures
  Quantiles <- matrix(nrow = length(dep_variable), ncol = length(quantiles))
  Scenario_Quantiles <- matrix(nrow = length(dep_variable), ncol = length(quantiles))
  
  # Loop through each quantile and compute Qreg
  for (i in seq_along(quantiles)) {
    q <- quantiles[i]
    QReg_result <- QReg(dep_variable, factors = factors, scenario, h=h, QTAU = q, min = min)
    Quantiles[, i] <- QReg_result$Pred_q  
    Scenario_Quantiles[, i] <- QReg_result$Scenario_Pred_q  
  }
  
  return(list(Quantiles = Quantiles, Scenario_Quantiles = Scenario_Quantiles))
 
}





