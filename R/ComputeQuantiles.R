

ComputeQuantiles <- function(dep_variable, factors , h = 1, edge = 0.05, scenario = NULL, min = TRUE) {
 
  
  # Prepare quantiles
  quantiles <- c(0.00, 0.25, 0.50, 0.75, 1)
  quantiles[1] <- quantiles[1]+edge # adjust left edge
  quantiles[5] <- quantiles[5]-edge # adjust right edge
  
  # Output structures
  Quantiles <- matrix(nrow = length(dep_variable), ncol = length(quantiles))
  Scenario_Quantiles <- if (!is.null(scenario)) matrix(nrow = length(dep_variable), ncol = length(quantiles)) else NULL
  coeff_df <- NULL
  pvalue_df <- NULL
  
  # Loop through each quantile and compute Qreg
  for (i in seq_along(quantiles)) {
    q <- quantiles[i]

    QReg_result <- QReg(dep_variable, factors = factors, MLDFM_result$Factors_list, scenario = scenario, h=h, QTAU = q, min = min)

    if (!is.null(scenario)) {
      Scenario_Quantiles[, i] <- QReg_result$Scenario_Pred_q
    }

    Quantiles[, i] <- QReg_result$Pred_q
    coeff_df <- cbind(coeff_df, round(QReg_result$Coeff,3))
    pvalue_df  <- cbind(pvalue_df, round(QReg_result$Pvalue,2))
  }

  colnames(coeff_df) <- paste0("q", quantiles)
  colnames(pvalue_df) <- paste0("q", quantiles)
  
  
  if (!is.null(scenario)) {
    return(list(Quantiles = Quantiles, Coeff = coeff_df, Pvalue = pvalue_df,  Scenario_Quantiles = Scenario_Quantiles))
  }

  return(list(Quantiles = Quantiles, Coeff = coeff_df, Pvalue = pvalue_df))

}





