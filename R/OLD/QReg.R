# QReg


QReg <- function(dep_variable, factors, scenario=scenario, h=1,  QTAU=0.05, min = TRUE) {
  
  t<- dim(factors)[1]
  r <- dim(factors)[2]
  
  
  # prepare regression data
  Y <- dep_variable
  LagY<-shift(Y,h)
  
  shifted_factors <- matrix(NA, nrow = nrow(factors), ncol = r)
  for (i in 1:r) {
    shifted_factors[,i] <- shift( factors[,i],h)
  }

  
  reg_data <- data.frame(Y = Y, LagY = LagY)
  reg_data <- cbind(reg_data, shifted_factors)  
  
  
  names(reg_data)[1:2] <- c("Y", "LagY")
  new_factor_names <- paste("factor", 1:r, sep = "")  
  names(reg_data)[3:(2 + r)] <- new_factor_names  
  factor_names_concat <- paste(new_factor_names, collapse = " + ")
  
  # build formula
  formula_str <- paste("Y ~ LagY", factor_names_concat, sep = " + ")
  formula <- as.formula(formula_str)

  
  # qreg
  fit_q <- rq(formula, tau = QTAU, data = reg_data)
    
  coefficients <- coef(fit_q)
    
  Pred_q <- as.numeric(coefficients[1])
  Pred_q <- Pred_q + as.numeric(coefficients[2]) * Y[] 
  for (i in 1:r) {
    Pred_q <- Pred_q + as.numeric(coefficients[i+2]) * factors[, i]
  }
    
  
  # qreg scenario
  Scenario_Pred_q <- vector(mode = "numeric", length = t)
  
  for (tt in 1:t){
    ellips <- scenario[[tt]]
    points <- nrow(ellips)
    pred <- vector(mode = "numeric", length = points)
    
    for(p in 1:points){
      pred[p] <- as.numeric(coefficients[1])
      pred[p] <- pred[p] + as.numeric(coefficients[2]) * Y[tt] 
      
      for (j in 1:r) {
        pred[p] <- pred[p] + as.numeric(coefficients[j+2]) * ellips[p,j]
      }
    }
    
    
    if(min == TRUE){
      value <- min(pred)
    }else{
      value <- max(pred)
    }
    Scenario_Pred_q[tt] <- value
  }
    
  # return
  return(list(Pred_q = Pred_q, Scenario_Pred_q = Scenario_Pred_q))

}


# shift function
shift <- function(x, lag) {
  n <- length(x)
  xnew <- rep(NA, n)
  if (lag < 0) {
    xnew[1:(n-abs(lag))] <- x[(abs(lag)+1):n]
  } else if (lag > 0) {
    xnew[(lag+1):n] <- x[1:(n-lag)]
  } else {
    xnew <- x
  }
  return(xnew)
}
