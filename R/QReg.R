# QReg


QReg <- function(dep_variable, factors, h=1,  QTAU=0.05, edge = 0.05) {
  
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
  
  # prepare quintiles list
  quintiles <- c(0.00, 0.25, 0.50, 0.75, 1)
  quintiles[1] <- quintiles[1]+edge # adjust left edge
  quintiles[5] <- quintiles[5]-edge # adjust right edge
  quintiles <- c(quintiles, QTAU) # add QTAU
  
  # final df
  All_q <- data.frame(matrix(ncol = 0, nrow = t))
  
  
  for (q in quintiles) {
    fit_q <- rq(formula, tau = q, data = reg_data)
    
    coefficients <- coef(fit_q)
    
    Pred_q <- as.numeric(coefficients[1])
    Pred_q <- Pred_q + as.numeric(coefficients[2]) * Y[] 
    for (i in 1:r) {
      Pred_q <- Pred_q + as.numeric(coefficients[i+2]) * factors[, i]
    }
    
    # save results
    column_name <- sprintf("Q.%02d", q * 100)  
    All_q[[column_name]] <- Pred_q
  }
  
  All_q <- as.data.frame(lapply(All_q, as.numeric))
  In_q <- head(All_q, (t-h)) # in sample
  Out_q <- tail(All_q, h) # out of sample
  
  return(list(All_q = All_q, In_q = In_q, Out_q = Out_q))

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
