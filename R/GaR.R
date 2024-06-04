# Growth at Risk


GaR <- function(dep_variable, factors, QTAU=0.05) {
  
  t<- dim(factors)[1]
  r <- dim(factors)[2]
  lag <- 1
  
  # prepare regression data
  Y <- dep_variable
  LagY<-shift(Y,lag)
  
  shifted_factors <- matrix(NA, nrow = nrow(factors), ncol = r)
  for (i in 1:r) {
    shifted_factors[,i] <- shift( factors[,i],lag)
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
  
  
  ################################################
  ########  QUANTILE QTAU  ##########
  ################################################
  
  # run Qregression
  GaRfit <- rq(formula, tau = QTAU, data = reg_data)
  
  # prediction 
  coefficients <- coef(GaRfit)
  PredGaR <- as.numeric(coefficients[1])
  PredGaR <- PredGaR + as.numeric(coefficients[2]) * Y[]
  for (i in 1:r) {
    PredGaR <- PredGaR + as.numeric(coefficients[i+2]) * factors[, i]
  }
  
  #print(PredGaR)
  #print(PredGaR[dim(reg_data)[1]])
  
  # Box test
  rho <- function(u,tau=QTAU)u*(tau - (u < 0))
  Bfit<-rq(Y~1 ,tau = QTAU) # baseline model -> Solution may be non unique
  X_squared <- Box.test(GaRfit$residuals, lag = 4, type = c("Box-Pierce", "Ljung-Box"), fitdf = 0)
  V1 <- sum(rho(GaRfit$resid, GaRfit$tau))
  V2 <- sum(rho(Bfit$resid, Bfit$tau))
  R1<-1-V1/V2
  
  
  ################################################
  ######## ALL QUINTILES ##########
  ################################################
  
  quintiles <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  AllQuintiles_df <- data.frame(Y = reg_data$Y)  
  
  for (q in quintiles) {
    GaRfit_q <- rq(formula, tau = q, data = reg_data)
    
    coefficients <- coef(GaRfit_q)
    
    PredGaR_q <- as.numeric(coefficients[1])
    PredGaR_q <- PredGaR_q + as.numeric(coefficients[2]) * LagY[] 
    for (i in 1:r) {
      PredGaR_q <- PredGaR_q + as.numeric(coefficients[i+2]) * shifted_factors[, i]
    }
    
    # save results
    column_name <- sprintf("GaR.%02d", q * 100)  
    AllQuintiles_df[[column_name]] <- PredGaR_q
  }
  
  
  ################################################
  ############ DENSITY GaR   #####################
  ################################################
  
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
