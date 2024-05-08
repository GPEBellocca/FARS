#Growth at Risk


GaR <- function(data, dep_variable, factors, QTAU=0.05) {
  
  
  r <- dim(factors)[2]
  t<-nrow(data)
  N<-ncol(data)
  
  
  # shift factors
  lag <- 1
  for (i in 1:2) {
    factors[,i] <- shift( factors[,i],lag)
  }
 
  
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
