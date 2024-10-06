# Density
library(sn)

Density <- function(dep_variable, Factors, h = 1,   edge = 0.05, est_points = 512, random_samples = 5000) {
 
  
  # prepare quantiles
  quintiles <- c(0.00, 0.25, 0.50, 0.75, 1)
  quintiles[1] <- quintiles[1]+edge # adjust left edge
  quintiles[5] <- quintiles[5]-edge # adjust right edge
  
  
  All_q_matrix <- matrix(nrow = length(dep_variable), ncol = length(quintiles))
  
  # Loop through each quantile and comoute Qreg
  for (i in seq_along(quintiles)) {
    q <- quintiles[i]
    Pred_q <- QReg(dep_variable, Factors, h=h, QTAU = q)
    All_q_matrix[, i] <- Pred_q  
  }
  
  # extract number of observations
  n_obs = nrow(All_q_matrix) 
 
  # Initialize variables
  density<-c() # density array
  density_matrix <- matrix(NA, nrow = n_obs, ncol = est_points) # density matrix
  distribution=matrix(0,n_obs,random_samples) # skew-t distribution 
  
  for (tt in 1:n_obs){
    
    # Initial values
    iqn=qnorm(0.75)-qnorm(0.25) # Interquartile range of standard normal distr
    l0=All_q_matrix[tt,3]  # Location
    s0=(All_q_matrix[tt,4] - All_q_matrix[tt,2]) / iqn # Scale
    sh0=0 # Shape
    
    # Lower and Upper bounds
    LB = c(   -10+l0,     1,   -100) #Omega must positive >=1
    UB = c(   +20+l0,    50,    100)
    
    # optim minimize the squared differences between the q of the obs and the q of a skew-t distribution
    # qst compute the q of skew-t distribution
    skewt<-optim(c(l0, s0, sh0),fn=function(x){
      sum((as.numeric(All_q_matrix[tt,])-qst(quintiles,xi=x[1],omega=x[2],alpha=x[3]))^2)
    }, lower=LB,upper=UB,  method="L-BFGS-B")
    
    # ????
    #skt_q<-qst(quintiles,xi=skewt$par[1],omega=skewt$par[2],alpha=skewt$par[3])
    
    # generate n random sample from skew-t distribution 
    skt<-rst(n=random_samples, xi=skewt$par[1], omega=skewt$par[2], alpha= skewt$par[3], dp=NULL) #delet nu=4
   
    # store samples 
    distribution[tt,]<-skt
  
    # compute density of generated samples
    fit<-dst(seq(-30,10,length.out = est_points),xi=skewt$par[1],omega=skewt$par[2],alpha=skewt$par[3])
    
    #fit<-density(skt,from=-30,to=10,nseq=est_points)
    density<-c(density,fit)
    density_matrix[tt, ] <- fit
    
  }
  
  
  return(list(density = density, density_matrix = density_matrix, distribution = distribution))
  
}





