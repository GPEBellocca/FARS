# Density
library(sn)

Density <- function(All_q_matrix,  edge = 0.05, est_points = 512, random_samples = 5000) {
 
  
  # prepare quantiles
  quantiles <- c(0.00, 0.25, 0.50, 0.75, 1)
  quantiles[1] <- quantiles[1]+edge # adjust left edge
  quantiles[5] <- quantiles[5]-edge # adjust right edge
  
  
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
    
    
    LB = c(   -10+l0,     1,   -100) 
    UB = c(   +20+l0,    50,    100)
    
   
    skewt<-optim(c(l0, s0, sh0),fn=function(x){
      sum((as.numeric(All_q_matrix[tt,])-qst(quantiles,xi=x[1],omega=x[2],alpha=x[3]))^2)
    }, lower=LB,upper=UB,  method="L-BFGS-B")
    
    
    
    # generate n random sample from skew-t distribution 
    skt<-rst(n=random_samples, xi=skewt$par[1], omega=skewt$par[2], alpha= skewt$par[3], dp=NULL) 
   
    # store samples 
    distribution[tt,]<-skt
  
    # compute density of generated samples
    fit<-dst(seq(-30,10,length.out = est_points),xi=skewt$par[1],omega=skewt$par[2],alpha=skewt$par[3])
    
    # fit density
    density<-c(density,fit)
    density_matrix[tt, ] <- fit
    
  }
  
  
  return(list(density = density, density_matrix = density_matrix, distribution = distribution))
  
}





