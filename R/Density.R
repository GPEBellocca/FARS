# Density
library(sn)

Density <- function(All_q, est_points = 512, random_samples = 5000) {
 
  
  # drop QTAU
  All_q <- All_q[,1:5]
  
  # extract quintiles list
  quintiles <- as.numeric(sub("Q.", "", colnames(All_q)))/100
  
  # convert quintiles data frame to matrix
  All_q_matrix <- as.matrix(All_q)
  
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
    s0=(All_q_matrix[tt,4] - All_q_matrix[tt,2]) / iqn; # Scale
    sh0=0 # Shape
    
    # Lower and Upper bounds
    LB = c(   -10+l0,     1,   -100); #Omega must positive >=1
    UB = c(   +20+l0,    50,    100);
    
    # optim minimize the squared differences between the q of the obs and the q of a skew-t distribution
    # qst compute the q of skew-t distribution
    skewt<-optim(c(l0, s0, sh0),fn=function(x){
      sum((as.numeric(All_q_matrix[tt,])-qst(quintiles,xi=x[1],omega=x[2],alpha=x[3]))^2)
    }, lower=LB,upper=UB,  method="L-BFGS-B")
    
    # ????
    #skt_q<-qst(quintiles,xi=skewt$par[1],omega=skewt$par[2],alpha=skewt$par[3])
    
    # generate n random sample from skew-t distribution 
    skt<-rst(n=random_samples, xi=skewt$par[1], omega=skewt$par[2], alpha= skewt$par[3], nu=4, dp=NULL)
   
    # store samples 
    distribution[tt,]<-skt
  
    # compute density of generated samples
    # dst(x,xi=skewt$par[1],omega=skewt$par[2],alpha=skewt$par[3])
    fit<-density(skt,from=-30,to=10,nseq=512)
    density<-c(density,fit$y)
    density_matrix[tt, ] <- fit$y
    
  }
  
  
  return(list(density = density_matrix, distribution = distribution))
  
}



