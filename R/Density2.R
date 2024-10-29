# Density
library(sn)
#library(geoR)

Density2 <- function(All_q_matrix,  edge = 0.05, est_points = 512, random_samples = 5000) {
 
  
  # prepare quantiles
  quintiles <- c(0.00, 0.25, 0.50, 0.75, 1)
  quintiles[1] <- quintiles[1]+edge # adjust left edge
  quintiles[5] <- quintiles[5]-edge # adjust right edge
  
  
  # extract number of observations
  n_obs = nrow(All_q_matrix) 
 
  # Return variables 
  density<-c() # density array
  density_matrix <- matrix(NA, nrow = n_obs, ncol = est_points) # density matrix
  distribution=matrix(0,n_obs,random_samples) # skew-t distribution 
  
  for (tt in 1:n_obs){
    
    # nlm procedure
    
    
    # minimization function


    miniGaR<-function(xx){((All_q_matrix[tt,5]-qst(quintiles[5], xi=xx[1], omega=max(xx[2], 0.0001), alpha=xx[3]))+
                             (All_q_matrix[tt,4]-qst(quintiles[4], xi=xx[1], omega=max(xx[2], 0.0001), alpha=xx[3]))+
                             (All_q_matrix[tt,3]-qst(quintiles[3], xi=xx[1], omega=max(xx[2], 0.0001), alpha=xx[3]))+
                             (All_q_matrix[tt,2]-qst(quintiles[2], xi=xx[1], omega=max(xx[2], 0.0001), alpha=xx[3]))+
                             (All_q_matrix[tt,1]-qst(quintiles[1], xi=xx[1], omega=max(xx[2], 0.0001), alpha=xx[3])))^2}

# 
#   
#     
#     miniGaR<-function(xx){((All_q_matrix[tt,5]-qst(quintiles[5], xi=xx[1], omega=abs(xx[2]), alpha=xx[3]))+
#                              (All_q_matrix[tt,4]-qst(quintiles[4], xi=xx[1], omega=abs(xx[2]), alpha=xx[3]))+
#                              (All_q_matrix[tt,3]-qst(quintiles[3], xi=xx[1], omega=abs(xx[2]), alpha=xx[3]))+
#                              (All_q_matrix[tt,2]-qst(quintiles[2], xi=xx[1], omega=abs(xx[2]), alpha=xx[3]))+
#                              (All_q_matrix[tt,1]-qst(quintiles[1], xi=xx[1], omega=abs(xx[2]), alpha=xx[3])))^2}


# 
#       miniGaR<-function(xx){((All_q_matrix[tt,5]-qst(quintiles[5], xi=xx[1], omega=(xx[2]), alpha=xx[3]))+
#                          (All_q_matrix[tt,4]-qst(quintiles[4], xi=xx[1], omega=(xx[2]), alpha=xx[3]))+
#                          (All_q_matrix[tt,3]-qst(quintiles[3], xi=xx[1], omega=(xx[2]), alpha=xx[3]))+
#                          (All_q_matrix[tt,2]-qst(quintiles[2], xi=xx[1], omega=(xx[2]), alpha=xx[3]))+
#                          (All_q_matrix[tt,1]-qst(quintiles[1], xi=xx[1], omega=(xx[2]), alpha=xx[3])))^2}
# 

    
    # parameters
    iqn=qnorm(quintiles[4])-qnorm(quintiles[2]) # interquartile range 1.34898
    
    xi <- mean(All_q_matrix[,3])
    omega <- -(mean(All_q_matrix[,4])-mean(All_q_matrix[,2])/ iqn)
    #omega <- -(mean(All_q_matrix[,4])-mean(All_q_matrix[,2]))/ iqn
    alpha <- 1.5
    p <- c(xi,omega,alpha)
    
    
    # model fit
    fit4<-nlm(miniGaR,p=p)
    
    # extract result form nlm
    xi <- fit4$estimate[1]
    omega <-fit4$estimate[2]
    alpha <- fit4$estimate[3]
    
    
    # compute distribution
    skt<-rst(n=random_samples, xi=xi, omega=omega, alpha=alpha, nu=4, dp=NULL)
    
    # store distribution
    distribution[tt,]<-skt
    
    # compute density
    fit<-density(skt,from=-40,to=10,n=est_points)
    
    # store density
    density<-c(density,fit$y)
    density_matrix[tt, ] <- fit$y
    
  }
  
  
  return(list(density = density, density_matrix = density_matrix, distribution = distribution))
  
}





