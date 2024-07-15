# Multi-Level Dynamic Factor Model - 3 blocks


beta_ols<-function(X,Y){
  t(Y) %*% X %*% solve(t(X) %*% X)
}

MutiLevelFactor3Blocks<-function(Yorig,r,block_ind,tol){
  #r_glob=c(1),r_reg=c(1,1,1)
  #Check function
  # Yorig <- openxlsx::read.xlsx("DataBase.xlsx",sheet = "fulldata",cols =2:625)
  # r1<-1:63
  # r2<-64:311
  # r3<-312:519
  # tol=0.000001
  
  
  r1 <- 1:block_ind[1]
  r2 <- (block_ind[1]+1):block_ind[2]
  r3 <- (block_ind[2]+1):block_ind[3]
  
  Yorig <- scale(Yorig) # Standarize
  R<-list() # where everything is saved
  names_Yorig<-colnames(Yorig)       # Get names
  
  
  ##################################################################
  ##                    Get groups for level 2                    ##
  ##################################################################
  

  R1 = Yorig[,r1] 
  NR1 = ncol(R1)
  R2 = Yorig[,r2] 
  NR2 = ncol(R2)
  R3 = Yorig[,r3] 
  NR3 = ncol(R3)

  
  ###########################################################################
  ###########################################################################
  ###                                                                     ###
  ###                 OBTAIN ESTIMATES OF INITIAL FACTORS                 ###
  ###                                                                     ###
  ###########################################################################
  ###########################################################################
  
  ##################################################################
  ##                    Global Factor                    ##
  ##################################################################
  
  # Order data and inputs for CCA
  y = cbind(R1,R2,R3)
  Nregio=c(NR1,NR2,NR3)
  sum(c(NR1,NR2,NR3))
  G0 = blockfact0(y, c(NR1,NR2,NR3),1,c(1,1,1))
  R[["GO"]]=scale(G0);  #CCA GLOBAL FACTOR
  t = nrow(y)
  yG0 = y-G0 %*% t(beta_ols(G0,y)) #Resid Global
  colnames(y)
  
  ##################################################################
  ##                    Semipervasive                    ##
  ##################################################################
  
  y_r12=yG0[,1:(NR1+NR2)] 
  y_r13=cbind(yG0[,1:NR1],yG0[,(NR1+NR2+1):(NR1+NR2+NR3)])
  y_r23=yG0[,(NR1+1):(NR1+NR2+NR3)]

  #F12 = blockfact0(y=y_r12,Nregio=c(NR1,NR2),r_glob=1,r_reg=c(1,1))  
  F13 = blockfact0(y_r13,c(NR1,NR3),1,c(1,1))  
  #F23 = blockfact0(y_r23,c(NR2,NR3),1,c(1,1))  
  
  y_r1=R1 
  y_r2=R2 
  y_r3=R3 
  
  lamR1 = beta_ols(F13,y_r1) 
  uR1L2 = y_r1-F13 %*% t(beta_ols(F13,y_r1))
  
  # lamR2 = beta_ols(cbind(F12,F23),y_r2) 
  # uR2L2 = y_r2-cbind(F12,F23) %*% t(beta_ols(cbind(F12,F23),y_r2))
  
  uR2L2 = y_r2
  
  lamR3 = beta_ols(F13,y_r3) 
  uR3L2 = y_r3-F13 %*% t(beta_ols(F13,y_r3))
  
  # lamR3 = beta_ols(cbind(F13,F23),y_r3) 
  # uR3L2 = y_r3-cbind(F13,F23) %*% t(beta_ols(cbind(F13,F23),y_r3))

  u=cbind(uR1L2,uR2L2,uR3L2)
  # 
  Nregio=c(NR1,NR2,NR3)
  r_reg=c(1,1,1)
  g = length(Nregio)
  RegInd = c(0, cumsum(Nregio))
  fhatreg = c() 
  for(i in 1:g){
    evec<-eigrs2(t(u[,(RegInd[i]+1):RegInd[i+1]]) %*% u[,(RegInd[i]+1):RegInd[i+1]])$evec
    fhatreg<-cbind(fhatreg,u[,(RegInd[i]+1):RegInd[i+1]]%*%evec[,(ncol(evec)-r_reg[i]+1):ncol(evec)])
  }
  
  fhatreg<-fhatreg / kronecker(matrix(1,nrow = nrow(fhatreg),ncol=1),t(matrix(sqrt(diag(t(fhatreg)%*%fhatreg)))))
  
  # Initial values  for while
  uu1 = sum(diag(t(u) %*%  u))
  uu0 = 10000000
  conteo = 0
  
  while ((log(uu0)-log(uu1))> tol){
    conteo = conteo+1
    lamG = beta_ols(G0,y)
    uB = y-G0 %*% t(beta_ols(G0,y))
    
    
    #ur12=uB[,1:(NR1+NR2)]
    ur13=cbind(uB[,1:NR1],uB[,(NR1+NR2+1):(NR1+NR2+NR3)]) 
    #ur23=uB[,(NR1+1):(NR1+NR2+NR3)]
    
    
    # lam12=beta_ols(F12,ur13)
    lam13=beta_ols(F13,ur13)
    # lam23=beta_ols(F23,ur13)
    
    
    lamR1 = beta_ols(F13,uB[,1:NR1]) 
    ur1 = uB[,1:NR1]-F13 %*% t(beta_ols(F13,uB[,1:NR1]))
    
    # lamR2 = beta_ols(cbind(F12,F23),uB[,(NR1+1):(NR1+NR2)]) 
    # ur2 = uB[,(NR1+1):(NR1+NR2)]-cbind(F12,F23) %*% t(beta_ols(cbind(F12,F23),uB[,(NR1+1):(NR1+NR2)]))
    
    
    ur2 = uB[,(NR1+1):(NR1+NR2)]
    
    # lamR3 = beta_ols(cbind(F13,F23),uB[,(NR1+NR2+1):(NR1+NR2+NR3)]) 
    # ur3 = uB[,(NR1+NR2+1):(NR1+NR2+NR3)]-cbind(F13,F23) %*% t(beta_ols(cbind(F13,F23),uB[,(NR1+NR2+1):(NR1+NR2+NR3)]))
    
    lamR3 = beta_ols(F13,uB[,(NR1+NR2+1):(NR1+NR2+NR3)]) 
    ur3 = uB[,(NR1+NR2+1):(NR1+NR2+NR3)]-F13 %*% t(beta_ols(F13,uB[,(NR1+NR2+1):(NR1+NR2+NR3)]))

    lamr1 = beta_ols(fhatreg[,1],ur1) 
    lamr2 = beta_ols(fhatreg[,2],ur2) 
    lamr3 = beta_ols(fhatreg[,3],ur3)
    
    lam = cbind(lamG,c(lam13[1:NR1],rep(0,NR2),lam13[(NR1+1):(NR1+NR3)]),
                #c(lam12[1:(NR1+NR2)],rep(0,NR3)),
                #c(rep(0,NR1),lam23[1:(NR2+NR3)]),
                c(lamr1,rep(0,NR2+NR3)),
                c(rep(0,NR1),lamr2,rep(0,NR3)),
                c(rep(0,NR1+NR2),lamr3))
    
    #################################################################
    ##                      Orthogonalization                      ##
    #################################################################
    
    fhat = beta_ols(lam,t(y))
    G0 = as.matrix(fhat[,1])
    F13 = as.matrix(fhat[,2])
    #F12 = as.matrix(fhat[,3])
    #F23 = as.matrix(fhat[,4])
    fhatreg = fhat[,3:5]
    
    F13 = F13-G0 %*% t(beta_ols(G0,F13))
    
    
    # R1x = fhatreg[,1]-c(F12,F13) %*% t(beta_ols(c(F12,F13),fhatreg[,1]))
    R1x = fhatreg[,1]-F13 %*% t(beta_ols(F13,fhatreg[,1]))
    fhatreg[,1] = R1x
    # R2x = fhatreg[,2]-c(F12,F23) %*% t(beta_ols(c(F12,F23),fhatreg[,2]))
    #fhatreg[,2] = R2x
    fhatreg[,2] = fhatreg[,2]
    # R3x = fhatreg[,3]-cbind(F13,F23) %*% t(beta_ols(cbind(F13,F23),fhatreg[,3]))
    R3x = fhatreg[,3]-F13 %*% t(beta_ols(F13,fhatreg[,3]))
    fhatreg[,3] = R3x
    
    
    #################################################################
    ##                          Residuals                          ##
    #################################################################
    
    uu0 = uu1
    e = (y - fhat %*%  t(lam)) 
    uu1 = sum(diag(t(e) %*%  e))
    cat("Iteración:",conteo,"\n")
  }
  
  ###########################################################################
  ###########################################################################
  ###                                                                     ###
  ###                      ORTHO., LOADINGS, FACTORS                      ###
  ###                                                                     ###
  ###########################################################################
  ###########################################################################
  
  # GLOBAL 
  V = (lam[,1]%*% (t(G0) %*%  G0)%*% t(lam[,1]))/t
  evec = eigrs2(V)$evec
  evec=t(lam[,1]) %*% evec[,ncol(evec)]
  G0 = G0 %*% evec
  
  V = (lam[,2]%*% (t(F13) %*%  F13)%*% t(lam[,2]))/t
  evec = eigrs2(V)$evec
  evec=t(lam[,2]) %*% evec[,ncol(evec)]
  F13 = F13 %*% evec
  
  # V = (lam[,3]%*% (t(F12) %*%  F14)%*% t(lam[,3]))/t
  # evec = eigrs2(V)$evec
  # evec=t(lam[,1]) %*% evec[,ncol(evec)]
  # F12 = F12 %*% evec
  # 
  # V = (lam[,4]%*% (t(F23) %*%  F23)%*% t(lam[,4]))/t
  # evec = eigrs2(V)$evec
  # evec=t(lam[,2]) %*% evec[,ncol(evec)]
  # F23 = F23 %*% evec
  
  
  # for(i in 1:1){
  #   evec<-eigrs2(t(fhatreg[,((i-1)*r_reg+1):(i*r_reg)]) %*% fhatreg[,((i-1)*r_reg+1):(i*r_reg)])$evec
  #   fhatreg[,((i-1)*r_reg+1):(i*r_reg)] <-fhatreg[,((i-1)*r_reg+1):(i*r_reg)]%*%evec
  # }
  # 

  #ORTHOGONALIZATION
  
  G0  = G0 -cbind(F13,fhatreg[,1],fhatreg[,2],fhatreg[,3]) %*% t(beta_ols(cbind(F13,fhatreg[,1],fhatreg[,2],fhatreg[,3]),G0))
  F13 = F13-cbind(G0,fhatreg[,1],fhatreg[,2],fhatreg[,3]) %*% t(beta_ols(cbind(G0,fhatreg[,1],fhatreg[,2],fhatreg[,3]),F13))
  F1   = fhatreg[,1] - cbind(G0,F13,fhatreg[,2],fhatreg[,3]) %*% t(beta_ols(cbind(G0,F13,fhatreg[,2],fhatreg[,3]),fhatreg[,1]))
  F2   = fhatreg[,2] - cbind(G0,F13,F1,fhatreg[,3]) %*% t(beta_ols(cbind(G0,F13,F1,fhatreg[,3]),fhatreg[,2]))
  F3   = fhatreg[,3] - cbind(G0,F13,F1,F2) %*% t(beta_ols(cbind(G0,F13,F1,F2),fhatreg[,3]))
  
  # G0  = G0 -cbind(F13,fhatreg[,1],fhatreg[,2],fhatreg[,3]) %*% t(beta_ols(cbind(F13,fhatreg[,1],fhatreg[,2],fhatreg[,3]),G0))
  # F13 = F13-cbind(G0,F24,F34,fhatreg[,1],fhatreg[,2],fhatreg[,3]) %*% t(beta_ols(cbind(G0,F24,F34,fhatreg[,1],fhatreg[,2],fhatreg[,3]),F13))
  # F12 = F24-cbind(G0,F13,F34,fhatreg[,1],fhatreg[,2],fhatreg[,3]) %*% t(beta_ols(cbind(G0,F13,F34,fhatreg[,1],fhatreg[,2],fhatreg[,3]),F24))
  # F23 = F34-cbind(G0,F13,F24,fhatreg[,1],fhatreg[,2],fhatreg[,3]) %*% t(beta_ols(cbind(G0,F13,F24,fhatreg[,1],fhatreg[,2],fhatreg[,3]),F34))
  # F1   = fhatreg[,1] - cbind(G0,F13,F24,F34,fhatreg[,2],fhatreg[,3]) %*% t(beta_ols(cbind(G0,F13,F24,F34,fhatreg[,2],fhatreg[,3]),fhatreg[,1]))
  # F2   = fhatreg[,2] - cbind(G0,F13,F24,F1,fhatreg[,3]) %*% t(beta_ols(cbind(G0,F13,F24,F1,fhatreg[,3]),fhatreg[,2]))
  # F3   = fhatreg[,3] - cbind(G0,F13,F24,F34,F1,F2) %*% t(beta_ols(cbind(G0,F13,F24,F34,F1,F2),fhatreg[,3]))
  # 
  
  ur13=cbind(y[,1:NR1],y[,(NR1+NR2+1):(NR1+NR2+NR3)])
  ur12=y[,1:(NR1+NR2)] 
  ur23=y[,(NR1+1):(NR1+NR2+NR3)]
  
  ur1=R1 
  ur2=R2 
  ur3=R3 
  
  lamG = beta_ols(G0,y) 
  lamr13 = beta_ols(F13,ur13)
  # lamr12 = beta_ols(F12,ur12) 
  # lamr23 = beta_ols(F23,ur23) 
  lamr1 = beta_ols(F1,ur1)
  lamr2 = beta_ols(F2,ur2)
  lamr3 = beta_ols(F3,ur3)
  
  
  lam = cbind(lamG,
              c(lamr13[1:NR1],rep(0,NR2),lamr13[(NR1+1):(NR1+NR3)]),
              #c(lam12[1:(NR1+NR2)],rep(0,NR3)),
              #c(rep(0,NR1),lam23[1:(NR2+NR3)]),
              c(lamr1,rep(0,NR2+NR3)),
              c(rep(0,NR1),lamr2,rep(0,NR3)),
              c(rep(0,NR1+NR2),lamr3))
  
  
  ###########################################################################
  ###########################################################################
  ###                                                                     ###
  ###                     Save                                            ###
  ###                                                                     ###
  ###########################################################################
  ###########################################################################
  
  R[["Factors"]] <- cbind(G0,F13,F1,F2,F3)
  R[["Loadings"]] <- lam
  
  colnames(R[["Factors"]])<-c("G0","F13","F1","F2","F3")
  
  colnames(R[["Loadings"]])<-c("G0","F13","F1","F2","F3")
  
  return(R)
}

# No differences with Matlab code
# factors<-openxlsx::read.xlsx("Factors3B.xlsx",colNames = F)
# round(cor(R$Factors,factors),5)
# lambdas<-openxlsx::read.xlsx("Lambda3B.xlsx",colNames = F)
# round(cor(R$Loadings,lambdas),5)
# lambdas[,1:4]<-lambdas[,1:4]*-1
# factors[,1:4]<-factors[,1:4]*-1
# round(tail(R$Loadings-lambdas,5),10)
# plot(c(R$Factors[,5]),t="l")
# lines(factors[,5],col="red",type = "c")
# round(cor(R$Factors,factors),1)

