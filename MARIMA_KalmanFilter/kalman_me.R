kalman_me <- function(Y,A,C,S1,S2,mu0,v0,steps){
  
  ####################### FUNCTION DESCRIPTION ########################################################################
  #
  # DESCRIPTION: This is an implementation of a Kalman Filter used for multivariate time series modelling
  #
  # INPUT
  #   Y : nxp matrix witn n samples (time periods) by p observations (i.e. mlat and mlong) per sampling interval
  #   A : sxs matrix with s being number of components in sate vector [xt vxt yt vyt]
  #   C : pxs matrix with p observations (i.e. mlat and mlong) by s components from state vector [xt vxt yt vyt]
  #   S1: sxs diagonal variance matrix of system      model error
  #   S2: pxp diagonal variance matrix of observation model error
  #   mu0: initial estimate of state-space model at time=1 given X0. 
  #        an sx1 vector. s is number of states in state vector [xt vxt yt vyt]      
  #   v0 : initial covariance of Sxx at time=1 given Sxx0 (a 2x2 matrix)
  #   steps : number of predictions after last observation in Y (a scalar)
  #
  #
  # OUTPUT
  #   REC:  list of reconstruced variables: X, Sxx, Syy, K
  #   PRED: list of predicted    variables: X, Sxx, Syy
  #
  ####################################################################################################################
  
  
  
  # INITIALIZATION
  Xhat <- mu0
  Sxx  <- v0
  Syy  <- C %*% Sxx %*% t(C) + S2
  
  N = dim(Y)[1]       #number of observations
  nstates = dim(A)[1] #number of states
  
  # STORAGE
  #(For storing variables per iteration)
  X.rec    <- array(dim=c(N+steps,nstates)) #A single X is a row vector
  X.pred   <- array(dim=c(N+steps,nstates))
  Sxx.rec  <- array(dim=c(dim(Sxx), N+steps))
  Sxx.pred <- array(dim=c(dim(Sxx), N+steps))
  Syy.rec  <- array(dim=c(dim(Sxx), N+steps))
  Syy.pred <- array(dim=c(dim(Sxx), N+steps))
  K.all    <- array(dim=c(dim(Sxx %*% t(C) %*% solve(Syy)), N+steps))

  
  
  for (i in 1:(N+steps)) {
    
    if (i<N+1) {
      # RECONSTRUCTION 
      K    <- Sxx %*% t(C) %*% solve(Syy)            #kalman gain
      Xhat <- Xhat + K%*%t(Y[i,]) - K%*%C%*%Xhat     #Xhat update, using kalman gain (reduces noise)
      Sxx  <- Sxx - K%*%Syy%*%t(K)                   #Sxx  update
      
      #storing
      K.all[,,i]    <- K 
      X.rec[i,]     <- Xhat
      Sxx.rec[,,i]  <- Sxx
      Syy.rec[,,i]  <- Syy
    }#endif
  
    # PREDICTION
    Xhat <- A%*%Xhat
    Sxx  <- A%*%Sxx%*%t(A) + S1
    Syy  <- C%*%Sxx%*%t(C) + S2
    
    #storing
    X.pred[i,]      <- Xhat
    Sxx.pred[,,i]   <- Sxx
    Syy.pred[,,i]   <- Syy
    
  }#endfor 
  
  # OUTPUT
  output  <- list(X.rec=X.rec,  Sxx.rec=Sxx.rec,  Syy.rec=Syy.rec,  K=K.all, X.pred=X.pred, Sxx.pred=Sxx.pred, Syy.pred=Syy.pred)
  return(output)
  
}#endfunction

