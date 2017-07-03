# TSA Assignment 4 #

############################### Script Setup & Data ##############################################################

rm(list=ls())
setwd("~/Desktop/TSA/Assignment4")
library(expm)
library(timeSeries)
library(tseries)
library(forecast)
library(lmtest)
library(PBSmodelling)
library(marima)

#Data Partitioning
data <- read.table("A4_gps_log.csv",header = T,sep = ",")
data.train <- data[1:950,]
data.test  <- data[951:1000,]

###############################  PLOTTING  ########################################################################

# Plot entire original series
dev.off()
par(mfrow=c(2,1))
#mlat
plot(data.train[,1], data.train[,2], ylab='mlat', xlab ="Time (s)", type = "l", xlim=c(0,1000), ylim=c(-6,7), lwd=3)
lines(data.test[,1], data.test[,2], col = "green", lwd=3)
legend(0,0, c("Training", "Test (for predictions)"), col = c("black", "green"), lty = c(1,1), bty ="n")

#mlong
plot(data.train[,1], data.train[,3], ylab='mlong', xlab ="Time (s)", type = "l", xlim=c(0,1000), ylim=c(-15,4), lwd=3)
lines(data.test[,1], data.test[,3], col = "green", lwd=3)
legend(0,0, c("Training", "Test (for predictions)"), col = c("black", "green"), lty = c(1,1), bty ="n")

#mlat and mlong 
plot(data.train[,3], data.train[,2], ylab='mlat', xlab ="mlong", type = "l", ylim=c(-6,7), xlim=c(-15,4), lwd=3)
lines(data.test[,3], data.test[,2], col = "green", lwd=3)
legend(-15,4, c("Training", "Test (for predictions)"), col = c("black", "green"), lty = c(1,1), bty ="n")
dev.off()

#CrossCorrelation Plots
acf(data.train[,-1])
pacf(data.train[,-1])

#Difference Data      (1st difference)
diff.mlat <- diff(data.train[,2])
diff.mlong <- diff(data.train[,3])
diff.data <- ts(cbind(diff.mlat, diff.mlong))

#CrossCorrelation Plots of Differenced Data
acf(diff.data)
pacf(diff.data)

#Difference Data AGAIN (2nd difference)
diff2.mlat <- diff(diff(data.train[,2]))
diff2.mlong <- diff(diff(data.train[,3]))
diff2.data <- ts(cbind(diff2.mlat, diff2.mlong))

#CrossCorrelation Plots of Differenced Data
acf(diff2.data)
pacf(diff2.data)

############################### Fitting MARIMA Models ##############################################################

# Differencing data for marima analysis
diff_data  <- define.dif(data.train[,2:3], difference=c(1,1,2,1))
diff2_data <- define.dif(data.train[,2:3], difference=c(1,2,2,2))

#Fitting MARIMA Model 0 (0,2,3)
ar<-c(1)
ma<-c(2)
# Define the proper model:
Model0 <- define.model(kvar=2, ar=ar, ma=ma, indep = NULL) 
# Difference the data
dif.y <-define.dif(t(data.train[,2:3]), c(1,2,2,2))
# Now call marima:
Marima0 <- marima(dif.y$y.dif,means=1,
                  ar.pattern=Model0$ar.pattern, ma.pattern=Model0$ma.pattern, Check=FALSE, Plot="log.det", penalty=0.0)
# print AR and MA estimates
Marima0
# MA f-values are quite high; what does that mean?
#MA f-values (squared t-values): Lag=3
#     x1=y1  x2=y2
# y1  33.56 41.58
# y2  30.15 16.08

## Looking at acf and ccf of residuals 
acf(t(Marima0$residuals))
pacf(t(Marima0$residuals))


#Fitting MARIMA Model 1 (0,1,1)
ar<-c(0)
ma<-c(1)
# Define the proper model:
Model1 <- define.model(kvar=2, ar=ar, ma=ma, indep = NULL) 
# Difference the data
dif.y <-define.dif(t(data.train[,2:3]), c(1,1,2,1))
# Now call marima:
Marima1 <- marima(dif.y$y.dif,means=1,
                  ar.pattern=Model1$ar.pattern, ma.pattern=Model1$ma.pattern, Check=FALSE, Plot="log.det", penalty=0.0)
# print AR and MA estimates
Marima1
# MA f-values are quite high; what does that mean?
#MA f-values (squared t-values): Lag=1
#     x1=y1  x2=y2
# y1 301.45  10.34
# y2  16.24 237.77

## Looking at acf and ccf of residuals 
acf(t(Marima1$residuals))
pacf(t(Marima1$residuals))


#Fitting MARIMA Model 2
ar<-c(1)
ma<-c(1)
# Define the proper model:
Model2 <- define.model(kvar=2, ar=ar, ma=ma) 
# Difference the data
dif.y <-define.dif(t(data.train[,2:3]), c(1,1,2,1))
# Now call marima:
Marima2 <- marima(dif.y$y.dif,means=1,
                  ar.pattern=Model2$ar.pattern, ma.pattern=Model2$ma.pattern, Check=FALSE, Plot="log.det", penalty=0.0)
# print AR and MA estimates
Marima2
#AR f-values (squared t-values):Lag=1
#     x1=y1  x2=y2
# y1 4.2230 1.5879
# y2 0.4252 5.8607
#
#MA f-values (squared t-values): Lag=1
#     x1=y1 x2=y2
# y1 48.44  6.38
# y2  6.03 23.65
## Looking at acf and ccf of residuals 
acf(t(Marima2$residuals[,-1]))
pacf(t(Marima2$residuals[,-1]))


#Fitting MARIMA Model 3
ar<-c(1)
ma<-c(2)
# Define the proper model:
Model3 <- define.model(kvar=2, ar=ar, ma=ma) 
# Difference the data
difference = matrix(c(1,1,2,1), nrow=2)   #difference = c(1,1,2,1)
dif.y <-define.dif(t(data.train[,2:3]), difference=difference)
# Now call marima:
Marima3 <- marima(dif.y$y.dif,means=1,
                  ar.pattern=Model3$ar.pattern, ma.pattern=Model3$ma.pattern, Check=FALSE, Plot="log.det", penalty=0.0)
# print AR and MA estimates
Marima3

# AR f-values (squared t-values): Lag=1
#      x1=y1    x2=y2
# y1 310.3038  10.3898
# y2  15.2963 241.5915

# MA f-values (squared t-values):Lag=2
#   x1=y1 x2=y2
# y1 52.03 16.90
# y2 14.46 25.04
## Looking at acf and ccf of residuals 
acf(t(Marima3$residuals[,-1]))
pacf(t(Marima3$residuals[,-1]))


#Reducing MARIMA Model 3
Marima3.reduced <- marima(dif.y$y.dif,means=1,
                  ar.pattern=Model3$ar.pattern, ma.pattern=Model3$ma.pattern, Check=FALSE, Plot="log.det",  penalty=2.0)
# print AR and MA estimates
Marima3.reduced

Marima3.reduced <- marima(dif.y$y.dif,means=1,
                          ar.pattern=Model3$ar.pattern, ma.pattern=Model3$ma.pattern, Check=FALSE, Plot="log.det",  penalty=3.84)
# AR f-values (squared t-values): Lag=1
#    x1=y1    x2=y2
# y1 310.3038  10.3898
# y2  15.2963 241.5915

# MA f-values (squared t-values): Lag=2
#    x1=y1 x2=y2
# y1 52.03 16.90
# y2 14.46 25.04
## Looking at acf and ccf of residuals 
acf(t(Marima3.reduced$residuals[,-1]))
pacf(t(Marima3.reduced$residuals[,-1]))


#Fitting MARIMA Model 4
ar<-c(1)
ma<-c(3)
# Define the proper model:
Model4 <- define.model(kvar=2, ar=ar, ma=ma) 
# Difference the data
difference = c(1,1,2,1)
dif.y <-define.dif(t(data.train[,2:3]), difference=difference)
# Now call marima:
Marima4 <- marima(dif.y$y.dif,means=1,
                  ar.pattern=Model4$ar.pattern, ma.pattern=Model4$ma.pattern, Check=FALSE, Plot="log.det", penalty=0.0)
# print AR and MA estimates
Marima4

# AR f-values (squared t-values): Lag=1
#    x1=y1    x2=y2
# y1 252.2008   0.7261
# y2   3.8175 213.6531

# MA f-values (squared t-values): Lag=3
#    x1=y1 x2=y2
# y1  9.89  2.64
# y2  0.64  4.42
acf(t(Marima4$residuals[,-1]))
pacf(t(Marima4$residuals[,-1]))



############################### MODEL VALIDATION ##############################################################

ks.test(Marima3.reduced$residuals[1,], y="pnorm", alternative = c("two.sided"))
ks.test(Marima3.reduced$residuals[2,], y="pnorm", alternative = c("two.sided"))


# QQ Plot
#Both variates
qqnorm(Marima3.reduced$residuals)
qqline(Marima3.reduced$residuals,col=2)

#Mlat
qqnorm(Marima3$residuals[1,])
qqline(Marima3$residuals,col=2)
#Mlong
qqnorm(Marima3$residuals[2,])
qqline(Marima3$residuals,col=2)

# Cumulative Periodogram  PLot
#     - do residuals look like white noise? 

cpgram(Marima3.reduced$residuals[1,]) # not quite!
cpgram(Marima3.reduced$residuals[2,]) # not quite!



############################### PREDICTIONS ##############################################################

## Predicting 50 years ahead
pred.data <- as.matrix(cbind(t(data.train[,2:3]), matrix(NA,nrow=2, ncol=50)))
pred <-  arma.forecast(pred.data, nstart=ncol(t(data.train[,2:3])), nstep=50, marima=Marima3.reduced, check=TRUE, dif.poly = dif.y$dif.poly)


#first, difference validation set
mlat.diff.test <- diff(data.test[,2])
mlong.diff.test <- diff(data.test[,3])
data.test.diff <- ts(cbind(mlat.diff.test, mlong.diff.test))



# Plot entire original series + forecasts
dev.off()
par(mfrow=c(3,1))
#mlat
plot(data.train[,1], data.train[,2], ylab='mlat', xlab ="Time (s)", type = "l", xlim=c(800,1000), ylim=c(-6,-2), lwd=3)
lines(data.test[,1], data.test[,2], col = "green", lwd=3)
lines(data.test[,1], pred$forecasts[1,951:1000], col = "red", lwd=3)
pred.int <- pred$forecasts[1,951:1000] + cbind(rep(0, 50), -1, 1)*qnorm(0.975)*sqrt(pred$pred.var[1,1,])
matlines(951:1000, pred.int, lty=c(1,2,2), col=2, lwd=2 )
legend(800,-4, c("Training", "Observations","Predictions"), col = c("black", "green", "red"), lty = c(1,1,1), bty ="n")

#mlong
plot(data.train[,1], data.train[,3], ylab='mlong', xlab ="Time (s)", type = "l", xlim=c(800,1000), ylim=c(0,4), lwd=3)
lines(data.test[,1], data.test[,3], col = "green", lwd=3)
lines(data.test[,1], pred$forecasts[2,951:1000], col = "red", lwd=3)
pred.int <- pred$forecasts[2,951:1000] + cbind(rep(0, 50), -1, 1)*qnorm(0.975)*sqrt(pred$pred.var[1,1,])
matlines(951:1000, pred.int, lty=c(1,2,2), col=2, lwd=2 )
legend(800,4, c("Training", "Observations", "Predictions"), col = c("black", "green", "red"), lty = c(1,1,1), bty ="n")



############################### KALMAN FILTER ############################################################################

kalman_me <- function(Y,A,C,S1,S2,mu0,v0,steps){
  
  ####################### FUNCTION DESCRIPTION ########################################################################
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




###################  Kalman Implementation ##########################################################
# Kalman Input
Y   = data.train[,2:3]
A   = matrix( c(1,0,0,0,1,1,0,0,0,0,1,0,0,0,1,1), ncol=4)
C   = matrix( c(1,0,0,0,0,1,0,0), ncol=4)
S1  = diag(c(0.0003, 1e-5, 0.0003, 1e-5))
S2  = diag(c(1e-4,1e-4))
mu0 = c(data.train[1,2], 0, data.train[1,3], 0)
v0  = S1
steps = 50
  
# Kalman Implementation
k.output1 <- kalman_me( Y=Y, A=A, C=C, S1=S1, S2=S2, mu0=mu0, v0=v0, steps=steps)

# Plotting Kalman Predictions
plot(data.train[,1], data.train[,2], ylab='mlat', xlab ="Time (s)", type = "l", xlim=c(800,1000), ylim=c(-8,-2), lwd=3)
lines(data.test[,1], data.test[,2], col = "green", lwd=3)
lines(data.test[,1], k.output1$X.pred[951:1000,1], col="red")
cl  <- k.output1$X.pred[951:1000,1] + qnorm(0.975)*sqrt(k.output1$Sxx.pred[1,1,951:1000])
cl2 <-  k.output1$X.pred[951:1000,1] - qnorm(0.975)*sqrt(k.output1$Sxx.pred[1,1,951:1000])
lines(data.test[,1],cl)
lines(data.test[,1], cl2)
legend(850,-5, c("Training", "Observations", "Kalman Predictions"), col = c("black", "green","red"), lty = c(1,1,1), bty ="n")

plot(data.train[,1], data.train[,3], ylab='mlong', xlab ="Time (s)", type = "l", xlim=c(800,1000), ylim=c(0,5), lwd=3)
lines(data.test[,1], data.test[,3], col = "green", lwd=3)
lines(data.test[,1], k.output1$X.pred[951:1000,3], col="red")
cl  <- k.output1$X.pred[951:1000,3] + qnorm(0.975)*sqrt(k.output1$Sxx.pred[1,1,951:1000])
cl2 <-  k.output1$X.pred[951:1000,3] - qnorm(0.975)*sqrt(k.output1$Sxx.pred[1,1,951:1000])
lines(data.test[,1],cl)
lines(data.test[,1], cl2)
legend(800,4, c("Training", "Observations", "Kalman Predictions"), col = c("black", "green","red"), lty = c(1,1,1), bty ="n")

plot(1:950, k.output1$X.rec[1:950,2], type="l", ylab="velocity mlat", ylim=c(0.08, -0.15))
cl  <- k.output1$X.rec[1:950,2] + qnorm(0.975)*sqrt(k.output1$Sxx.pred[2,2,951:1000])
cl2 <-  k.output1$X.rec[1:950,2] - qnorm(0.975)*sqrt(k.output1$Sxx.pred[2,2,951:1000])
lines(1:950,cl, col="red")
lines(1:950, cl2, col="red")
legend(0,-0.15, c("Velocity mlat", "95% CL"), col = c("black", "red"), lty = c(1,1), bty ="n")

plot(1:950, k.output1$X.rec[1:950,4], type="l", ylab="velocity mlong", ylim=c(0.12, -0.15))
cl  <- k.output1$X.rec[1:950,4] + qnorm(0.975)*sqrt(k.output1$Sxx.pred[4,4,951:1000])
cl2 <-  k.output1$X.rec[1:950,4] - qnorm(0.975)*sqrt(k.output1$Sxx.pred[4,4,951:1000])
lines(1:950,cl, col="red")
lines(1:950, cl2, col="red")
legend(0,-0.15, c("Velocity mlong", "95% CL"), col = c("black", "red"), lty = c(1,1), bty ="n")
