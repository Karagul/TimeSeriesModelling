rm(list=ls())
setwd("~/Desktop/TSA/assignment1")
library(expm)
library(forecast)
data <- read.table("apartment_prices.csv",header = T,sep = "")


##############################################################
## 1.1: Plot apartment price as function of time
Y <- ts(log(data[ 1:92, 3]))
Yfull <- ts(log(data[1:97, 3]))
plot(data$time, Yfull, pch = 16, type = "l", xlab = "Year", ylab = "Apartment Prices log(kr)" )

## 1.2: Does the global mean value and standard deviation of the prices give a reasonable representation of the data?
mean = mean(Y)


## 1.3: Formulate GLM: Y = X*theta + error
X <- cbind(rep(1, 92), data$time[1:92]) 
Xfull <- cbind(rep(1, 97), data$time[1:97]) 
theta_hat <- solve(t(X)%*%X, t(X)%*%Y)
Y_hat = Xfull%*%theta_hat

## check
lreg <- lm( Y ~ 0 + X )
summary(lreg)
plot(lreg)

#Plot GLM with predictions
plot(data$time, Yfull, type = "l", xlab = "Year", ylab = "Apartment Prices log(kr)" )
lines(data$time[1:92], lreg$fitted.values, type = "l", col = "red")
lines(data$time[92:97],  lreg$coefficients[2]*seq(2014.75, 2016, 0.25) + lreg$coefficients[1], type = "l" , col = "blue")

## NOTE: OLS is only appropriate when the variance of the errors is constant and 
##       the residuals are mutually independent



##########################################################
## 1.4 : EXPONENTIAL SMOOTHING MODEL
lambda <- 0.8  
Y_hatExpS <- rep(0, 93)

#We initialize the 1st and 2nd forecasts by using Y[1] as starting point (S0)
Y_hatExpS[1] <- Y[1]      #We use 1st observation as 1st forecast 
Y_hatExpS[2] <- ( 1-lambda ) * ( Y[1] )  +  lambda * Y_hatExpS[1]

#We iterate to find the remaining forecasts using eq. 3.74 from book
for (i in 3:93) {
  Y_hatExpS[i] <-  ( 1-lambda ) * ( Y[i-1] )  +  lambda * Y_hatExpS[i-1]
}

#Predictions: Our last observation is Y[92], so we use eq. 3.75 for predicting 2015 and 2016
for (i in 93:97) {
  Y_hatExpS[i] <-  ( 1-lambda ) * ( Y[92] )  +  lambda * Y_hatExpS[i-1]
}

#Check with Holtwinters function
ExSmooth <- HoltWinters( Yfull[1:92] , alpha = (1 - lambda), beta = F, gamma = F, l.start = Yfull[1]  )
forecastsEX <- forecast::forecast.HoltWinters( ExSmooth, h=5, level = .95 )
#PLOT: Exponential Smoothing Model, lambda = 0.8, with predictions and 95% CL
plot(forecastsEX, main = "Exponential Smoothing and Predictions with 95% CL")
lines(ExSmooth$fitted[,1], lty = 2 , col = "green")
lines(92:97, Yfull[92:97], type = "l")


####################################################################
### 1.5: LOCAL LINEAR TREND MODEL: 
lambda <- 0.8
p <- 2

## USEFUL FUNCTIONS to set up problem
#1
x_n  <- function(N) {
  a <- matrix( ncol=2, nrow=N )
  for (i in 1:N) {
    a[i,] <- fgen(-N +i)
  }
  return(a) #returns Nx2 matrix of ones (col1) and -N going down to 0 (col2)
}

#2
invSIGMA_n <- function(N, lambda) {
  invS <- matrix(ncol = 1, nrow = N)
  for (i in 1:N) {
    invS[i] <- lambda^(N-i)
  }
  return( Diagonal(N , invS) ) #returns NxN diagonal matrix; it's the inverse of SIGMA
} 

#3: returns 1x2 matrix
fgen <- function(N)  matrix( c(1,N), nrow=1 )  

#4 : [as eq 3.100] : To initialize F_n
F_n <- function(N, lambda) {
  F <- matrix(data = 0, nrow = 2, ncol = 2)
  for (i in 0: (N-1) ) {
    F <- F + lambda^i * t(fgen(-i)) %*% fgen(-i)
  }
  return(F) #returns 2x2 matrix
}

#5 : [as eq 3.100] : To initialize h_n
h_n <- function(N, lambda) {
  h = matrix(data = 0, nrow = 2, ncol = 1)
  for (i in 0:(N-1) ) {
    h <- h + lambda^i * t(fgen(-i)) %*% Y[N-i]  
  }
  return(h) #returns 2x1 matrix
}

#6 : Theta Estimate : [as eq 3.99]
theta.hatLT <- function(N, lambda) { solve(F_n(N, lambda)) %*% h_n(N, lambda) }

#7 : Predict Y.. i.e. Y_hat [as eq 3.101]
Y.hatLT <- function(N, l, lambda) { fgen(l) %*% theta.hatLT(N, lambda) }

#8 : as in slides 
total_memory_n <- function(N, lambda) {
  return((1 - lambda^(N)) / (1 - lambda))
} 


 

###  PREDICTIONS: YHAT  
Y_hatLT <- rep(NA, 97)
for (i in 2:91) { #We start at 2 bc Fn[1] is singular --> no inverse to calculate theta 
  Y_hatLT[i+1] <- Y.hatLT(i,1, lambda)
}
for (i in 92:96) { #Here we start making predictions L steps ahead from last observation: Y[92]
  Y_hatLT[i+1] <- Y.hatLT(92,i-91, lambda)
}

### ESTIMATES: Using F_n(0) instead of F_n(1)
EstimatesLT <- rep(NA, 92)
for (i in 2:91) { 
  EstimatesLT[i+1] <- Y.hatLT(i,0, lambda)
}

### UNCERTAINTY: Global and Local Sigma estimator, Variance of Prediction error, and 95% confidence interval:

## Global estimator of sigma2 [as in slide 26 in L#4]
residuals <- as.vector(Y[3:92] - Y_hatLT[3:92])
global.est  <- rep(NA, 90)
for (i in 1:90) {
  global.est[i] <- residuals[i]^2 / ( 1 + fgen(1) %*% solve(F_n(i+2, lambda)) %*% t(fgen(1)) ) 
}
sigma2.global <- (sum(global.est) / (90-2))
sigma.global <- sqrt(sigma2.global) #1626

## Local Estimator of sigma2 [as slide 26 in L#4]
sigma2.local<- ( t(residuals)%*%invSIGMA_n(90, lambda)%*%residuals) /(5-2) #/(T-p)
sigma.local <- sqrt(sigma2.local)

## Variance of Prediction Error : [as eq 3.102]
var_PredictionErrorLT.local <- rep(NA, 5)
var_PredictionErrorLT.global <- rep(NA, 5)
for (i in 92:96) { 
  var_PredictionErrorLT.local[i-91] <- sigma2.local  * (1 + fgen(i-91) %*% solve(F_n(i, lambda)) %*% t(fgen(i-91)))  
  var_PredictionErrorLT.global[i-91] <- sigma2.global  * (1 + fgen(i-91) %*% solve(F_n(i, lambda)) %*% t(fgen(i-91)))  
}

## 95% CL
UpperCL <- matrix(NA,nrow=5, ncol=2) #1st col for local, 2nd for global
LowerCL <- matrix(NA,nrow=5, ncol=2)
for (i in 93:97) {
  #With Local Sigma2 Estimator
  UpperCL[i-92,1] <- Y_hatLT[i]  +  qt( 0.975, df = 5 - p) * sqrt(var_PredictionErrorLT.local[i-92]) 
  LowerCL[i-92,1] <- Y_hatLT[i]  +  qt( 0.025, df = 5 - p) * sqrt(var_PredictionErrorLT.local[i-92]) 
  
  #With Global Sigma2 Estimator
  UpperCL[i-92,2] <- Y_hatLT[i]  +  qt( 0.975, df = 92 - p) * sqrt(var_PredictionErrorLT.global[i-92]) 
  LowerCL[i-92,2] <- Y_hatLT[i]  +  qt( 0.025, df = 92 - p) * sqrt(var_PredictionErrorLT.global[i-92]) 
}
  

## PLOT: Local Trend Model, lambda = 0.8, with estimates (f_n(0)), and predictions with 95% CL
plot(data$time[1:97], Yfull, type = "l", main = "Local Linear Trend Model, Estimates, and Predictions with 95% CL", ylim = c(8.7,11))
points(data$time[1:97], Y_hatLT[1:97], lty = 1 , col = "green", pch = "+")
lines(data$time[1:97], EstimatesLT[1:97], type = "l", col = "red")
lines(data$time[93:97], UpperCL[1:5,1], type = "l", col = "blue")
lines(data$time[93:97], LowerCL[1:5,1], type = "l", col = "blue")
lines(data$time[93:97], UpperCL[1:5,2], type = "l", col = "yellow")
lines(data$time[93:97], LowerCL[1:5,2], type = "l", col = "yellow")
legend( 1992,11,  c("1-Step Predictions", "Estimates", "95% CL:local estimator", "95% CL:global estimator"),
        lty = c(1,1,1,1), col =c("green",'red',"blue", "yellow") , bty = "n")



#########################################################################
## 1.6: FIND OPTIMAL LAMBDA VALUE: Formulate SSE as a function of lambda
#(Hint: disregard the first 20 1-step predictions as burn in period)
LAMBDA <- seq(0.01,1, by=0.01)
SSE <- matrix(0, nrow = length(LAMBDA), ncol = 1)
Y_hat.temp <- rep(0, 92) #We dont use 2015/16 for lambda training
temp <- rep(0, 92)

for (j in 1:length(LAMBDA)) {
  
  #Predictions under different lambda values
  #Change range to 2:51 if you want to find optimal lambda with data prior to 2005, else 20:91
  for (i in 2:51) {         
    Y_hat.temp[i+1] <- Y.hatLT(i,1, LAMBDA[j])
    temp[i+1]   <- (Y[i+1] - Y_hat.temp[i+1])^2  #squared residuals
  }
  SSE[j,1] <- sum(temp) 
}
plot(LAMBDA, SSE)

## Optimal Lambda with 20 obs burn-in period
index <- which( min(SSE) == SSE )
lambda_opt <- LAMBDA[index] #lambda_opt = 0.53
## Optimal Lambda using data prior to 2005
lambda_opt2005 <- LAMBDA[index] #0.69





#################### EXTRA: Not Part of Assignment ###############################################
#Predictions with new found optimal lambdas
Y_hatLT.opt <- rep(NA, 97)
Y_hatLT.opt2005 <- rep(NA, 97)

for (i in 20:92) {  
  Y_hatLT.opt[i+1] <- Y.hatLT(i,1, lambda_opt)
  Y_hatLT.opt2005[i+1] <- Y.hatLT(i,1, lambda_opt2005)
  
}
for (i in 93:97) { #Here we start making predictions L steps ahead from last observation: Y[92]
  Y_hatLT.opt[i] <- Y.hatLT(92,i-92, lambda_opt)
  Y_hatLT.opt2005[i] <- Y.hatLT(92,i-92, lambda_opt2005)
  
}

#Uncertainty: Variance of Prediction error, and 95% confidence interval:
var_PredictionErrorLT.opt <- rep(NA, 97)
var_PredictionErrorLT.opt2005 <- rep(NA, 97)

for (i in 20:92) { 
  var_PredictionErrorLT.opt[i+1] <- varErrorLT(i,1, lambda_opt)
  var_PredictionErrorLT.opt2005[i+1] <- varErrorLT(i,1, lambda_opt2005)
  
}
for (i in 93:96) { 
  var_PredictionErrorLT.opt[i+1] <- varErrorLT(i,1, lambda_opt)
  var_PredictionErrorLT.opt2005[i+1] <- varErrorLT(i,1, lambda_opt2005)
  
}

UpperCL.opt <- rep(NA, 97)
LowerCL.opt <- rep(NA, 97)
UpperCL.opt2005 <- rep(NA, 97)
LowerCL.opt2005 <- rep(NA, 97)
for (i in 20:93) { #We start at 3 bc we have no Yhat before that
  UpperCL.opt[i] <- Y_hatLT.opt[i]  +  qt( 0.975, df = i - p) * sqrt(var_PredictionErrorLT.opt[i])
  LowerCL.opt[i] <- Y_hatLT.opt[i]  +  qt( 0.025, df = i - p) * sqrt(var_PredictionErrorLT.opt[i])
  UpperCL.opt2005[i] <- Y_hatLT.opt2005[i]  +  qt( 0.975, df = i - p) * sqrt(var_PredictionErrorLT.opt2005[i])
  LowerCL.opt2005[i] <- Y_hatLT.opt2005[i]  +  qt( 0.025, df = i - p) * sqrt(var_PredictionErrorLT.opt2005[i])
}
for (i in 94:97) { #We calculate the last 4 using the last availalbe variance
  UpperCL.opt[i] <- Y_hatLT.opt[i]  +  qt( 0.975, df = i - p) * sqrt(var_PredictionErrorLT.opt[93])
  LowerCL.opt[i] <- Y_hatLT.opt[i]  +  qt( 0.025, df = i - p) * sqrt(var_PredictionErrorLT.opt[93])
  UpperCL.opt2005[i] <- Y_hatLT.opt2005[i]  +  qt( 0.975, df = i - p) * sqrt(var_PredictionErrorLT.opt2005[93])
  LowerCL.opt2005[i] <- Y_hatLT.opt2005[i]  +  qt( 0.025, df = i - p) * sqrt(var_PredictionErrorLT.opt2005[93])
}

########################
#PLOT
png(filename = "~/Desktop/TSA/assignment1/Ltrend_Lopt_Lopt2005.png")
plot( data$time[1:97], Yfull, xlim=c(1992, 2016), ylim=c(8.7, 10.9))
#lines( data$time[1:97], Y_hatExpS[1:97], col="blue", xlim=c(1992, 2016))
lines( data$time[1:97], UpperCL[1:97], col="red", xlim=c(1992, 2016))
lines( data$time[1:97], LowerCL[1:97], col="red", xlim=c(1992, 2016))
lines( data$time[1:97], Y_hatLT[1:97], col="red")
lines( data$time[1:97], UpperCL.opt[1:97], col="green", xlim=c(1992, 2016))
lines( data$time[1:97], LowerCL.opt[1:97], col="green", xlim=c(1992, 2016))
lines( data$time[1:97], Y_hatLT.opt[1:97], col="green")
lines( data$time[1:97], UpperCL.opt2005[1:97], col="blue", xlim=c(1992, 2016))
lines( data$time[1:97], LowerCL.opt2005[1:97], col="blue", xlim=c(1992, 2016))
lines( data$time[1:97], Y_hatLT.opt2005[1:97], col="blue")
dev.off()













