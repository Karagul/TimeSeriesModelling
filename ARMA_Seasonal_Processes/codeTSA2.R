# ___________ Script Setup _____________ #

rm(list=ls())
setwd("~/Desktop/TSA/assignment2")
library(expm)
library(timeSeries)
library(tseries)
library(forecast)

# _________ Question 2.2: Univeristy Activity _______________ #

# Variable Setup
t    <- (51:60)
mu <- 1000
Y1    <- c(989,724,1013,1280,1092,811,952,1277,1111,848) 
Y1 <- Y1 - mu
psi1 <- 0.8; psi2 <- 0.7; psi3 <- 0.9

# Predictions Y_t+1, Y_t+2
pred1 <- psi1*Y1[10] - psi2*Y1[9] + psi3*Y1[7] - psi1*psi3*Y1[6] + psi2*psi3*Y1[5] 
pred2 <- psi1*pred1  - psi2*Y1[10] + psi3*Y1[8] - psi1*psi3*Y1[7] + psi2*psi3*Y1[6]
Y1 <- Y1 + mu
pred1 <- pred1+mu
pred2 <- pred2+mu
# 95% Confidence Intervals
upper1 <- pred1 + qt(0.975, 56)*50
lower1 <- pred1 + qt(0.025, 56)*50
upper2 <- pred2 + qt(0.975, 56)*50*sqrt(1+psi1^2)
lower2 <- pred2 + qt(0.025, 56)*50*sqrt(1+psi1^2)
#Plot
plot(t, Y1, type = "l", xlab="year", ylab="number of passed course modules", ylim=c(600,1400),xlim=c(50,63))
lines(61:62, c(pred1,pred2), type="o", lty=2, col="blue" )
lines(60:61, c(Y1[10],pred1), lty=2, col="blue")
points(61:62, c(upper1,upper2), pch="+", col="red")
points(61:62, c(lower1,lower2), pch="+", col="red")




# ______Question 2.3: Random walk ______________ #
# 2.3.3: Simulate a white noise process of 1000 values with mean zero and s^2 = sqrt(2). 
# Then figure out a fast and compact way to calculate the process Yt based on this.

Wnoise <- rnorm(n = 1000, mean = 0, sd = 2^0.25)
Yt <- rep(NA, 1000)
#Initialize 1st value
Yt[1] <- 0.25 + Wnoise[1]
#Then, we let the random walk do the walk
for (i in 2:1000) {
  Yt[i] <- Yt[i-1] + Wnoise[i] 
}
plot(Wnoise, type="l", xlab="time")

# 2.3.4: Simulate 10 realizations of the process Yt: nonstationary as the acf shows a slow decay
W <- matrix(NA, nrow = 1000, ncol = 10)
Y10 <- matrix(NA, nrow = 1000, ncol = 10)
#create 10 white noises
for (i in 1:10) { W[ ,i] <- rnorm(n=1000, mean=0, sd=2^0.25) }
#initialize 1st value
for (i in 1:10) { Y10[1,i] <- 0.25 + W[1,i] }
#let the random walks, walk
for (j in 1:10) {
  for (i in 2:1000) { Y10[i,j] <- Y10[i-1,j] + W[i,j] }
}
#plot of all 10 random walks
plot(Y10[,1], type="l", ylim=c(-100,80), ylab="Yt", xlab="time")
lines(Y10[,2], col="blue")
lines(Y10[,3], col="red")
lines(Y10[,4], col="yellow")
lines(Y10[,5], col="pink" )
lines(Y10[,6], col="brown")
lines(Y10[,7], col="green")
lines(Y10[,8], col="magenta")
lines(Y10[,9], col="purple")
lines(Y10[,10], col="orange")

#from list to numeric output
ACF <- vector("list", 10)
for (i in 1:10) { ACF[i] <- acf(Y10[,i], plot = FALSE) }
ACFp <- sapply(1:length(ACF), function(i) as.numeric(ACF[[i]]))
#plot their estimated autocorrelation functions 
plot(ACFp[,1], ylim=c(0,1), type="l", ylab="ACF", xlab="lag")
lines(ACFp[,2], col="blue")
lines(ACFp[,3], col="red")
lines(ACFp[,4], col="yellow")
lines(ACFp[,5], col="pink" )
lines(ACFp[,6], col="brown")
lines(ACFp[,7], col="green")
lines(ACFp[,8], col="magenta")
lines(ACFp[,9], col="purple")
lines(ACFp[,10], col="orange")

# ______ Question 2.4: Simulating seasonal processes ______________ #

# Simulate the following models (where monthly data are assumed)
# Notice we show 60 lags to see what happens in 5 12month seasons
#1: AR(1) process w/out seasonality  (or2)
arma1 <- arima.sim(n = 300, list(ar = c(rep(0,11),-0.6), ma = c(0, 0)) )
plot(arma1)
acf(arma1, lag.max = 60)
pacf(arma1, lag.max = 60)
#2                                    (or1)
arma2 <-   arima.sim(n = 300, list(ar = c(-0.8,0), ma = c(0, 0)) )
plot(arma2)
acf(arma2, lag.max = 60)
pacf(arma2, lag.max = 60)
#3                                   (or3)
arma3 <- arima.sim(n = 300, list(ar = c(-0.8,0), ma = c(rep(0,11), 0.5)) )
plot(arma3)
acf(arma3, lag.max = 60)
pacf(arma3, lag.max = 60)
#4                                     (or5)
arma4 <- arima.sim(n = 300, list(order = c(0,0,12), ma = c(0.5,rep(0,10),0.4)))
plot(arma4)
acf(arma4, lag.max = 60)
pacf(arma4, lag.max = 60)
#5                                     (or4)
arma5 <- arima.sim(n = 300, list(order = c(0,0,12), ar = c(-0.8,rep(0,10),-0.6))) #bug
arma52 <- arima.sim(n = 300, list(ar = c(rep(0,11),-0.9*-0.7) ) ) #bug?
plot(arma52)
acf(arma5, lag.max = 60)
pacf(arma5, lag.max = 60)
#6                                     (or6)
arma6 <-  arima.sim(n = 300, list(ar = c(rep(0,11),-0.6), ma=c(0.5)))
plot(arma6)
acf(arma6, lag.max = 60)
pacf(arma6, lag.max = 60)


