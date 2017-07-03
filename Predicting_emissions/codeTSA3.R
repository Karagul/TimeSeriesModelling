# TSA Assignment 3 #

#__________________ Script Setup & Data _____________________________#

rm(list=ls())
setwd("~/Desktop/TSA/Assignment3")
library(expm)
library(timeSeries)
library(tseries)
library(forecast)
library(lmtest)
library(PBSmodelling)

data <- read.table("hcoe17.csv",header = T,sep = ";")
xts  <- ts(rev(data[1:735,3]))  #reversing data!  #, frequency = 24
lnxts <- log(xts) #log transformation of data to subdue increase in variance

# Subsetting Training Data (ie excluding last 48 hr observations for testing)
lnxts.train <- lnxts[1:687] 


# __________________ Visualizing the Data ____________________________#

# Plot entire original series
plot(xts[1:687], ylab='NOx', xlab ="Time (Hrs): Starting Jan 7 2017 00:00", type = "l", xlim=c(0,750))
lines(688:735, xts[688:735], col = "green")
legend(0,110, c("Training", "Test (for predictions)"), col = c("black", "green"), lty = c(1,1), bty ="n")

# Plot entire log series
plot(lnxts[1:687], ylab='ln(NOx)', xlab ="Time (Hrs): Starting Jan 7 2017 00:00", type = "l", xlim=c(0,750), ylim=c(0,5))
lines(688:735, lnxts[688:735], col = "green")
legend(0,5, c("Training", "Test (for predictions)"), col = c("black", "green"), lty = c(1,1), bty ="n")


# Plot weekly series, starting on Saturday morning 12:01am (aka Friday after midnight)
plot(lnxts[1:168], type ="l", col = "black", #1st week
xaxt = "n", ylab = "ln(NOx)", xlab="", ylim = c(0,7)) #
axis(side=1, at=seq(1,169,24), labels=c("Saturday", "Sunday", "Monday", "Tuesday",
                                        "Wednesday", "Thursday", "Friday", "Saturday"),
     pos=0, lty=c(1,1,1,1,1,1,1,1), las=2)
lines(lnxts[169:336], col = "blue")          #2nd week
lines(lnxts[337:504], col = "green")         #3rd week
lines(lnxts[505:672], col = "red")           #4th week
#lines(xts[673:735], col = 'yellow')
legend(0,6, c("Week 1", "Week 2", "Week 3", "Week 4"), col = c("black", "blue", "green", "red"), lty = c(1,1,1,1), bty ="n")

#___________________ ACF and PACF _____________________________#

acf(xts, lag.max = 200, plot = TRUE)
abline(v=seq(0,200,24), lty = c(2,2,2,2,2,2,2,2,2))

acf(lnxts, lag.max = 200, plot = TRUE)
abline(v=seq(0,200,24), lty = c(2,2,2,2,2,2,2,2,2))

pacf(xts,lag.max = 200, plot = TRUE)
abline(v=seq(0,200,24), lty = c(2,2,2,2,2,2,2,2,2))

pacf(lnxts,lag.max = 200, plot = TRUE)
abline(v=seq(0,200,24), lty = c(2,2,2,2,2,2,2,2,2))



#__________ USEFUL FUNCTIONS (Tests and Plots) _____________#


# Useful Plotting Function:
diagtool <- function(residuals){
  par(mfrow=c(3,1), mar=c(3,3,1,1), mgp=c(2,0.7,0))
  plot(residuals, type="l")
  acf(residuals)
  pacf(residuals)
}


# Likelihood Ratio Test Function:
#     To test whether the more complex model is significanlty better than the simpler model
ratiotest <- function(fit1, fit2, df){
  pchisq(-2* ( fit1$loglik - fit2$loglik ), df=df, lower.tail = FALSE)
}


# F-test Function:
Ftest <- function(fit1, fit2){
  s1 <- sum(fit1$residuals^2)
  s2 <- sum(fit2$residuals^2)
  n1 <- 3
  n2 <- 4
  
  pf( (s1-s2)/(n2-n1) / (s2/(length(fit1$residuals)-n2)), df1 = n2 - n1, df2 = (length(fit1$residuals)-n2), lower.tail = FALSE)
}

# Signs test using Binomial Model:
signs <- function(analyse0){
  # sign test, mean, and sd
  res <- analyse0$residuals
  n.res <- length(res)
  (n.res-1)/2 # mean:
  sqrt((n.res-1)/4)  ### sd:
  (n.res-1)/2 + 1.96 * sqrt((n.res-1)/4) * c(-1,1) ### 95% interval:
  ### test:
  (N.sign.changes <- sum( res[-1] * res[-n.res]<0 ))
  binom.test(N.sign.changes, n.res-1)
}

# p-values of coefficients function
# gives the p-values of the hypothesis test: H0 (Null Hypothesis) coeff = 0 -vs- H1 coeff != 0 
p_values <- function(aa){
  (1-pnorm(abs(aa$coef)/sqrt(diag(aa$var.coef))))*2
}



#___________INITIAL MODEL EXPLORATION: SEEKING STATIONARITY______________#

tsdisplay(xts, lag.max = 170,
          main="Original Time Series Data - No Transform - No ARIMA applied", xlab="Time (hrs)")

tsdisplay(lnxts, lag.max = 200,
          main="Log Transformation - No ARIMA components applied", xlab="Time (hrs)")

tsdisplay(diff(lnxts,24), lag.max = 200,
          main="Seasonally differenced: (0,0,0)x(0,1,0)_24", xlab="Time(hrs)")

tsdisplay(diff(diff(lnxts,24)), lag.max = 200,
          main="Non-seasonal & Seasonally differenced: (0,1,0)x(0,1,0)_24", xlab="Time (hrs)")


# COMMENTS: The 2 differenced ACF and PACF suggest adding an ARIMA(0,0,0)(0,0,1)24 component
# on top of the differences; resulting in  an ARIMA(0,1,0)(0,1,1)24.


#__________ MODEL EXPLORATION _____________#
#       (Trying different setups)

#.......Non-Differenced Models...........#     AIC WINNERS: Models 4,5 and 14s
model0 <- arima(x = lnxts.train, order = c(3, 0, 0))
model0; p_values(model0)  # AIC=80.38, all coeff. significant

model1 <- arima(x = lnxts.train, order = c(3, 0, 0), seasonal = list(order = c(1, 0, 0), period = 24))
model1; p_values(model1)   # AIC=62, all coeff. significant                                                           

model2 <- arima(x = lnxts.train, order = c(2, 0, 0), seasonal = list(order = c(1, 0, 0), period = 24))
model2; p_values(model2)  # AIC=76.86, all coeff. significant

model3 <- arima(x = lnxts.train, order = c(3, 0, 1), seasonal = list(order = c(1, 0, 0), period = 24))
model3; p_values(model3)  # AIC=61.79, 3 coeff. not significant (remove?)

model4 <- arima(x = lnxts.train, order = c(3, 0, 1), seasonal = list(order = c(1, 0, 1), period = 24))
model4; p_values(model4)  # AIC=-9.15, 3 coeff. not significant (remove?), sar almost 1

model5 <- arima(x = lnxts.train, order = c(2, 0, 1), seasonal = list(order = c(1, 0, 1), period = 24))
model5; p_values(model5)  # lowest AIC=-10.51 of non-differenced models, all signfinicant, sar almost 1, which suggests seasonal differencing

model14 <- arima(x=lnxts.train, order=c(3,0,0), seasonal=list(order=c(1,0,2),period=24))
model14; p_values(model14) # lowest AIC=-8.18, ar2 insignificant

model14s <- arima(x=lnxts.train, order=c(3,0,0), seasonal=list(order=c(1,0,2),period=24),transform.pars = FALSE, fixed= c(NA,0,NA,NA,NA,0,NA))
model14s;                 # lowest AIC=-10.74, 

# Model Suggested by auto.arima: (2,0,1)x(0,0,0)
auto_model <- auto.arima(lnxts.train, stepwise=FALSE, approximation=FALSE)
auto_model; p_values(auto_model)   # AIC=78.22, all coeff. significant


#.......Seasonal-Differenced Models...........#     AIC WINNERS: Models 6, 6s and sD 

model6 <- arima(x = lnxts.train, order = c(3, 0, 0), seasonal = list(order = c(0, 1, 1), period = 24))
model6; p_values(model6);  # aic = 15.75, ar2 coeff. insig

model6s <- arima(x=lnxts.train, order=c(3,0,0), seasonal=list(order=c(0,1,1),period=24),transform.pars = FALSE, fixed= c(NA,0,NA,NA))
model6s;                   # aic = 14.69   # REMOVING INSIGNIFICANT AR2 coefficient

model7 <- arima(x = lnxts.train, order = c(2, 0, 0), seasonal = list(order = c(0, 1, 1), period = 24))
model7; p_values(model7)  # aic = 19.42, all sig, 

model8 <- arima(x = lnxts.train, order = c(2, 0, 1), seasonal = list(order = c(1, 1, 1), period = 24))
model8; p_values(model8)  # aic = 17.81, 1 coeff. insig, 

model10 <- arima(x = lnxts.train, order = c(2, 0, 0), seasonal = list(order = c(0, 2, 1), period = 24)) 
model0; p_values(model10) # aic = 80.38, all sig

model11 <- arima(x = lnxts.train, order = c(2, 0, 1), seasonal = list(order = c(0, 2, 1), period = 24)) 
model11; p_values(model11) # aic = 482.96, all sig

# Model Suggested by auto.arima: (2,0,4)x(0,1,0)
auto_model_sD <- auto.arima(lnxts.train, lambda=0, d=0, D=1, max.order=9, stepwise=FALSE, approximation=FALSE)
auto_model_sD; p_values(auto_model_sD)   # AIC= -1160.05!!! , 1 coeff. insignificant (ma3)


#.......Seasonal-Differenced + Non-Seasonal Differenced Models...........#     AIC WINNERS: Models 13 and sDd 

model9 <- arima(x = lnxts.train, order = c(2, 1, 1), seasonal = list(order = c(1, 1, 1), period = 24))
model9; p_values(model9)  # aic = 67.39, ar2 coeff insignificant 

model13 <- arima(x = lnxts.train, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 24)) 
model13; p_values(model13) # aic = 61.95, all sig

model12 <- arima(x = lnxts.train, order = c(0, 1, 0), seasonal = list(order = c(0, 1, 1), period = 24)) 
model12; p_values(model12) # aic = 64.99, all sig

# Model Suggested by auto.arima: (5,1,0)x(0,1,0)
auto_model_sDd <- auto.arima(lnxts.train, lambda=0, d=1, D=1, max.order=9, stepwise=FALSE, approximation=FALSE)
auto_model_sDd; p_values(auto_model_sDd)   # AIC=-1096.57 !!! , 1 coeff. insignificant (ar3)



#________________________ RESIDUAL ANALYSIS _________________________#
#          (We filter out models that dont pass this test)
#    (Then we compute further detailed tests in the next section)

tsdisplay(residuals(model0))
tsdisplay(residuals(model1))
tsdisplay(residuals(model2))
tsdisplay(residuals(model3))
tsdisplay(residuals(model4), lag.max = 60) # Good! Looks almost like white noise
    tsdiag(model4, gof.lag = 60)
tsdisplay(residuals(model5), lag.max = 60) # Looks the same as above
    tsdiag(model5, gof.lag = 60)
tsdisplay(residuals(model6), lag.max = 60) # Great! Nothing outside the bounds
    tsdiag(model6, gof.lag = 60)
tsdisplay(residuals(model7), lag.max = 60)
tsdisplay(residuals(model8), lag.max = 60)
tsdisplay(residuals(model9))
tsdisplay(residuals(model10), lag.max = 60) # 24 lag out of bounds, 
    tsdiag(model10, gof.lag = 60)
tsdisplay(residuals(model11), lag.max = 60)
    tsdiag(model11, gof.lag = 60)
tsdisplay(residuals(model12))
tsdisplay(residuals(auto_model))
tsdisplay(residuals(auto_model_sD))
tsdisplay(residuals(auto_model_sDd))
    tsdiag(auto_model_sDd, gof.lag = 60)


#______________ Model Diagnostics & Tests ___________________# 
#      (on models that passed initial inspection)    


# QQ Plot
qqnorm(model4$residuals)
qqline(model4$residuals,col=2)

qqnorm(model5$residuals)
qqline(model5$residuals,col=2)

qqnorm(model6$residuals)
qqline(model6$residuals,col=2)

qqnorm(model10$residuals)
qqline(model10$residuals,col=2)

qqnorm(model11$residuals)
qqline(model11$residuals,col=2)

# Signs (+ -) Autocorrelation Test
signs(model4)
signs(model5)
signs(model6)
signs(model10)
signs(model11)

# Cumulative Periodogram  PLot
#     - do residuals look like white noise? 
cpgram(lnxts.train)
cpgram(model4$residuals) # yes!
cpgram(model5$residuals) # yes!
cpgram(model6$residuals) # yes!
cpgram(model10$residuals) # yes!
cpgram(model11$residuals) # yes!

par(mfrow=c(2,2))
qqnorm(model5$residuals)
qqline(model5$residuals,col=2)

qqnorm(model10$residuals)
qqline(model10$residuals,col=2)

cpgram(model5$residuals) # yes!
cpgram(model10$residuals) # yes!
dev.off()

# Ratio Likelihood Tests: Is the more complex models significantly better? A: Nope!
ratiotest(model5, model4, 1)    #not sig. p =  0.4206363
ratiotest(model10, model11, 1)  #not sig. p =  0.07512814

#_____________________RMSE: COMPARING ALL MODELS_______________________#

getrmse <- function(x,h,fit)
{
  train.end <- time(x)[length(x)-h]
  test.start <- time(x)[length(x)-h+1]
  train <- window(x,end=train.end)
  test <- window(x,start=test.start)
  fc <- forecast(fit,h=h)
  return(accuracy(fc,test)[2,"RMSE"]) #the smaller the RMSE, the better prediction performance
}


#..... Non-differenced Models.........#
getrmse(lnxts,h=48, model1)   #model1: 1.026942
getrmse(lnxts,h=48, model2)   #model2: 1.028274
getrmse(lnxts,h=48, model3)   #model3: 1.02536
getrmse(lnxts,h=48, model4)   #model4: 0.9664464  Great!
getrmse(lnxts,h=48, model5)   #model5: 0.9654734  Great as well!   BEST!
getrmse(lnxts,h=48, auto_model)  #model_auto: 1.031515  
getrmse(lnxts,h=48, model14)  #model_sigs:  0.9674714
getrmse(lnxts,h=48, model14s) # 0.9674714

#.....Seasonal Differenced Models.......#
getrmse(lnxts,h=48, model6)   #model6: 0.9678455
getrmse(lnxts,h=48, model6s)  #model6s 0.9679578
getrmse(lnxts,h=48, model7)   #model7: 0.967771
getrmse(lnxts,h=48, model8)   #model8: 0.9673714
getrmse(lnxts,h=48, model10)  #model10: 0.8925458                 BEST! 
getrmse(lnxts,h=48, model11)  #model11: 0.888723                  BESTEST EVER!
getrmse(lnxts,h=48, auto_model_sD)  #model_sD: 0.9739474

#.....Seasonal & Non-seasonal Differenced Models.......#
getrmse(lnxts,h=48, model9)   #model9: 0.9735253
getrmse(lnxts,h=48, model13)  #model13: 0.9710542 
getrmse(lnxts,h=48, model12)  #model12: 0.9866931  
getrmse(lnxts,h=48, auto_model_sDd)  #model_sDd: 1.314824  worst!

#_______________________________________________________________***  RMSE WINNERS: Models 5, 10 and 11 ***_____


#_____________Forecasts_________________#
pred_model4 <-forecast.Arima(model4, h=48, level = 95) 
pred_model11 <- forecast.Arima(model11, h=48, level = 95) 
pred_model6 <-forecast.Arima(model6, h=48, level = 95) 
pred_model14 <- forecast.Arima(model14, h=48, level = 95) 
pred_model14s <-forecast.Arima(model14s, h=48, level = 95) 
pred_model5 <- forecast.Arima(model5, h=48, level = 95) 
pred_model10 <-forecast.Arima(model10, h=48, level = 95) 
pred_model11 <- forecast.Arima(model11, h=48, level = 95)

plot(pred_model4, ylab="NOx", xlab="Time(hrs) Model 4")       #model4
lines(688:735, lnxts[688:735], col = "green")

plot(forecast(model6, h = 48), ylab="NOx", xlab="Time(hrs) Model 6")  #model6
lines(688:735, lnxts[688:735], col = "green")

plot(forecast(model14), ylab="NOx", xlab="Time(hrs) Model 14") #model14
lines(688:735, lnxts[688:735], col = "green")

plot(forecast(model14s), ylab="NOx", xlab="Time(hrs) Model 14s") #model14s
lines(688:735, lnxts[688:735], col = "green")

plot(forecast(model5), ylab="NOx", xlab="Time(hrs) Model 5")  #model5
lines(688:735, lnxts[688:735], col = "green")

plot(forecast(model10), ylab="NOx", xlab="Time(hrs) Model 10") #model10    #BEST!
lines(688:735, lnxts[688:735], col = "green")

plot(forecast(model11), ylab="NOx", xlab="Time(hrs) Model 11") #model11    #AND BEST!
lines(688:735, lnxts[688:735], col = "green")


# Transform forecasts back to original scale 
mean4 = exp(pred_model4$mean)
lower4 = exp(pred_model4$lower)
upper4 = exp(pred_model4$upper)

mean5 = exp(pred_model5$mean)
lower5 = exp(pred_model5$lower)
upper5 = exp(pred_model5$upper)

mean6 = exp(pred_model6$mean)
lower6 = exp(pred_model6$lower)
upper6 = exp(pred_model6$upper)

mean10 = exp(pred_model10$mean)
lower10 = exp(pred_model10$lower)
upper10 = exp(pred_model10$upper)

mean11 = exp(pred_model11$mean)
lower11 = exp(pred_model11$lower)
upper11 = exp(pred_model11$upper)

mean14s = exp(pred_model14s$mean)
lower14s = exp(pred_model14s$lower)
upper14s = exp(pred_model14s$upper)


# PLOT IN ORIGINAL SCALE

# Plot in original scale: MOdel 4
plot(xts[1:687], ylab='NOx', xlab ="Time (Hrs): MODEL 4: (3,0,1)x(1,0,1): AIC=-9.15  RMSE = 0.9664464 ", type = "l", xlim=c(0,750))
lines(688:735, xts[688:735], col = "green")
lines(688:735, mean4, col ='red')
lines(688:735, lower4, col ='orange')
lines(688:735, upper4, col ='orange')
legend(0,110, c("Training", "Test (for predictions)", "Predictions", "95% CI"), col = c("black", "green", "red", "orange"), lty = c(1,1,1,1), bty ="n")


# Plot in original scale: MOdel 5
plot(xts[1:687], ylab='NOx', xlab ="Time (Hrs): MODEL 5: (2,0,1)x(1,0,1): AIC=-10.51  RMSE = 0.9654734  ", type = "l", xlim=c(0,750))
lines(688:735, xts[688:735], col = "green")
lines(688:735, mean5, col ='red')
lines(688:735, lower5, col ='orange')
lines(688:735, upper5, col ='orange')
legend(0,110, c("Training", "Test (for predictions)", "Predictions", "95% CI"), col = c("black", "green", "red", "orange"), lty = c(1,1,1,1), bty ="n")

# Plot in original scale: MOdel 6
plot(xts[1:687], ylab='NOx', xlab ="Time (Hrs): MODEL 6: (3,0,1)x(0,1,1): AIC=14.69  RMSE = 0.9678455 ", type = "l", xlim=c(0,750))
lines(688:735, xts[688:735], col = "green")
lines(688:735, mean6, col ='red')
lines(688:735, lower6, col ='orange')
lines(688:735, upper6, col ='orange')
legend(0,110, c("Training", "Test (for predictions)", "Predictions", "95% CI"), col = c("black", "green", "red", "orange"), lty = c(1,1,1,1), bty ="n")

# Plot in original scale: Model 10
plot(xts[1:687], ylab='NOx', xlab ="Time (Hrs): MODEL 10: (2,0,0)x(0,2,1): AIC=80.38 RMSE = 0.8925458  ", type = "l", xlim=c(0,750))
lines(688:735, xts[688:735], col = "green")
lines(688:735, mean10, col ='red')
lines(688:735, lower10, col ='orange')
lines(688:735, upper10, col ='orange')
legend(0,110, c("Training", "Test (for predictions)", "Predictions", "95% CI"), col = c("black", "green", "red", "orange"), lty = c(1,1,1,1), bty ="n")


# Plot in original scale: Model 11
plot(xts[1:687], ylab='NOx', xlab ="Time (Hrs): MODEL 11: (2,0,1)x(0,2,1): AIC=482.96 RMSE = 0.888723   ", type = "l", xlim=c(0,750))
lines(688:735, xts[688:735], col = "green")
lines(688:735, mean11, col ='red')
lines(688:735, lower11, col ='orange')
lines(688:735, upper11, col ='orange')
legend(0,110, c("Training", "Test (for predictions)", "Predictions", "95% CI"), col = c("black", "green", "red", "orange"), lty = c(1,1,1,1), bty ="n")


# Plot in original scale: MOdel 14s
plot(xts[1:687], ylab='NOx', xlab ="Time (Hrs): MODEL 14s: (3,0,0)x(1,0,2): AIC=-8.18 RMSE = 0.9674714  ", type = "l", xlim=c(0,750))
lines(688:735, xts[688:735], col = "green")
lines(688:735, mean14s, col ='red')
lines(688:735, lower14s, col ='orange')
lines(688:735, upper14s, col ='orange')
legend(0,110, c("Training", "Test (for predictions)", "Predictions", "95% CI"), col = c("black", "green", "red", "orange"), lty = c(1,1,1,1), bty ="n")

# POINT PREDICTIONS (1,24,48 hr)

mean5[1]; mean5[24]; mean5[48]
lower5[1]; lower5[24]; lower5[48]
upper5[1]; upper5[24]; upper5[48];


