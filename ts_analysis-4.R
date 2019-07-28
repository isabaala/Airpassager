require(forecast)
require(itsmr)
require(tseries)
traffic=cdg

attach(traffic)
str(traffic)

#rescale CDG_ORLY
head(kpassenger)

#create the time series
ts.traffic=ts(data = kpassenger, start = c(2000,01), frequency = 12)

#plot the time series
plot(ts.traffic, ylab="numbers of air passenger in thousands of millions",
     xlab="Years", main="Evolution of monthly air traffic from 2000 to 2018", 
     xlim=c(2000, 2020), ylim=c(4500, 11000))

plot(aggregate(ts.traffic), ylab="air traffic (in thousands of passengers)", main="Overall Trend")
boxplot(ts.traffic~cycle(ts.traffic), names=c("Jan", "Feb","Mar","Apr","May","Jun","Jul","Aug",
                                              "Sept","Oct","Nov","Dec"), 
        xlab=" ", ylab="air traffic (in thousands of passengers)", main="Seasonal Variations")


#show the seasonality, we only plot the 1-60lags
acf(ts.traffic, lag.max = 60)

#remove trend and seasonality
#the model is additive so we'll go first for a linear regression usig OLS
#we suppose a linear trend so we'll try with a simple polynomial regression format
n <- length(ts.traffic)
t <- 1:n
seasonality <- 12

#create our seasonal dummy variables
seasonal.dummy<-rep(1:seasonality , times = n/seasonality) 

# ----------------------------Test with the Raw data (no log-transform)------------------------------------------

#OLS outputs
#the linear trend is statistically significant and different from 0, thus
#there exists a time effect.
reg.t <- lm(ts.traffic ~ t + as.factor(seasonal.dummy))
summary(reg.t)

#we can check for the quadratic trend :
t2 <- t^2
reg.t2 <- lm(ts.traffic ~ t + t2 + as.factor(seasonal.dummy))
summary(reg.t2)

#cubic trend :
#adding a cubic trend shows that the coefficients is not very significant
t3 <- t^3
reg.t3 <-lm(ts.traffic ~ t + t2 + t3 + as.factor(seasonal.dummy))
summary(reg.t3)


#plot the residuals : we will plot the residuals only for the linear and quadratic trend
#the residuals still show a seasonal pattern
res.t2 <- reg.t2$residuals
plot(res.t2, main = "Residuals of the quadratic regression", ylab="Residuals")
lines(res.t2)
abline(h=0, col="red", lwd=2)
#test the stationarity
Box.test(x = res.t2,lag = 20,type = "Ljung-Box")
adf.test(res.t2)
#------------------------Dealing with the data log-transformed---------------------------------------


#Since the polynomial term in the trend was significant, it could suggest that
#the data could be exponentially function of the time. We can fit a linear model
#with the log of the data :

log.ts <- log(ts.traffic)

plot(log.ts, ylab="log(airtraffic)", main="log-transformed air passenger traffic",xlim=c(2000, 2020))

log.reg.t <- lm(log.ts ~ t + as.factor(seasonal.dummy))
summary(log.reg.t)

#check the residuals :
log.res <- log.reg.t$residuals
plotc(log.res )
Box.test(x = log.res,lag = 20,type = "Ljung-Box")
adf.test(log.res)

#It is still not better so we will go for a moveing average filter (brockwell-davis method) in order
#to get rid of the seasonality more efficiently. 
bd.model <- as.numeric(log(ts.traffic))
n.bd <- length(bd.model)

ma.bd <- bd.model
for(j in 7:(n.bd-6)) {
  ma.bd[j] <- (0.5*bd.model[j-6] + bd.model[j-5] + bd.model[j-4] + bd.model[j-3] + bd.model[j-2] + bd.model[j-1] 
               + bd.model[j] + bd.model[j+1]+bd.model[j+2]+bd.model[j+3]+bd.model[j+4]
               +bd.model[j+5]+0.5*bd.model[j+6])/12
}

plot(bd.model[7:222], main="Air passenger traffic, with adjusted trend")
lines(bd.model[7:222])
lines(ma.bd[7:222],col="red", lwd=2)

#same coding using the filter(): check if our results are consists
wg<-c(.5, rep(1,11),.5)/12
tr.es<-filter(bd.model,filter=wg, sides=2)
plotc(bd.model)
lines(tr.es, col="blue")

xtilde <- matrix(bd.model - ma.bd,ncol=12,byrow=TRUE)
w <- rep(0,12)
w[1] <- mean( xtilde[-1,1]  )
w[2] <- mean( xtilde[-1,2]  )
w[3] <- mean( xtilde[-1,3]  )
w[4] <- mean( xtilde[-1,4]  )
w[5] <- mean( xtilde[-1,5]  )
w[6] <- mean( xtilde[-1,6]  )
w[7] <- mean( xtilde[-n.bd,7]  )
w[8] <- mean( xtilde[-n.bd,8]  )
w[9] <- mean( xtilde[-n.bd,9]  )
w[10] <- mean( xtilde[-n.bd,10]  )
w[11] <- mean( xtilde[-n.bd,11]  )
w[12] <- mean( xtilde[-n.bd,12]  )
w <- w - mean(w)
w <- rep(w,times = n.bd/12)
bd.model.ns <- bd.model - w
plotc(bd.model.ns)

#Comparing the coefficient to STL decomposition:
stl.ts <- stl(log.ts, s.window = "periodic")
stl.ts$time.series[1:12, 1:3]

cbind(stl.ts$time.series[1:12, 1:3],w[1:12])

#Now the seasonality is removed, we will go for a differentiation method in order
#to get rid of the trend. 

#we differentiate once to remove the trend
dYt<-diff(bd.model.ns)

#then we check if the TS is centered around 0, if not we'll centering it
mean(dYt)
mu <- mean(dYt)
dYt <- dYt - mu
mean(dYt)

#The plot shows that there maybe some seasonality left
plotc(dYt)

# ---------------------------Check for stationarity------------------------------------

#tests for stationarity
#Except the Box test all three test suggests that we have a stationary time series
Box.test(dYt, lag=24, type = "Ljung-Box")
adf.test(dYt, alternative = "stationary")
pp.test(dYt, alternative = "stationary")
kpss.test(dYt, null = "L")

#estimating the parameter of the model : acf/pacf
acf(dYt) # more or less exponential decrease, with only2 or 3 significant spikes.
pacf(dYt) #1 significant spikes and maybe one other

# --------------------------Fit the model1---------------------------------------------
auto.arima(dYt, d=0, D=0, max.p=4, max.q=2, max.P=0, max.Q=2, max.order=6, stepwise = T,
           trace=T, seasonal = T)
fit1<-Arima(y=dYt, order=c(4,0,1), seasonal = c(0,0,0), include.constant = F)
fit1
fit1$coef
fit1$sigma2

plotc(fit1$residuals) # we see a little seasonality remaining
acf(fit1$residuals) #good

#the model appears to be valid, we continue with the forecasting
tsdiag(fit1)

# -------------------------------Forescasting------------------------------------------------

#prediction of the stationary time series
h <- 60
mypred2 <- predict(object=fit1,n.ahead = h)
n<-length(dYt)

plot(x=t[-1],y=dYt,xlim = c(0,n+h),type="l",lwd=2 , xlab="t")
lines(x=n:(n+h) , y = c(dYt[n-1],mypred2$pred) , col="red" , lwd=2 )

#prediction of the actual time series
n<-length(log.ts)
mypred_mu <- mypred2$pred+mu 
mypred_cum <- cumsum(mypred_mu)+w[1:h] 
finalpred2 <- bd.model.ns[n]+mypred_cum 

plot(x = 1:n , y = log.ts , type="l" , lwd = 2 , 
     xlim=c(80,n+h) , ylim = c(8.5,9.5) , xlab="t" , ylab="log Air Traffic")
lines(x = n:(n+h) , y = c(log.ts[n],finalpred2) , 
      col = "red" , lwd = 2 , lty = 1) # new prediction


# --------------------- Comparaison with another model-------------------------------------
#another model
ft<-auto.arima(log.ts)
tsdiag(auto.arima(log.ts))


