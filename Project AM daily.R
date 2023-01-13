# https://fred.stlouisfed.org/series/LNS11300002
# Deforestation Warnings Brazil
library(astsa)
library(tseries)
library(forecast)
library(tseries)
library(timeSeries)
library(tsbox)
library(lubridate)
library(dplyr)
library(nortest)
library(TSA)
library(fGarch)


dev.off()
df <- read.csv("/Users/mpign/OneDrive/Documentos/12-2020/Spring USU/Econometrics II/Books and data/R and data/Project/deter-amz-daily-3-25-2022-2_08_27 AM com filtro.csv",header = TRUE)
df
head(df)
df$ï..viewDate<-mdy(df$ï..viewDate) # date format

head(df)
str(df)

# Sum the deforestation daily
df2<-df %>%
  group_by(ï..viewDate) %>%
  summarise(Tot.deforestation = sum(areaMunKm))

head(df2)

df3<-as.data.frame(df2,row.names = 1)
View(df3)
head(df3)
str(df3)

plot(df3,type="l")

#transform in time series
df4<-ts_ts(ts_long(df3))
head(df4)
str(df4)
plot(df4)

str(df4)

# organize ts
Wrate<-ts(df3$Tot.deforestation,start=c(2015,213),frequency = 365)
print(Wrate)
plot(Wrate)
str(Wrate)
install.packages("imputeTS")
library(imputeTS)
Wrate <- na_interpolation(Wrate)
sum(is.na(Wrate))
#Wrate<-na.remove(Wrate) # exclude NAs. Problem?
#Wrate=ts
plot(Wrate) #lost some days

##
#Wrate<-substituteNA(df4, type = "zeros")
#Wrate<-na.omit(df4, method = "z")
#str(Wrate)
#Wrate1<-ts_ts(ts_long(Wrate))
#str(Wrate1)
#Wrate<-ts(Wrate,start=c(2015,213),frequency = 365)
#str(Wrate)
#plot(Wrate)
#Wrate3<-ts_ts(ts_long(Wrate2))
#str(Wrate3)
#Wrate1<-ts_ts(ts_long(Wrate))
#str(Wrate1)
#plot(Wrate1,type="l")
##


plot(Wrate,ylab="Daily Deforestation Warnings")
abline(reg=lm(Wrate~time(Wrate)))

dev.off()
hist(Wrate)
time(Wrate)
par(mfrow=c(2,1))
acf(Wrate)
pacf(Wrate)
acf2(Wrate)
adf.test(Wrate,alternative="stationary",k=0)
adf.test(Wrate,alternative="stationary")

y1<-window(Wrate,end=c(2021,6))
y2<-window(Wrate,start=c(2021,7))
str(y1)


dev.off()

## trying the log
#plot(log(y1))
#dWrate<-diff(y1,differences=1,lag=1)
#plot(dWrate)
#acf2
#acf(dWrate)
#pacf(dWrate)
##

# estimating model: best ARIMA(4,1,1)
mWrate<-auto.arima(y1,ic="bic",seasonal=TRUE,stationary=FALSE,trace=TRUE)
print (mWrate)
predict(mWrate, n.head)
str(mWrate)
z1<-mWrate$coef/sqrt(diag(mWrate$var.coef))
print(z1)
p1<-2*(1-pnorm(abs(mWrate$coef)/sqrt(diag(mWrate$var.coef))))
print(p1)               
confint(mWrate,level=0.95)  
str(mWrate)

#model vs actual data
dev.off()
plot(mWrate$x)
lines(mWrate$fitted,col="red")

# Diagnostics

#In-sample accuracy
is<-accuracy(mWrate)
print(is)
#str(mWrate)

#Q-statistics - Autocorrelation test

tsdiag(mWrate,gof.lag=36)
LB1<-Box.test(mWrate$residuals,lag=36,type="Ljung")
LB1

#ly<-log(y1)
#lmWrate<-auto.arima(ly,ic="bic",seasonal=TRUE,stationary=FALSE,trace=TRUE)
#acf2(lmWrate$residuals)


acf2(mWrate$residuals)
acf(mWrate$residuals)
pacf(mWrate$residuals)

qqnorm(mWrate$residuals)
qqline(mWrate$residuals)
shapiro.test(mWrate$residuals)
jarque.bera.test(mWrate$residuals)

ad.test(mWrate$residuals)
cvm.test(mWrate$residuals)
lillie.test(mWrate$residuals)


McLeod.Li.test(mWrate,y1)

# Forecasting 
m.fcast<-forecast(mWrate,h=300)
plot(m.fcast)
lines(y2)

dev.off()

# Out of sample accuracy
os<-accuracy(m.fcast,y2)
os

##Problem: data cant estimate Holt Winters
Hm3<-HoltWinters(y1)
plot(forecast(Hm3,h=5))

Hm3 <- tbats(y1)
fc <- forecast(Hm3, h=55)
plot(fc)

#Hm3
#str(Hm3)
#print(Hm3$fitted)
#print(c(Hm3$alpha,Hm3$beta,Hm3$gamma))
plot(Hm3,col=1:2)
legend("topleft", c("data","fitted"), col=c(1,2), lty=1, bty="n")
str(Hm3)
His<-accuracy(Hm3$fitted.values,y1)
His
par(1,2)
is
His

##
tsdiag(mWrate,gof.lag=36)

LB1<-Box.test(mWrate$residuals,lag=36,type="Ljung")
LB1

LB2<-Box.test(mWrate$residuals^2,lag=36,type="Ljung")
LB2
McLeod.Li.test(mWrate,m3$residuals^2)



vm1<-garchFit(~arma(5,2)+garch(1,0),data=diff(y1),cond.dist="norm")

vm2<-garchFit(~arma(4,1)+garch(2,0),data=diff(y1),cond.dist="norm",trace=FALSE)

vm5<-garchFit(~arma(4,1)+garch(5,0),diff(y1),cond.dist="norm",trace=FALSE)

vm8<-garchFit(~arma(4,1)+garch(8,0),diff(y1),cond.dist="norm",trace=FALSE)

vm10<-garchFit(~arma(4,1)+garch(10,0),diff(y1),cond.dist="norm",trace=FALSE)

vm11<-garchFit(~arma(4,1)+garch(1,1),diff(y1),cond.dist="norm",trace=TRUE)

vm12<-garchFit(~arma(4,1)+garch(1,2),diff(y1),cond.dist="norm",trace=FALSE)

vm21<-garchFit(~arma(4,1)+garch(2,1),diff(y1),cond.dist="norm",trace=FALSE)

vm55<-garchFit(~arma(4,1)+garch(5,5),diff(y1),cond.dist="norm",trace=FALSE)

vm85<-garchFit(~arma(4,1)+garch(8,5),diff(y1),cond.dist="norm",trace=FALSE)

vm101<-garchFit(~arma(4,1)+garch(10,1),diff(y1),cond.dist="norm",trace=FALSE)

vm1010<-garchFit(~arma(4,1)+garch(10,10),diff(y1),cond.dist="norm",trace=FALSE)

# str(vm1)

vm1@fit$ics

vm2@fit$ics

vm5@fit$ics

vm8@fit$ics

vm10@fit$ics


vm11@fit$ics

vm12@fit$ics

vm21@fit$ics

vm55@fit$ics

vm85@fit$ics

vm101@fit$ics

vm1010@fit$ics

#Residuals
#sres1<-ts(vm11@residuals/vm11@sigma.t,start=c(2015,213),end = c(2021,6),frequency = 365)
sres1<-vm11@residuals/vm11@sigma.t
sres12<-sres1^2
#str(vm1)

# diagnostis residuals

plot(sres1,type="l")
plot(sres12,type="l")

par(mfrow=c(2,2))
acf(sres1)
pacf(sres1)
acf(sres12)
pacf(sres12)
dev.off()

LBvm1<-Box.test(sres12,lag=12,type="Ljung")
LBvm1

LBvm12<-Box.test(sres1,lag=12,type="Ljung")
LBvm12

qqnorm(sres1, main="Oil ARIMA+GARCH Residual Quantile plot")
qqline(sres1)


qqnorm(sres12, main="Oil ARIMA+GARCH Residual Quantile plot")
qqline(sres12)

# Forescast

vm11.predict<-predict(vm1,n.ahead=300,plot=TRUE)
str(vm11)

forecast(vm11,h=5)


#getMethod("predict","fGARCH")
#str(vm1)
#str(vm1.predict)

vm11m<-ts(vm11.predict$meanForecast,start=c(2015,213),end = c(2021,6),frequency = 365)
#str(vm1m)
vm11u<-ts(vm11.predict$meanForecast+2*vm11.predict$standardDeviation,start = c(2010, 11),end = c(2010, 25),frequency = 52)
vm11d<-ts(vm11.predict$meanForecast-2*vm11.predict$standardDeviation,start = c(2010, 11),end = c(2010, 25),frequency = 52)

#dev.off()
#dev.set(dev.next())
#par(mfrow=c(2,2))
#while (!is.null(dev.list()))  dev.off()

#plot(diff(y1))
#polygon(x=c(index(vm11u),rev(index(vm11u))),y=c(vm11d,rev(vm11u)),col="slategray1",border="slategray4")
#lines(diff(y2),col="red")
#lines(vm11m,col="blue")

##
#dfr<-ts(df4,start=c(2015,213),end=c(2022,69),frequency = 365)
#plot(dfr)
#str(dfr)
#library(zoo)
#z <- read.zoo(df3)
#plot(z, type = "l")
#head(z)
#library(lubridate)
#df4<-as.Date(df3$ï..viewDate,"%d/%m/%Y")
#head(df4)
###
#a_vec <- rep(0, (vm11))      
#u2 <- length(ar)      
#a_vec[1] = ar[1] + ma[1]       
#if ((n.ahead - 1) > 1) {
#  for (i in 2:(n.ahead - 1)) {
#    a_vec[i] <- ar[1:min(u2, i - 1)] * a_vec[(i - 1):(i - u2)]             
#  }
#}
####


