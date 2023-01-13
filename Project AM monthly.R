# https://fred.stlouisfed.org/series/LNS11300002
# Deforestation Warnings Brazil monthly
library(astsa)
library(tseries)
library(forecast)
library(dplyr)

dev.off()
df <- read.csv("/Users/mpign/OneDrive/Documentos/12-2020/Spring USU/Econometrics II/Books and data/R and data/Project/deter-amz-aggregated-3-25-2022-1_00_36 AM DesmatREal.csv",header = TRUE)
head(df)
str(df)
df<-df$Date <- paste(df$ï..year, "-", df$month)
head(df)

df[1:2] <- list(NULL) 
head(df)
df[2:4] <- list(NULL) 
head(df)
df2<-df[,c(2,1)]
head(df2)

df2$Date<-ym(df2$Date)

head(df2)
str(df2)


df3<-df2 %>%
  group_by(Date) %>%
  summarise(Tot.deforestation = sum(area))

head(df3)

dev.off()

df3<-as.data.frame(df2,row.names = 1)
View(df3)
head(df3)
str(df3)

plot(df3,type="l")

library(tsbox)
df4<-ts_ts(ts_long(df3))
str(df4)
head(df4)
plot(df4)

----
  
  dfr<-ts(df4,start=c(2015,213),end=c(2022,69),frequency = 365)
plot(dfr)
str(dfr)


library(zoo)
z <- read.zoo(df3)

plot(z, type = "l")
head(z)

library(lubridate)
df4<-as.Date(df3$ï..viewDate,"%d/%m/%Y")
head(df4)

----
  plot(ts(df4))

Wrate<-ts(df4,start=c(2015,8),frequency = 12)
Wrate<-na.remove(Wrate)

plot(Wrate,ylab="Daily Deforestation Warnings")
abline(reg=lm(Wrate~time(Wrate))


dev.off()
hist(Wrate)
time(Wrate)
par(mfrow=c(2,1))
acf(Wrate)
pacf(Wrate)
acf2(Wrate)
acf2(Wrate)
adf.test(Wrate,alternative="stationary",k=0)
adf.test(Wrate,alternative="stationary")

y1<-window(Wrate,end=c(2021,6))
y2<-window(Wrate,start=c(2021,7))


dev.off()

plot(log(y1))

dWrate<-diff(y1,differences=1,lag=1)
plot(dWrate)
acf2
acf(dWrate)
pacf(dWrate)

mWrate<-auto.arima(y1,ic="bic",seasonal=TRUE,stationary=FALSE,trace=TRUE)
print (mWrate)
str(mWrate)
z1<-mWrate$coef/sqrt(diag(mWrate$var.coef))
print(z1)
p1<-2*(1-pnorm(abs(mWrate$coef)/sqrt(diag(mWrate$var.coef))))
print(p1)               
confint(mWrate,level=0.95)   

# Diagnostics

#In-sample accuracy
is<-accuracy(mWrate)
print(is)
str(mWrate)

#Q-statistics - Autocorrelation test

tsdiag(mWrate,gof.lag=36)
LB1<-Box.test(mWrate$residuals,lag=36,type="Ljung")
LB1

ly<-log(y1)
lmWrate<-auto.arima(ly,ic="bic",seasonal=TRUE,stationary=FALSE,trace=TRUE)
acf2(lmWrate$residuals)

acf2(mWrate$residuals)
qqnorm(mWrate$residuals)
qqline(mWrate$residuals)
shapiro.test(mWrate$residuals)
jarque.bera.test(mWrate$residuals)

library(nortest)
ad.test(mWrate$residuals)
cvm.test(mWrate$residuals)
lillie.test(mWrate$residuals)

library(TSA)
McLeod.Li.test(mWrate,y1)

# Forecasting 
m.fcast<-forecast(mWrate,h=9)
plot(m.fcast)
lines(Wrate)

# Out of sample accuracy
os<-accuracy(m.fcast,y2)
os

       
       