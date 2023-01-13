#title: "Project Econometrics"
#author: "Marcelo Pignatari"
#date: 12/15/22


#libraries
library("psych")
library("dplyr")
library("corrplot")
library("stargazer")
library("censReg")
library("sampleSelection")
library("AER")
library("truncreg")
library(dplyr)
#install.packages("npsf")
library(dplyr)
library(tseries)
library(lmtest)

#loading the data
df<- read.csv("C:/Users/mpign/OneDrive/Documentos/12-2020/Fall 2022 USU/Econometrics III/Project/rend15.csv")
#,stringsAsFactors = TRUE)
#rend 4 limit 8000
head(df)
df$X <- NULL
df$Unnamed..0 <- NULL
df$Unnamed..0.1 <- NULL
df$Unnamed..0.1.1 <- NULL
head(df)

#naming the variables
summary(df)
#describe(df)

df <- rename(df, gender = 'V0302_4')
df <- rename(df, age = 'V8005')
df <- rename(df, indigenous = 'V0404_0')
df <- rename(df, black = 'V0404_4')
df <- rename(df, yellow = 'V0404_6')
df <- rename(df, brown = 'V0404_8')
df <- rename(df, educ = 'V4803')
df <- rename(df, rend = 'V4719')
names(df)
#head(df)
#df$logrend<-log((df$rend)+1)

#summary stats
stargazer(df,type="text")

#histogram
hist(df$rend,freq=FALSE)
lines(density(df$rend))

#OLS 
OLS<-lm(rend~gender+age+indigenous+black+yellow+brown+educ+marr,data=df)
qqnorm(OLS$residuals)
qqline(OLS$residuals)
#bptest
bptest(OLS)

#FGLS
uhat1 <- resid(OLS)
logu1 <- log(uhat1^2)
reg_step1 <- lm(logu1~gender+age+indigenous+black+yellow+brown+educ+marr,df)
ghat1 <- predict(reg_step1)
hhat1 <- exp(ghat1)
w1<- 1/hhat1
MRM1_FGLS <- lm(rend~gender+age+indigenous+black+yellow+brown+educ+marr, weights = w1, df )
stargazer(MRM1_FGLS,type="text")
#diagnostics tests
bptest(MRM1_FGLS)
qqnorm(MRM1_FGLS$residuals)
qqline(MRM1_FGLS$residuals)
plot(MRM1_FGLS$residuals)
vif(MRM1_FGLS)
jarque.bera.test(residuals(MRM1_FGLS))

#Tobit 
Tobit<-censReg(rend~gender+age+indigenous+black+yellow+brown+educ+marr,data=df,left=0)
str(Tobit)
mg1<-margEff(Tobit)
mgT1<-summary(margEff(Tobit))
mgT1
stargazer(mgT1,type="text")
#getting the residuals
y<-df$rend
x<-cbind(df$gender,df$age,df$indigenous,df$black,df$yellow,df$brown,df$educ,df$marr)
y_hat<-(b0+b1*x[,1]+b2*x[,2]+b3*x[,3]+b4*x[,4]+b5*x[,5]+b6*x[,6]+b7*x[,7]+b8*x[,8])*pnorm((b0+b1*x[,1]+b2*x[,2]+b3*x[,3]+b4*x[,4]+b5*x[,5]+b6*x[,6]+b7*x[,7]+b8*x[,8])/sigma)+sigma*dnorm((b0+b1*x[,1]+b2*x[,2]+b3*x[,3]+b4*x[,4]+b5*x[,5]+b6*x[,6]+b7*x[,7]+b8*x[,8])/sigma)
resT<- y-y_hat
#diagnostic tests
bptest(Tobit)
qqnorm(resT)
qqline(resT)
vif(Tobit)
jarque.bera.test(resT)
#attempt to get the standard errors mannualy: Error: cannot allocate vector of size 81.0 Gb
#X<-model.matrix(Tobit)
#vce <- solve(t(X) %*% X) %*% (t(X) %*% ((sigma)*diag(104281)) %*% X) %*% solve(t(X) %*% X)

#Heckman Two-Step
rend<-df$rend
gender<-df$gender
age<-df$age
indigenous<-df$indigenous
black<-df$black
yellow<-df$yellow
brown<-df$brown
educ<-df$educ
marr<-df$marr
inlf<-df$lfp
logrend<-df$logrend
#STEP 1
probit<-glm(inlf~gender+age+indigenous+black+yellow+brown+educ+marr, family=binomial(link="probit"),data = df)
#STEP 2
lambda<-dnorm(cbind(1,gender,age,indigenous,black,yellow,brown,educ,marr)%*%(probit$coef))/pnorm(cbind(1,gender,age,indigenous,black,yellow,brown,educ,marr)%*%(probit$coef))
twostep<-lm(rend[rend>0]~gender[rend>0]+age[rend>0]+indigenous[rend>0]+black[rend>0]+yellow[rend>0]+brown[rend>0]+educ[rend>0]+marr[rend>0]+lambda[rend>0])
stargazer(twostep, header=FALSE, type="text",column.sep.width = "3pt")
#Diagnostic tests
bptest(twostep)
vif(twostep)
qqnorm(residuals(twostep))
qqline(residuals(twostep))
jarque.bera.test(twostep$residuals)


#Semi-parametric
x<-cbind(df$gender,df$black,df$brown,df$educ)
clad.fn<-function(param,y,x){
  b0<-param[1]
  b1<-param[2]
  b2<-param[3]
  b3<-param[4]
  b4<-param[5]
  n=length(df$rend)
  fn=0
  for(i in 1:n){
    dev=(df$rend[i]-max(b0+b1*x[i,1]+b2*x[i,2]+b3*x[i,3]+b4*x[i,4],0))
    fn=fn+abs(dev)
  }
  return((1/n)*fn)
}
# MINIMIZATION ROUTINE
clad.es<-nlm(f=clad.fn,
             p=c(0,0,0,0,0),
             y=df$rend,
             x=x,
             hessian=TRUE)
print(clad.es)
# RESULTS
est<-clad.es$estimate
cov<-solve(clad.es$hessian)
se<-sqrt(diag(cov))
zv<-est/se
summ<-cbind(est, se, zv)
colnames(summ)<-c('Estimate', 'Std Error', 'z value')
print(summ)
#Diagnostic Tests: Error
#bptest(clad.es)
#vif(clad.es)

#Graphs
titleName <- expression("FGLS Normal Q-Q Plot")

par(mfrow = c(1, 3))
qqnorm(MRM1_FGLS$residuals, main = titleName)
qqline(MRM1_FGLS$residuals)

titleName <- expression("Tobit Normal Q-Q Plot")
qqnorm(resT, main = titleName)
qqline(resT)

titleName <- expression("Heckman Normal Q-Q Plot")
qqnorm(twostep$residuals, main = titleName)
qqline(twostep$residuals)
