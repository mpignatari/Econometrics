# Term project R Group1 


library(tidyverse)
library(stargazer)
library(car)
library(effects)
library(corrgram)
library(wooldridge)
library(lmtest)

###############################################################################
# Data : We are going to use the dataset gpa2!

data('gpa2')
help(gpa2)
range(gpa2$sat)
summary(gpa2)

# If you have Na's, you have to clean the dataset first by using the following code. 

df <- gpa2
df <- df[complete.cases(df),]

# if you want to exclude some variables, you can use select function from dplyr package. 

df2 <- select(df, c(sat,athlete,hsperc,female,hsize,black,white,hsrank))

head(df2)
summary(df2)
row(df2)
nrow(df2)
# Find the most influential variables by eyeballing the correlogram
df2<-df2[,c(2,1,3:8)]
corrgram(df2, order=FALSE, lower.panel=panel.shade,upper.panel=panel.pie, text.panel=panel.txt,main="Correlations between variables")

# let's look at the structure of our data set:
str(df2)

# Model 1
MRM1<-lm(sat~.,df2)
count(df2,black)
bptest(MRM1)
hist(resid(MRM1), freq = FALSE)
lines(density(resid(MRM1), adjust = 2), lty = "dotted", col="red", lwd=3) 
coeftest(MRM1)
coeftest(MRM1, vcov= hccm(MRM1, type="hc0"))
vif(MRM1)
# Let's try log transformation and see if that works. 
MRM_log <- lm(log(sat) ~ athlete+hsperc+female+hsize+black+white+hsrank,df)
stargazer(MRM1,MRM_log, type = "text")
hist(resid(MRM_log), freq = FALSE)
lines(density(resid(MRM_log), adjust = 2), lty = "dotted", col="red", lwd=3) 
bptest(MRM_log) # still bad news! we reject the null: Homoskedastic 

# Model 2
MRM_inter <- lm(sat ~ athlete+hsperc+female+hsize+black+white+hsrank+athlete:black,df)
stargazer(MRM_inter, type = "text")
bptest(MRM_inter)
vif(MRM_inter)

# Model 3
MRM_inter2 <- lm(sat ~ athlete+hsperc+female+hsize+black+white+hsrank+athlete:black+athlete:female,df)
stargazer(MRM_log_inter, type = "text") 
bptest(MRM_inter2)
vif(MRM_inter2)

# Model 4
MRM_log_inter <- lm(log(sat) ~ athlete+hsperc+female+hsize+black+white+hsrank+athlete:black+athlete:female,df)
stargazer(MRM_log_inter, type = "text")
bptest(MRM_log_inter)
vif(MRM_log_inter) 

stargazer(MRM1,MRM_inter,MRM_inter2,MRM_log_inter,df2,type="text")

###############################################################################
### Engines

### Robust standard errors
coeftest(MRM1)
coeftest(MRM1, vcov= hccm(MRM_log, type="hc0"))

coeftest(MRM_inter)
coeftest(MRM_inter, vcov= hccm(MRM_inter, type="hc0"))

coeftest(MRM_inter2)
coeftest(MRM_inter2, vcov= hccm(MRM_inter2, type="hc0"))


coeftest(MRM_log_inter)
coeftest(MRM_log_inter, vcov= hccm(MRM_log_inter, type="hc0"))

###LMP Model and WLS
summary(df2)
# constructing the dummy dependent variable: EAA (earning above average)
df3<- mutate(df2, EEA=ifelse((sat>mean(sat)),1,0))

#Model 1
LPM <- lm(EEA ~ .,df3)
bptest(LPM)

y_hat <- predict(LPM)
summary(y_hat) # this suggests that LPM is not the best approach. 
h     <- y_hat * (1-y_hat)
range(h) # So we need to force h>0
h<- ifelse(h<0,0.01,h)
summary(h)
w<- 1/h
LPM_wls <- lm(EEA ~., weights = w ,df3)

#Model 2
LPM2 <- lm(EEA ~ athlete+hsperc+female+hsize+black+white+hsrank+athlete:black,df3)
bptest(LPM2)

y_hat <- predict(LPM2)
summary(y_hat) # this suggests that LPM is not the best approach. 
h     <- y_hat * (1-y_hat)
range(h) # So we need to force h>0
h<- ifelse(h<0,0.01,h)
summary(h)
w<- 1/h
LPM_wls2 <- lm(EEA ~athlete+hsperc+female+hsize+black+white+hsrank+athlete:black,weights = w ,df3)


#Model 3
LPM3 <- lm(EEA ~ athlete+hsperc+female+hsize+black+white+hsrank+athlete:black+athlete:female,df3)
bptest(LPM3)

y_hat <- predict(LPM3)
summary(y_hat) # this suggests that LPM is not the best approach. 
h     <- y_hat * (1-y_hat)
range(h) # So we need to force h>0
h<- ifelse(h<0,0.01,h)
summary(h)
w<- 1/h
LPM_wls3 <- lm(EEA ~athlete+hsperc+female+hsize+black+white+hsrank+athlete:black+athlete:female,weights = w ,df3)


stargazer(LPM_wls,LPM_wls2,LPM_wls3,type = 'text')


### FGLS

# Model 1
uhat1 <- resid(MRM1)
logu1 <- log(uhat1^2)
reg_step1 <- lm(logu1~athlete+hsperc+female+hsize+black+white+hsrank,df)
ghat1 <- predict(reg_step1)
hhat1 <- exp(ghat1)
w1<- 1/hhat1
MRM1_FGLS <- lm(sat ~ athlete+hsperc+female+hsize+black+white+hsrank, weights = w1, df )
stargazer(MRM1_FGLS,type="text")
bptest(MRM1_FGLS)


#log
MRM_log <- lm(log(sat) ~ athlete+hsperc+female+hsize+black+white+hsrank,df)

uhatl1 <- resid(MRM_log)
logul1 <- log(uhatl1^2)
reg_stepl1 <- lm(logul1~athlete+hsperc+female+hsize+black+white+hsrank,df)
ghatl1 <- predict(reg_stepl1)
hhatl1 <- exp(ghatl1)
wl1<- 1/hhatl1
MRM_log_FGLS <- lm(log(sat) ~ athlete+hsperc+female+hsize+black+white+hsrank, weights = wl1, df )
stargazer(MRM_log_FGLS,type="text")
bptest(MRM_log_FGLS)

MRM_log <- lm(log(sat) ~ athlete+hsperc+female+hsize+black+white+hsrank,df)
stargazer(MRM1,MRM_log, type = "text")
hist(resid(MRM_log), freq = FALSE)
lines(density(resid(MRM_log), adjust = 2), lty = "dotted", col="red", lwd=3) 
bptest(MRM_log) # still bad news! we reject the null: Homoskedastic 

stargazer(MRM1_FGLS,MRM_log_FGLS, type = "text")

# Model 2
uhat2 <- resid(MRM_inter)
logu2 <- log(uhat2^2)
reg_step2 <- lm(logu2~athlete+hsperc+female+hsize+black+white+hsrank+athlete:black,df)
ghat2 <- predict(reg_step2)
hhat2 <- exp(ghat2)
w2<- 1/hhat2
MRM_inter_FGLS <- lm(sat ~ athlete+hsperc+female+hsize+black+white+hsrank+athlete:black, weights = w2, df )
stargazer(MRM_inter_FGLS,type="text")
bptest(MRM_inter_FGLS)

##Model 3
uhat3 <- resid(MRM_inter2)
logu3 <- log(uhat3^2)
reg_step3 <- lm(logu3~athlete+hsperc+female+hsize+black+white+hsrank+athlete:black+athlete:female,df)
ghat3 <- predict(reg_step3)
hhat3 <- exp(ghat3)
w3<- 1/hhat3
MRM_inter2_FGLS <- lm(sat ~ athlete+hsperc+female+hsize+black+white+hsrank+athlete:black+athlete:female, weights = w3, df )
stargazer(MRM_inter2_FGLS,type="text")
bptest(MRM_FGLS)

# Model 4
uhat4 <- resid(MRM_log_inter)
logu4 <- log(uhat4^2)
reg_step4 <- lm(logu4~athlete+hsperc+female+hsize+black+white+hsrank+athlete:black+athlete:female,df)
ghat4 <- predict(reg_step4)
hhat4 <- exp(ghat4)
w4<- 1/hhat4
MRM_log_inter_FGLS <- lm(log(sat) ~ athlete+hsperc+female+hsize+black+white+hsrank+athlete:black+athlete:female, weights = w3, df )
stargazer(LPM_wls3,MRM_inter2_FGLS,MRM_log_inter_FGLS,type="text")


stargazer(LPM_wls3,MRM_inter2_FGLS,MRM_log_inter_FGLS,type="text")

a<-coeftest(MRM_log_inter, vcov= hccm(MRM_log_inter, type="hc0"))
stargazer(MRM1_FGLS,MRM_inter_FGLS,MRM_inter2_FGLS,MRM_log_inter_FGLS,type = 'text',column.labels = c("level", "1 inter","2 inter","log 2 inter"))

###########################################################################
-----------------
qqnorm(resid(MRM1), pch = 1, frame = FALSE)
qqline(resid(MRM1), col = "red", lwd = 2)
qqnorm(resid(MRM_log), pch = 1, frame = FALSE)
qqline(resid(MRM_log), col = "red", lwd = 2)
###############################################################################
# Hypothesis testing 

# testing if gender matters? 

# using our Heteroskedastic log model
linearHypothesis(MRM_log_inter, c("female=0", "athlete:female=0"), vcov=hccm(MRM_log_inter,type="hc0")) # Ok, we reject the null. So it seems that gender matters. 

##############################################################################

# some additional codes FYI
# Redefining dummy variables

df <- mutate(happiness, Ishappy=ifelse((happy=="very happy"|happy=="pretty happy"),1,0) , income_gt25k=ifelse((income=="$25000 or more"),1,0), Isreligious=ifelse((attend=="never"),0,1))
df <- df[,c("Ishappy","income_gt25k", "Isreligious", "prestige","divorce","educ","babies","preteen","teens","mothfath16","black","female","blackfemale","unem10")]
corrgram(df, order=TRUE, lower.panel=panel.shade,upper.panel=panel.pie, text.panel=panel.txt,main="Correlations between variables")










