toydata <- read.table("C:/Users/Brian/Desktop/toydata3.txt", quote="\"", comment.char="")

#standardize response
mu1 = mean(toydata$V1)
mu2 = mean(toydata$V2)
mu3 = mean(toydata$V3)
toydata$V1 = (toydata$V1)/sd(toydata$V1)*sign(mu1)
toydata$V2 = (toydata$V2)/sd(toydata$V2)*sign(mu2)
toydata$V3 = (toydata$V3)/sd(toydata$V3)*sign(mu3)
#set min time to 0
toydata$V4 = toydata$V4 - min(toydata$V4)
#rename vars
names(toydata)[names(toydata)=="V1"] <- "Th"
names(toydata)[names(toydata)=="V2"] <- "SG"
names(toydata)[names(toydata)=="V3"] <- "FG"
names(toydata)[names(toydata)=="V4"] <- "Time"
names(toydata)[names(toydata)=="V5"] <- "V"
names(toydata)[names(toydata)=="V6"] <- "TP"
Lags = embed(toydata$ImCS,16)
Lags = as.data.frame(Lags[,c(2:16)])
#Lags <- Lags[,-c(1:4)]
toydata <- toydata[-c(1:15),]
toydata <- cbind.data.frame(toydata,Lags)
#create segmenting var for ARMA model
#nseg = 20;
#lseg = floor(length(toydata$ImCS)/nseg)
#Seg = rep(1:nseg,each = lseg)
#nleft = length(toydata$ImCS)-nseg*lseg
#Seg = c(Seg,rep(20,nleft))
#toydata$Seg <- Seg
#run naive model and get summary
require("mgcv")
mm <- gam(ImCS ~ s(ThetaPhase, bs = "cc")+s(Velocity,k=20)+s(Time,k=20),family=scat(theta=NULL,link="identity"),
          data = toydata)
m1 <- gam(ImCS ~ s(ThetaPhase, bs = "cc")+s(Velocity)+s(Position, bs = "cc")+s(Time,k=20),family=scat(theta=NULL,link="identity"),
         data = toydata)
m2 <- gam(ImCS ~ s(ThetaPhase, Velocity)+s(Time),family=scat(theta=NULL,link="identity"),
          data = toydata)
m2 <- gam(ImCS ~ s(ThetaPhase, bs = "cc")+s(Velocity)+s(Time)+V1*V2,
          data = toydata,method = "REML")
m3 <- gam(ImCS ~ s(ThetaPhase, bs = "cc")+s(Velocity)+s(Time)+(V1+V2+V3)^3,
          data = toydata,method = "REML")
m4 <- gam(ImCS ~ s(ThetaPhase, bs = "cc")+s(Velocity)+s(Time)+(V1+V2+V3+V4)^3,
          data = toydata,method = "REML")
m5 <- gam(ImCS ~ s(ThetaPhase, bs = "cc")+s(Velocity)+s(Time)+(V1+V2+V3+V4+V5)^3,
          data = toydata,method = "REML")
m6 <- gam(ImCS ~ s(ThetaPhase, bs = "cc")+s(Velocity)+s(Time)+(V1+V2+V3+V4+V5+V6)^3,
          data = toydata,method = "REML")
m7 <- gam(ImCS ~ s(ThetaPhase, bs = "cc")+s(Velocity)+s(Time)+(V1+V2+V3+V4+V5+V6+V7)^3,
          data = toydata,method = "REML")
m8 <- gam(ImCS ~ s(ThetaPhase, bs = "cc")+s(Velocity)+s(Time)+(V1+V2+V3+V4+V5+V6+V7+V8)^3,
          data = toydata,method = "REML")
m9 <- gam(ImCS ~ s(ThetaPhase, bs = "cc")+s(Velocity)+s(Time)+(V1+V2+V3+V4+V5+V6+V7+V8+V9)^3,
          data = toydata,method = "REML")
m10 <- gam(ImCS ~ s(ThetaPhase, bs = "cc")+s(Velocity)+s(Time)+(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10)^3,
          data = toydata,method = "REML")
m11 <- gam(ImCS ~ s(ThetaPhase, bs = "cc")+s(Velocity)+s(Time)+(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11)^3,
           data = toydata,method = "REML")
m12 <- gam(ImCS ~ s(ThetaPhase, bs = "cc")+s(Velocity)+s(Time)+(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12)^3,
           data = toydata,method = "REML")
summary
m14 <- gam(ImCS ~ s(ThetaPhase, bs = "cc")+s(Velocity)+s(Time)+(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14)^3,
           data = toydata,method = "REML")
m15 <- gam(ImCS ~ s(ThetaPhase, bs = "cc")+s(Velocity)+s(Time)+(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15)^3,
           data = toydata,method = "REML")
summary(mm)
#inspect smooths
layout(matrix(1:4, nrow = 2, ncol = 2))
plot(mm,scale = 0)
layout(1)
#inspect residual autocorrelation
layout(matrix(1:2, ncol = 2))
acf(resid(fmod), lag.max = 36, main = "ACF")
pacf(resid(fmod), lag.max = 36, main = "pACF")
layout(1)
#run ARMA1-3 gamms
ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")
## AR(1)
m1 <- gamm(ImCS ~ s(ThetaPhase, bs = "cc")+ s(Position, bs = "cc") + s(Time,bs="cr")+ s(Velocity,bs="cr"),
           data = toydata,correlation = corExp(form=~Time,nugget=TRUE))
## AR(2)
m2 <- gamm(ImCS ~ s(ThetaPhase, bs = "cc") + s(Position, bs = "cc") + s(Time,bs="cr")+ s(Velocity,bs="cr"),
           data = toydata,correlation = corARMA(form = ~ 1, p = 2),method = "REML",control = ctrl)
## AR(3)
m3 <- gamm(ImCS ~ s(ThetaPhase, bs = "cc") + s(Position, bs = "cc") + s(Time,bs="cr")+ s(Velocity,bs="cr"),
           data = toydata,correlation = corARMA(form = ~ 1, p = 3),method = "REML",control = ctrl)
#compare models
anova.gam(m1,m)

#BIC
n = length(toydata$SG)
m1<-gamlss(SG~pb(T,method="GAIC",k=log(n))+pb(V,method="GAIC",k=log(n))+cy(TP,method="GAIC",k=log(n)),
           sigma.formula=~pb(T,method="GAIC",k=log(n))+pb(V,method="GAIC",k=log(n))+cy(TP,method="GAIC",k=log(n)),
           nu.formula=~pb(T,method="GAIC",k=log(n))+pb(V,method="GAIC",k=log(n))+cy(TP,method="GAIC",k=log(n)),
           tau.formula=~pb(T,method="GAIC",k=log(n))+pb(V,method="GAIC",k=log(n))+cy(TP,method="GAIC",k=log(n)),
           data=toydata,family="JSU",control=gamlss.control(c.crit=0.1))
plot(m1,ts=TRUE)
wp(m1,ylim.all=0.6)
term.plot(m1,what="mu",pages=1,ylim="free",ask=FALSE)
term.plot(m1,what="sigma",pages=1,ylim="free",ask=FALSE)
term.plot(m1,what="nu",pages=1,ylim="free",ask=FALSE)
term.plot(m1,what="tau",pages=1,ylim="free",ask=FALSE)
#AIC
m2<-gamlss(SG~pbz(Time)+pbz(V)+cy(TP),
           sigma.formula=~pbz(Time)+pbz(V)+cy(TP),
           nu.formula=~pbz(Time)+pbz(V)+cy(TP),
           tau.formula=~pbz(Time)+pbz(V)+cy(TP),
           data=toydata,family="JSU",method=RS(80),control=gamlss.control(c.crit=0.1))
plot(m2,ts=TRUE)
wp(m2,ylim.all=0.6)
term.plot(m2,what="mu",pages=1,ylim="free",ask=FALSE)
term.plot(m2,what="sigma",pages=1,ylim="free",ask=FALSE)
term.plot(m2,what="nu",pages=1,ylim="free",ask=FALSE)
term.plot(m2,what="tau",pages=1,ylim="free",ask=FALSE)
#AIC with lags
m2<-gamlss(SG~pbz(Time)+pbz(V)+cy(TP)+la(SG,lags=30,from.lag=1),
           sigma.formula=~pbz(Time)+pbz(V)+cy(TP)+la(SG,lags=30,from.lag=1),
           nu.formula=~pbz(Time)+pbz(V)+cy(TP)+la(SG,lags=30,from.lag=1),
           tau.formula=~pbz(Time)+pbz(V)+cy(TP)+la(SG,lags=30,from.lag=1),
           data=toydata,family="SEP3",control=gamlss.control(c.crit=0.1),method=RS(80))
plot(m2,ts=TRUE)
wp(m2,ylim.all=0.6)
term.plot(m2,what="mu",pages=1,ylim="free",ask=FALSE)
term.plot(m2,what="sigma",pages=1,ylim="free",ask=FALSE)
term.plot(m2,what="nu",pages=1,ylim="free",ask=FALSE)
term.plot(m2,what="tau",pages=1,ylim="free",ask=FALSE)
#mgcv
m3<-gamlss(SG~ga(~s(Time)+s(V)+s(TP,bs="cc"),method="REML",select=TRUE),
           sigma.formula=~ga(~s(Time)+s(V)+s(TP,bs="cc"),method="REML",select=TRUE),
           nu.formula=~ga(~s(Time)+s(V)+s(TP,bs="cc"),method="REML",select=TRUE),
           tau.formula=~ga(~s(Time)+s(V)+s(TP,bs="cc"),method="REML",select=TRUE),
           data=toydata,family="SEP3",method=RS(80),control=gamlss.control(c.crit=0.1))
plot(m3,ts=TRUE)
wp(m3,ylim.all=0.6)
term.plot(m3,what="mu",pages=1,ylim="free",ask=FALSE)
term.plot(m3,what="sigma",pages=1,ylim="free",ask=FALSE)
term.plot(m3,what="nu",pages=1,ylim="free",ask=FALSE)
term.plot(m3,what="tau",pages=1,ylim="free",ask=FALSE)
#mixed model
m2<-gamlss(SG~pb(Time)+pb(V)+cy(TP)+re(random=~1,correlation=corAR1()),
           sigma.formula=~pb(Time)+pb(V)+cy(TP),
           nu.formula=~pb(Time)+pb(V)+cy(TP),
           tau.formula=~pb(Time)+pb(V)+cy(TP),
           data=toydata,family="JSU",method=RS(80),control=gamlss.control(c.crit=0.1))
plot(m2,ts=TRUE)
wp(m2,ylim.all=0.6)
term.plot(m2,what="mu",pages=1,ylim="free",ask=FALSE)
term.plot(m2,what="sigma",pages=1,ylim="free",ask=FALSE)
term.plot(m2,what="nu",pages=1,ylim="free",ask=FALSE)
term.plot(m2,what="tau",pages=1,ylim="free",ask=FALSE)