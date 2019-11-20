library(gamlss)
library(gamlss.add)
library(mgcv)
library(RevoUtilsMath)

s1data<-read.table("C:/Data/mouse23/2014-09-29_10-21-24/begin1/CS/rtable.csv",sep=",",header=TRUE)
s1data$f10 <- scale(s1data$f10,center = FALSE, scale = TRUE)
s1data$f24 <- scale(s1data$f24,center = FALSE, scale = TRUE)
s1data$f80 <- scale(s1data$f80,center = FALSE, scale = TRUE)
s1data$Running_Speed <- abs(s1data$Running_Speed)
s1data$Session <- factor(rep("one",length(s1data$Time)))

s2data<-read.table("C:/Data/mouse23/2014-09-29_10-21-24/begin2/CS/rtable.csv",sep=",",header=TRUE)
s2data$f10 <- scale(s2data$f10,center = FALSE, scale = TRUE)
s2data$f24 <- scale(s2data$f24,center = FALSE, scale = TRUE)
s2data$f80 <- scale(s2data$f80,center = FALSE, scale = TRUE)
s2data$Running_Speed <- abs(s2data$Running_Speed)
s2data$Session <- factor(rep("two",length(s2data$Time)))

s3data<-read.table("C:/Data/mouse23/2014-09-29_10-21-24/begin3/CS/rtable.csv",sep=",",header=TRUE)
s3data$f10 <- scale(s3data$f10,center = FALSE, scale = TRUE)
s3data$f24 <- scale(s3data$f24,center = FALSE, scale = TRUE)
s3data$f80 <- scale(s3data$f80,center = FALSE, scale = TRUE)
s3data$Running_Speed <- abs(s3data$Running_Speed)
s3data$Session <- factor(rep("three",length(s3data$Time)))

mydata <- rbind(s1data,s2data,s3data)
rm(s1data,s2data,s3data)
mydata <- mydata[mydata$Running_Speed <= 30,]

mf24s1<-gamlss(f80~1,sigma.formula=~1,
               nu.formula=~ga(~s(Time,by=Session,bs="cr")+s(Running_Speed,bs="cr")+s(Theta_Phase,bs="cc")+s(Position,bs="cc")
                              +ti(Time,Running_Speed,by=Session,bs="cr")+ti(Time,Theta_Phase,by=Session,bs=c("cr","cc"))+ti(Time,Position,by=Session,bs=c("cr","cc"))
                              +ti(Running_Speed,Theta_Phase,bs=c("cr","cc"))+ti(Running_Speed,Position,bs=c("cr","cc"))+ti(Theta_Phase,Position,bs=c("cc","cc")),method="REML"),
               tau.formula=~1, data=mydata,family="SEP3",method=RS(200),control=gamlss.control(c.crit=0.5))

mns <- gam(f80 ~ s(Time,by=Session)+s(Running_Speed)+s(Theta_Phase,bs="cc")+s(Position,bs="cc")
          +ti(Time,Running_Speed,by=Session)+ti(Time,Theta_Phase,by=Session,bs=c("tp","cc"))+ti(Time,Position,by=Session,bs=c("tp","cc"))
          +ti(Running_Speed,Theta_Phase,bs=c("tp","cc"))+ti(Running_Speed,Position,bs=c("tp","cc"))+ti(Theta_Phase,Position,bs=c("cc","cc")),method="REML",select=TRUE
          family=scat(theta=NULL,link="identity"),data = mydata, control = list(trace=TRUE))

plot(mf80s1,ts=TRUE)
wp(mf80s1,ylim.all=0.6)
plot(getSmo(mf80s1,what="mu"),pages=1,scale=0)
plot(getSmo(mf80s1,what="sigma"),pages=1,scale=0)
plot(getSmo(mf80s1,what="nu"),pages=1,scale=0)
plot(getSmo(mf80s1,what="tau"),pages=1,scale=0)
