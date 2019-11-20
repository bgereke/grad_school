#Licking
library(mgcv)
library(ggplot2)
library(stats)
library(pracma)

setwd("Z:/imaging/wheel_run_due/")

#load all data
files <- c("md229_04062017_wheelrundue_2017_04_06_04_47_06_986.txt",
           "md231_04062017_wheelrundue_2017_04_06_06_09_57_310.txt")
nd <- 1 #number of days per mouse
nm <- length(files)/nd #number of mice
d <- 1
m <- 1
sdf <- 2500 #smoothing spline degrees of freedom for speed estimation
tl <- 384 #track length
DATA<-read.table(files[1],sep=",",header=TRUE)
# DATA$beep_on <- NULL
names(DATA)[names(DATA)=="ard_timestamp"] <- "time"
DATA$mouse<-as.factor(m)
DATA$day<-as.factor(d)
m <- m + 1
#scale time to seconds
DATA$time <- (DATA$time - min(DATA$time))/1000
#interpolate position
time <- DATA$time[which(diff(DATA$position)>0)+1]
position <- DATA$position[which(diff(DATA$position)>0)+1]
sf<-splinefun(time,position,method="monoH.FC") 
DATA$position <- sf(DATA$time)
#get speed and acceleration
sm <- sm<-smooth.spline(DATA$time,sf(DATA$time),df=sdf,all.knots=TRUE)
DATA$speed <- predict(sm,DATA$time,deriv=1)$y
DATA$acc <- predict(sm,DATA$time,deriv=2)$y
#remove edges
DATA <- DATA[DATA$time>=min(time) & DATA$time<=max(time),]

#get licks
pks <- findpeaks(DATA$lick_voltage,zero="-")
pks <- pks[pks[,1]>3,]
DATA$licks <- 0
DATA$licks[pks[,2]] <- 1
#create history terms 
wmax<-80
History <- embed(DATA$licks,wmax+1)
History <- History[,2:ncol(History)]
#get time since reward
DATA$rtime <- DATA$time
rt <- DATA$time[which(diff(DATA$reward)>0)+1]
for(t in 1:nrow(DATA)){
  DATA$rtime[t] <- DATA$time[t] - max(rt[rt<DATA$time[t]])
}
DATA$rtime[is.infinite(DATA$rtime)] <- NA
DATA$rtime[is.na(DATA$rtime)] <- DATA$time[is.na(DATA$rtime)]

DATA <- DATA[-c(1:wmax),]
#concatenate with remaining files
for(t in 2:length(files)){
  data<-read.table(files[t],sep=",",header=TRUE)
  # data$beep_on <- NULL
  names(data)[names(data)=="ard_timestamp"] <- "time"
  data$mouse<-as.factor(m)
  data$day<-as.factor(d)
  if(m==nm)(d<-d+1)
  if(m<nm){m<-m+1}else(m<-1)
  
  #scale time to seconds
  data$time <- (data$time - min(data$time))/1000
  
  #interpolate position
  time <- data$time[which(diff(data$position)>0)+1]
  position <- data$position[which(diff(data$position)>0)+1]
  sf<-splinefun(time,position,method="monoH.FC") 
  data$position <- sf(data$time)
  
  #get speed and acceleration
  sm <- sm<-smooth.spline(data$time,sf(data$time),df=sdf,all.knots=TRUE)
  data$speed <- predict(sm,data$time,deriv=1)$y
  data$acc <- predict(sm,data$time,deriv=2)$y
  
  #remove edges
  data <- data[data$time>=min(time) & data$time<=max(time),]
  
  #get licks
  pks <- findpeaks(data$lick_voltage,zero="-")
  pks <- pks[pks[,1]>3,]
  data$licks <- 0
  data$licks[pks[,2]] <- 1
  
  #create history terms for each day
  lags <- embed(data$licks,wmax+1)
  lags <- lags[,2:ncol(lags)]
  History <- rbind(History,lags)
  
  #get time since reward
  data$rtime <- data$time
  rt <- data$time[which(diff(data$reward)>0)+1]
  for(t in 1:nrow(data)){
    data$rtime[t] <- data$time[t] - max(rt[rt<data$time[t]])
  }
  data$rtime[is.infinite(data$rtime)] <- NA
  data$rtime[is.na(data$rtime)] <- data$time[is.na(data$rtime)]
  
  data <- data[-c(1:wmax),]
  DATA <- rbind(DATA, data)
}

DATA$speed <- DATA$speed*60/64 #rpm
DATA$acc <- DATA$acc/64*60^2 #rpm^2
DATA$MouseDay<-with(DATA,interaction(mouse,day))
days <- unique(DATA$day)
History<-as.matrix(History)
Hbins<-matrix(1:wmax,nrow(History),ncol(History),byrow = TRUE)

DATA$rtime <- log(DATA$rtime)

ridge<-diag(wmax)
Lags <- matrix(0,nrow = nrow(History),ncol = nd*wmax)
Lags[DATA$day=="1",1:wmax] <- History[DATA$day=="1",]
Lags[DATA$day=="2",(wmax+1):(2*wmax)] <- History[DATA$day=="2",]
Lags[DATA$day=="3",(2*wmax+1):(3*wmax)] <- History[DATA$day=="3",]
# DATA <- DATA[-(Lags[,1]==1|Lags[,wmax+1]==1),] 
# Lags <- Lags[-(Lags[,1]==1|Lags[,wmax+1]==1),]
# Lags <- Lags[,-c(1,wmax+1)]

#wrap-to-lap
DATA$position <- DATA$position/tl
DATA$position <- DATA$position - floor(DATA$position)
DATA$position <- tl*DATA$position

#rotate from reward location
DATA$rotpos <- DATA$position - 177
DATA$rotpos[DATA$rotpos<0] <- DATA$rotpos[DATA$rotpos<0] - min(DATA$rotpos) + 74



mod <- bam(licks ~ day + History + mouse +
                   s(rotpos,by=day,bs="cr",k=15,id=1) + 
                   s(time,by=mouse,bs="cr") + 
                   s(speed,by=mouse,bs="cr") + 
                   s(acc,by=mouse,bs="cr"),
                   family=poisson(link="log"),
                   data = DATA)

nb <- 100
pbins <- seq(min(DATA$position),max(DATA$position),length.out=nb)
tbins <- seq(min(DATA$time),max(DATA$time),length.out=nb)
rtbins <- seq(min(DATA$rtime),max(DATA$rtime),length.out=nb)
sbins <- seq(0,max(DATA$speed),length.out=nb)
abins <- seq(-5000,5000,length.out=nb)
pdat <- with(DATA,data.frame(position=rep(pbins,times=nd),time=rep(tbins,times=nd),
                             rtime=rep(rtbins,times=nd),speed=rep(sbins,times=nd),
                             acc=rep(abins,times=nd),day=rep(days,each=nb)))
pdat$Lags<-matrix(1,nrow(pdat),ncol(Lags))
pred<-predict.bam(mod, pdat, type = "terms", se.fit = TRUE)

#plot position fit
df <- data.frame(position=rep(pbins,times=nd),
                 fit=c(pred$fit[,colnames(pred$fit)=="s(position):day1"][1:nb],
                       pred$fit[,colnames(pred$fit)=="s(position):day2"][(nb+1):(2*nb)],
                       pred$fit[,colnames(pred$fit)=="s(position):day3"][(2*nb+1):(3*nb)]),
                 day=rep(days,each=nb))
se <- c(pred$se.fit[,colnames(pred$se.fit)=="s(position):day1"][1:nb],
        pred$se.fit[,colnames(pred$se.fit)=="s(position):day2"][(nb+1):(2*nb)],
        pred$se.fit[,colnames(pred$se.fit)=="s(position):day3"][(2*nb+1):(3*nb)])
df$low <- exp(df$fit - 2*se)
df$high <- exp(df$fit + 2*se)
df$fit <- exp(df$fit)
ggplot(df, aes(x = position, y = fit)) +
  geom_rect(xmin=126,xmax=135,ymin=0,ymax=2,fill="red")+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray",alpha=0.75) +
  geom_line() +
  scale_x_continuous(breaks=seq(10,190,by=20),labels=seq(1,10),name="position (zones)",expand=c(0,0))+
  scale_y_continuous(name="Lick Rate (multiplier)",expand=c(0,0))+
  facet_wrap(~day,nrow=1) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA,size=1))

#plot time fit
df <- data.frame(time=rep(tbins,times=nd),
                 fit=c(pred$fit[,colnames(pred$fit)=="s(time):day1"][1:nb],
                       pred$fit[,colnames(pred$fit)=="s(time):day2"][(nb+1):(2*nb)],
                       pred$fit[,colnames(pred$fit)=="s(time):day3"][(2*nb+1):(3*nb)]),
                 day=rep(days,each=nb))
se <- c(pred$se.fit[,colnames(pred$se.fit)=="s(time):day1"][1:nb],
        pred$se.fit[,colnames(pred$se.fit)=="s(time):day2"][(nb+1):(2*nb)],
        pred$se.fit[,colnames(pred$se.fit)=="s(time):day3"][(2*nb+1):(3*nb)])
df$low <- exp(df$fit - 2*se)
df$high <- exp(df$fit + 2*se)
df$fit <- exp(df$fit)
ggplot(df, aes(x = time, y = fit)) +
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray",alpha=0.75) +
  geom_line() +
  scale_x_continuous(breaks=seq(0,600,by=60),labels=seq(0,10),name="time (min)",expand=c(0,0))+
  scale_y_continuous(name="Lick Rate (multiplier)",expand=c(0,0))+
  facet_wrap(~day,nrow=1) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA,size=1))

#plot rtime fit
df <- data.frame(rtime=rep(rtbins,times=nd),
                 fit=c(pred$fit[,colnames(pred$fit)=="s(rtime):day1"][1:nb],
                       pred$fit[,colnames(pred$fit)=="s(rtime):day2"][(nb+1):(2*nb)],
                       pred$fit[,colnames(pred$fit)=="s(rtime):day3"][(2*nb+1):(3*nb)]),
                 day=rep(days,each=nb))
se <- c(pred$se.fit[,colnames(pred$se.fit)=="s(rtime):day1"][1:nb],
        pred$se.fit[,colnames(pred$se.fit)=="s(rtime):day2"][(nb+1):(2*nb)],
        pred$se.fit[,colnames(pred$se.fit)=="s(rtime):day3"][(2*nb+1):(3*nb)])
df$low <- exp(df$fit - 2*se)
df$high <- exp(df$fit + 2*se)
df$fit <- exp(df$fit)
ggplot(df, aes(x = rtime, y = fit)) +
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray",alpha=0.75) +
  geom_line() +
  scale_x_continuous(name="log[time since reward (sec)]",expand=c(0,0))+
  scale_y_continuous(name="Lick Rate (multiplier)",expand=c(0,0))+
  facet_wrap(~day,nrow=1) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA,size=1))

#plot speed fit
df <- data.frame(speed=rep(sbins,times=nd),
                 fit=c(pred$fit[,colnames(pred$fit)=="s(speed):day1"][1:nb],
                       pred$fit[,colnames(pred$fit)=="s(speed):day2"][(nb+1):(2*nb)],
                       pred$fit[,colnames(pred$fit)=="s(speed):day3"][(2*nb+1):(3*nb)]),
                 day=rep(days,each=nb))
se <- c(pred$se.fit[,colnames(pred$se.fit)=="s(speed):day1"][1:nb],
        pred$se.fit[,colnames(pred$se.fit)=="s(speed):day2"][(nb+1):(2*nb)],
        pred$se.fit[,colnames(pred$se.fit)=="s(speed):day3"][(2*nb+1):(3*nb)])
df$low <- exp(df$fit - 2*se)
df$high <- exp(df$fit + 2*se)
df$fit <- exp(df$fit)
ggplot(df, aes(x = speed, y = fit)) +
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray",alpha=0.75) +
  geom_line() +
  scale_x_continuous(name="speed (rpm)",expand=c(0,0))+
  scale_y_continuous(name="Lick Rate (multiplier)",expand=c(0,0))+
  facet_wrap(~day,nrow=1) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA,size=1))

#plot acc fit
df <- data.frame(acc=rep(abins,times=nd),
                 fit=c(pred$fit[,colnames(pred$fit)=="s(acc):day1"][1:nb],
                       pred$fit[,colnames(pred$fit)=="s(acc):day2"][(nb+1):(2*nb)],
                       pred$fit[,colnames(pred$fit)=="s(acc):day3"][(2*nb+1):(3*nb)]),
                 day=rep(days,each=nb))
se <- c(pred$se.fit[,colnames(pred$se.fit)=="s(acc):day1"][1:nb],
        pred$se.fit[,colnames(pred$se.fit)=="s(acc):day2"][(nb+1):(2*nb)],
        pred$se.fit[,colnames(pred$se.fit)=="s(acc):day3"][(2*nb+1):(3*nb)])
df$low <- exp(df$fit - 2*se)
df$high <- exp(df$fit + 2*se)
df$fit <- exp(df$fit)
ggplot(df, aes(x = acc, y = fit)) +
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_vline(xintercept = 0) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray",alpha=0.75) +
  geom_line() +
  scale_x_continuous(name="acceleration (rpm^2)",expand=c(0,0))+
  scale_y_continuous(name="Lick Rate (multiplier)",expand=c(0,0))+
  facet_wrap(~day,nrow=1) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA,size=1))

#plot history fit
df <- data.frame(lags=rep(seq(2*16,wmax*16,length.out=wmax-1),times=nd),
                 fit=mod$coefficients[c(3:(wmax+1),(wmax+3):(2*wmax+1),(2*wmax+3):(3*wmax+1))],
                 day=rep(days,each=wmax-1))
se <- summary(mod)$se[2:(2*wmax-1)]
df$low <- exp(df$fit - 2*se)
df$high <- exp(df$fit + 2*se)
df$fit <- exp(df$fit)
ggplot(df, aes(x = lags, y = fit)) +
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray",alpha=0.75) +
  geom_line() +
  scale_x_continuous(name="lags (ms)",expand=c(0,0))+
  scale_y_continuous(name="Lick Rate (multiplier)",expand=c(0,0),limits=c(0,4.6))+
  facet_wrap(~day,nrow=1) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA,size=1))
