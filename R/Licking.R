#Licking
library(mgcv)
library(ggplot2)

setwd("C:/Users/Brian/Documents")

DATA<-read.table("rtable_Licking.csv",sep=",",header=TRUE)
DATA$Mouse<-as.factor(DATA$Mouse)
DATA$Day<-as.factor(DATA$Day)
DATA$MouseDay<-with(DATA,interaction(Mouse,Day))
days <- unique(DATA$Day)
nd <- length(days)

#create history terms for each day
wmax<-20
History <- embed(DATA$Licks[DATA$Day==days[1]],wmax+1)
History <- History[,2:ncol(History)]
DATA <- DATA[-c(1:wmax),]
for(d in 2:nd){
  
  lags <- embed(DATA$Licks[DATA$Day==days[d]],wmax+1)
  lags <- lags[,2:ncol(lags)]
  History <- rbind(History,lags)
  
  start <- min(which(DATA$Day==days[d]))
  DATA <- DATA[-c(start:(start+wmax-1)),]
  
}

History<-as.matrix(History)
Hbins<-matrix(log((1:wmax)+5),nrow(History),ncol(History),byrow = TRUE)

History <- History[DATA$RTime>0,]
DATA <- DATA[DATA$RTime>0,]
DATA$RTime <- log(DATA$RTime)

ridge<-diag(wmax)

mod <- bam(Licks ~ Day + History + s(RTime,by=Day,bs="cr") +
                   s(Position,by=Day,bs="cr") +
                   te(Time,Speed,bs="cr"),
                   # s(Position,Mouse,bs="fs",xt="cr",m=1),
                   family=poisson(link="log"),
                   data = DATA)

rmod <- bam(mod$residuals ~ s(Position,by=Day,bs="cr"),
           # family=poisson(link="log"),
           data = DATA)

nb <- 100
pbins <- seq(min(DATA$Position),max(DATA$Position),length.out=nb)
pdat <- with(DATA,data.frame(Day = rep(days,each=nb),Position=rep(pbins,times=nd)))
pdat$History<-matrix(1,nrow(pdat),ncol(History))
pred<-predict.bam(mod, pdat, type = "iterms", se.fit = TRUE)
want1<-grep("s\\(Position\\)\\:Day1",colnames(pred$fit))
want2<-grep("s\\(Position\\)\\:Day2",colnames(pred$fit))
want3<-grep("s\\(Position\\)\\:Day3",colnames(pred$fit))
want4<-grep("s\\(Position\\)\\:Day4",colnames(pred$fit))
want5<-grep("s\\(Position\\)\\:Day5",colnames(pred$fit))
idx1 <- pdat$Day == days[1]
idx2 <- pdat$Day == days[2]
idx3 <- pdat$Day == days[3]
idx4 <- pdat$Day == days[4]
idx5 <- pdat$Day == days[5]
pdf1 <- data.frame(fit=pred$fit[idx1,want1], 
                  low=pred$fit[idx1,want1]-2*pred$se.fit[idx1,want1], 
                  high=pred$fit[idx1,want1]+2*pred$se.fit[idx1,want1])
pdf2 <- data.frame(fit=pred$fit[idx2,want2], 
                   low=pred$fit[idx2,want2]-2*pred$se.fit[idx2,want2], 
                   high=pred$fit[idx2,want2]+2*pred$se.fit[idx2,want2])
pdf3 <- data.frame(fit=pred$fit[idx3,want3], 
                   low=pred$fit[idx3,want3]-3*pred$se.fit[idx3,want3], 
                   high=pred$fit[idx3,want3]+3*pred$se.fit[idx3,want3])
pdf4 <- data.frame(fit=pred$fit[idx4,want4], 
                   low=pred$fit[idx4,want4]-4*pred$se.fit[idx4,want4], 
                   high=pred$fit[idx4,want4]+4*pred$se.fit[idx4,want4])
pdf5 <- data.frame(fit=pred$fit[idx5,want5], 
                   low=pred$fit[idx5,want5]-5*pred$se.fit[idx5,want5], 
                   high=pred$fit[idx5,want5]+5*pred$se.fit[idx5,want5])
pdf <- rbind(pdf1,pdf2,pdf3,pdf4,pdf5)
pdf$Day <- pdat$Day
pdf$Position <- pdat$Position

mu1 = mod$coefficients[1]
mu2 = mod$coefficients[1] + mod$coefficients[2]
mu3 = mod$coefficients[1] + mod$coefficients[3]
mu4 = mod$coefficients[1] + mod$coefficients[4]
mu5 = mod$coefficients[1] + mod$coefficients[5]

pdf$mu <- rep(c(mu1,mu2,mu3,mu4,mu5),each=nb)

pdf$fit <- exp(pdf$fit+pdf$mu)/median(diff(DATA$Time))
pdf$low <- exp(pdf$low+pdf$mu)/median(diff(DATA$Time))
pdf$high <- exp(pdf$high+pdf$mu)/median(diff(DATA$Time))

ggplot(pdf, aes(x = Position, y = fit)) +
  geom_rect(xmin=58,xmax=73,ymin=0.5,ymax=22,fill="red")+
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray",alpha=0.75) +
  geom_line() +
  scale_x_continuous(breaks=c(13.5,39.5,65.5,91.5),labels=c("z1","z2","z3","z4"),name="position",expand=c(0,0))+
  scale_y_continuous(name="Lick Rate (Hz)",expand=c(0,0))+
  geom_hline(aes(yintercept=1),linetype="dashed")+
  geom_vline(aes(xintercept=26),linetype="dashed")+
  geom_vline(aes(xintercept=52),linetype="dashed")+
  geom_vline(aes(xintercept=78),linetype="dashed")+
  facet_wrap(~Day,nrow=1)+theme(aspect.ratio = 1)+
  theme(aspect.ratio = 1,
        panel.background = element_blank(),aspect.ratio = 1,
        # panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(colour="black",fill=NA,size=1))

acfdf2 <- data.frame(ACF = exp(mod$coefficients[6:25]), Lags = seq(1:wmax))
acfdf <- rbind(acfdf,acfdf2)
acfdf$model <- as.factor(rep(c("AllLicks","FirstLicks"),each=wmax))
dt<-median(diff(DATA$Time))*1000

ggplot(acfdf, aes(x = Lags, y = ACF)) +
  geom_line(aes(colour=model),size=1) +
  scale_x_continuous(breaks=seq(1,wmax,by=2),labels=round(seq(dt,dt*wmax,by=2*dt)),name="lags (ms)",expand=c(0,0))+
  scale_y_continuous(name="multiplier",expand=c(0,0))+
  geom_hline(aes(yintercept=1),linetype="dashed")+
  theme(aspect.ratio = 1,
        panel.background = element_blank(),aspect.ratio = 1,
        # panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(colour="black",fill=NA,size=1))
