library(mgcv)
library(MASS)
library(ggplot2)

setwd("G:/Rats/PCtable/")

#load data
atdata <- read.table("spktable.csv",sep=",",header=TRUE)
hist <- as.matrix(read.table("histtable.csv",sep=",",header=TRUE))


atdata$burstpos<-ordered(atdata$burstpos)
atdata$burstlen<-ordered(atdata$burstlen)
atdata$burstlen[atdata$burstlen>"5"]<-"5"
atdata$burstpos[atdata$burstpos>"5"]<-"5"
atdata$cell <- ordered(atdata$cell)
atdata$day <- as.factor(atdata$day) 
atdata$mouse <- as.factor(atdata$mouse)
atdata$session <- as.factor(atdata$session)
atdata$preisi = log10(atdata$preisi)
atdata$postisi = log10(atdata$postisi)
# atdata<-atdata[atdata$vel>0.1,]
atdata$daymouse<-with(atdata,interaction(day,mouse))

atcells <- unique(atdata$cell)
nc <- length(atcells)
for(c in 1:nc){
  # if(sum(data$cell==atcells[c])==0){
  #   atdata<-atdata[atdata$cell!=atcells[c],]
  #   next
  # }
  if(sum(atdata$cell==atcells[c])<300){
    hist<-hist[atdata$cell!=atcells[c],]
    atdata<-atdata[atdata$cell!=atcells[c],]
    next
  }
  if(sum(atdata$cell==atcells[c]&atdata$burstpos=="2")==0){
    hist<-hist[atdata$cell!=atcells[c],]
    atdata<-atdata[atdata$cell!=atcells[c],]
    next
  }
  if(sum(atdata$cell==atcells[c]&atdata$burstpos=="3")==0){
    hist<-hist[atdata$cell!=atcells[c],]
    atdata<-atdata[atdata$cell!=atcells[c],]
    next
  }
}

histbins <- 10^seq(log10(0.001+0.01),log10(4+0.01),length.out=51)-0.01
histbins <- histbins[-length(histbins)]
hbins<-1:ncol(hist)
hbins <- matrix(rep(histbins,times=nrow(hist)),nrow=nrow(hist),byrow=TRUE)


#replicate Quirk et al. Figure 2A
bdat1 <- data.frame(burstpos = 1:4, fit = rep(1,times=4),se = rep(0,time=4))
bdat1$fit[1] <- mean(atdata$namp[atdata$burstlen==1 & atdata$burstpos==1 & atdata$session==1 & atdata$time<240])
bdat1$fit[2] <- mean(atdata$namp[atdata$burstlen==2 & atdata$burstpos==2 & atdata$session==1 & atdata$time<240])
bdat1$fit[3] <- mean(atdata$namp[atdata$burstlen==3 & atdata$burstpos==3 & atdata$session==1 & atdata$time<240])
bdat1$fit[4] <- mean(atdata$namp[atdata$burstlen==4 & atdata$burstpos==4 & atdata$session==1 & atdata$time<240])
bdat1$se[1] <- se(atdata$namp[atdata$burstlen==1 & atdata$burstpos==1 & atdata$session==1 & atdata$time<240])
bdat1$se[2] <- se(atdata$namp[atdata$burstlen==2 & atdata$burstpos==2 & atdata$session==1 & atdata$time<240])
bdat1$se[3] <- se(atdata$namp[atdata$burstlen==3 & atdata$burstpos==3 & atdata$session==1 & atdata$time<240])
bdat1$se[4] <- se(atdata$namp[atdata$burstlen==4 & atdata$burstpos==4 & atdata$session==1 & atdata$time<240])

bdat2 <- data.frame(burstpos = 1:4, fit = rep(1,times=4),se = rep(0,time=4))
bdat2$fit[1] <- mean(atdata$namp[atdata$burstlen==1 & atdata$burstpos==1 & atdata$session==1 & atdata$time>360])
bdat2$fit[2] <- mean(atdata$namp[atdata$burstlen==2 & atdata$burstpos==2 & atdata$session==1 & atdata$time>360])
bdat2$fit[3] <- mean(atdata$namp[atdata$burstlen==3 & atdata$burstpos==3 & atdata$session==1 & atdata$time>360])
bdat2$fit[4] <- mean(atdata$namp[atdata$burstlen==4 & atdata$burstpos==4 & atdata$session==1 & atdata$time>360])
bdat2$se[1] <- se(atdata$namp[atdata$burstlen==1 & atdata$burstpos==1 & atdata$session==1 & atdata$time>360])
bdat2$se[2] <- se(atdata$namp[atdata$burstlen==2 & atdata$burstpos==2 & atdata$session==1 & atdata$time>360])
bdat2$se[3] <- se(atdata$namp[atdata$burstlen==3 & atdata$burstpos==3 & atdata$session==1 & atdata$time>360])
bdat2$se[4] <- se(atdata$namp[atdata$burstlen==4 & atdata$burstpos==4 & atdata$session==1 & atdata$time>360])

bdat <- rbind(bdat1,bdat2)
bdat$condition <- factor(rep(c("early","late"),each=4))

ggplot(data = bdat,aes(x=burstpos,y=fit,group=condition,color=condition)) +
  geom_errorbar(aes(ymin=fit-2*se,ymax=fit+2*se),width = 0) +
  geom_line() + 
  geom_point(size=2) +
  scale_x_continuous(name="burst length (# of spikes)")+
  scale_y_continuous(name="last spike/first spike")+
  scale_color_manual(values=c("black","gray")) +
  theme(aspect.ratio = 1,
        panel.spacing.x=unit(0.2,"lines"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA,size=1))

ggsave("mean_atten.eps",device="eps",width=4,height=4,units="in")


ggplot(data = atdata[atdata$burstpos=="4",], aes(x=time,y=namp)) +
  geom_smooth(method="gam",formula=y~s(x,bs="cr",k=6,fx=F)) +
  facet_wrap(~session,nrow=1)

ggplot(data = atdata[atdata$cell==atcells[100]&atdata$session=="1",], aes(x=pos,y=zamp)) +
  geom_point() +
  geom_smooth(method="gam",formula=y~s(x,bs="cr"))

ggplot(data = atdata[atdata$session=="1"&atdata$burstpos!="1",],
       aes(x=time,y=dpamp)) +
  geom_hline(aes(yintercept=1),linetype="dashed")+
  geom_smooth(method="gam",formula=y~s(x,k=10,fx=F),color="black",alpha=1) +
  scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),
                     name="time (min)",limits=c(0,600),expand=c(0,0))+
  scale_y_continuous(name="third spike/first spike")+
  facet_wrap(~daymouse,nrow=6)+
  theme(aspect.ratio = 1,
        panel.spacing.x=unit(0.2,"lines"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA,size=1))

ggsave("third_namp_timexses.eps",device="eps",width=4,height=4,units="in")

atdata$burstpos<-factor(atdata$burstpos,ordered=FALSE)
mod<-bam(namp ~ session + s(time,by=session,fx=T),
                data = atdata[atdata$burstpos=="3",])

mod <- bam(burstlen ~ session + s(time,by=session,fx=TRUE),
           data = atdata,family=poisson(link=log))