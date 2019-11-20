library(mgcv)
library(ggplot2)
library(MASS)
library(itsadug)
library(RColorBrewer)
library(wesanderson)
library(directlabels)
library(scales)
library(akima)

setwd("E:/For Brian/GAMS_full/")

#load my data
cwd<-getwd()
MYDATA<-readRDS("MYDATAPOW_w.rds")
# MYDATA$Mouse<-ordered(MYDATA$Mouse)
# MYDATA<-MYDATA[MYDATA$Mouse>3,]

#Set model prediction vars
numdays<-length(unique(MYDATA$Day))
days<-seq(1,numdays)
nb<-50
nf<-50
freq<-10^seq(log10(2),log10(100),length.out=nf)
mint<-0
maxt<-600
tbins<-seq(mint,maxt,length.out=nb)
minrs<-0
maxrs<-0.9
rsbins<-seq(minrs,maxrs,length.out=nb)
pbins <- seq(-pi,pi,length.out=nb)

#create new data
dgrid<-expand.grid(pbins,rsbins,tbins,unique(MYDATA$Session),unique(MYDATA$Strain))
colnames(dgrid)<-c("phase","runningspeed","time","session","strain")
newdat<-with(MYDATA,data.frame(Theta_Phase = dgrid$phase,Running_Speed=dgrid$runningspeed,
                               Time=dgrid$time,Session=dgrid$session,Strain=dgrid$strain))

#define indices and remove unused predictions
ppidx<-which(!duplicated(newdat[,c('Theta_Phase','Strain')]))
pridx<-which(!duplicated(newdat[,c('Running_Speed','Strain')]))
ptidx<-which(!duplicated(newdat[,c('Time','Session','Strain')]))
used<-c(ppidx,pridx,ptidx)
newdat<-newdat[sort(unique(used)),]
newdat$Day<-MYDATA$Day[1]
newdat$Trial<-MYDATA$Trial[1]
newdat$Mouse<-MYDATA$Mouse[1]
newdat$StrainSession<-with(newdat,interaction(Strain,Session))

#define new indices
SessionStrainIdx <- which(!duplicated(newdat[,c('StrainSession')]))
pi1hidx<-SessionStrainIdx[1]
pi1aidx<-SessionStrainIdx[4]
pi2hidx<-SessionStrainIdx[2]
pi2aidx<-SessionStrainIdx[5]
pi3hidx<-SessionStrainIdx[3]
pi3aidx<-SessionStrainIdx[6]
ppidx<-which(!duplicated(newdat[,c('Theta_Phase','Strain')]))
pridx<-which(!duplicated(newdat[,c('Running_Speed','Strain')]))
ptidx<-which(!duplicated(newdat[,c('Time','StrainSession')]))
pt1idx<-ptidx[newdat$Session[ptidx]=="1"]
pt2idx<-ptidx[newdat$Session[ptidx]=="2"]
pt3idx<-ptidx[newdat$Session[ptidx]=="3"]

#preallocate plot data for iterms
lags<-30
acf<-matrix(0,nrow=nf,ncol=lags)
acfrho<-acf
rhofit<-matrix(0,nrow=nf,ncol=1)
tpfit<-matrix(0,nrow=nb*nf*2,ncol=3)
rsfit<-matrix(0,nrow=nb*nf*2,ncol=3)
timefit<-matrix(0,nrow=nb*nf*3*2,ncol=4)
Rsq<-matrix(0,nrow=1,ncol=nf)
I1h<-matrix(0,nrow=2,ncol=nf)

#set known values
tpfit[,2]<-rep(rep(freq,each=nb),times=2) #frequency
tpfit[,3]<-rep(pbins,times=nf*2) #runnings speed
rsfit[,2]<-rep(rep(freq,each=nb),times=2) #frequency
rsfit[,3]<-rep(rsbins,times=nf*2) #runnings speed
timefit[,2]<-rep(unique(MYDATA$Session),each=nb*nf*2) #session
timefit[,3]<-rep(rep(freq,each=nb),times=3*2) #frequency
timefit[,4]<-rep(tbins,times=3*nf*2) #time

#prediction using "lpmatrix"
ns<-3000
isim<-matrix(0,nrow = 3*nf*2,ncol = ns)
psim<-matrix(0,nrow = nf*nb*2,ncol = ns)
rsim<-matrix(0,nrow = nf*nb*2,ncol = ns)
timesim<-matrix(0,nrow = 3*nf*nb*2,ncol = ns)

fnum<-1:nf

for(f in 1:nf)
{

  fname<-paste0("warped_test_f",fnum[f],".rds")
  fmod<-readRDS(fname)

  Rsq[f] <- 1 - fmod$deviance/fmod$null.deviance

  #get mean coefficient values and variance-covariance matrix
  coefs<-coef(fmod)
  vc<-vcov(fmod)
  I1h[1,f] <- vc[1,1]
  I1h[2,f] <- coefs[1]
  rhofit[f] <- fmod$AR1.rho

  #draw random coefficient values from the model
  sim<-mvrnorm(ns,mu=coefs,Sigma=vc)

  #get predictions
  pred<-predict.bam(fmod, newdat, type = "lpmatrix", newdata.guaranteed = TRUE)

  rm(coefs,vc)
  gc()

  #define indices
  i1hidx<-f
  i1aidx<-3*nf+f
  i2hidx<-nf+f
  i2aidx<-4*nf+f
  i3hidx<-2*nf+f
  i3aidx<-5*nf+f
  pidx<-c(((f-1)*nb+1):(f*nb),((f-1)*nb+1+nf*nb):(f*nb+nf*nb))
  ridx<-c(((f-1)*nb+1):(f*nb),((f-1)*nb+1+nf*nb):(f*nb+nf*nb))
  t1idx<-c(((f-1)*nb+1):(f*nb),((f-1)*nb+1+3*nb*nf):(f*nb+3*nb*nf))
  t2idx<-c(((f-1)*nb+1+nb*nf):(f*nb+nb*nf),((f-1)*nb+1+4*nb*nf):(f*nb+4*nb*nf))
  t3idx<-c(((f-1)*nb+1+2*nb*nf):(f*nb+2*nb*nf),((f-1)*nb+1+5*nb*nf):(f*nb+5*nb*nf))
  
  #intercepts
  isim[i1hidx,] <- pred[pi1hidx, 1] %*% t(sim[, 1]) #session 1 hyb
  isim[i1aidx,] <- pred[pi1aidx, 2] %*% t(sim[, 2]) + pred[pi1hidx, 1] %*% t(sim[, 1]) #session 1 alz
  isim[i2hidx,] <- pred[pi2hidx, 3] %*% t(sim[, 3]) + pred[pi1hidx, 1] %*% t(sim[, 1]) #session 2 hyb
  isim[i2aidx,] <- pred[pi2aidx, 4] %*% t(sim[, 4]) + pred[pi1hidx, 1] %*% t(sim[, 1]) #session 2 alz
  isim[i3hidx,] <- pred[pi3hidx, 5] %*% t(sim[, 5]) + pred[pi1hidx, 1] %*% t(sim[, 1]) #session 3 hyb
  isim[i3aidx,] <- pred[pi3aidx, 6] %*% t(sim[, 6]) + pred[pi1hidx, 1] %*% t(sim[, 1]) #session 3 alz
  
  #thetaphase
  want<-grep("s\\(Theta\\_Phase\\)",names(fmod$coefficients))
  psim[pidx,]<-pred[ppidx,want] %*% t(sim[,want])

  #runningspeed
  want<-grep("s\\(Running\\_Speed\\)",names(fmod$coefficients))
  rsim[ridx,]<-pred[pridx,want] %*% t(sim[,want])

  #time by session
  want<-grep("s\\(Time\\)",names(fmod$coefficients))
  want1 <- want[1:(length(want)/3)]
  want2 <- want[(length(want)/3+1):(2*length(want)/3)]
  want3 <- want[(2*length(want)/3+1):(length(want))]
  timesim[t1idx,] <- pred[pt1idx, want1] %*% t(sim[, want1]) #session 1
  timesim[t2idx,] <- pred[pt2idx, want2] %*% t(sim[, want2]) #session 2
  timesim[t3idx,] <- pred[pt3idx, want3] %*% t(sim[, want3]) #session 3

  rm(pred,sim)
  gc()
  
  print(paste("Completed: ",f," of ",nf),quote=FALSE)
  
}
rm(newdat)
gc()

#define indices for plot terms
isimidx<-1:nrow(isim) #intercepts
psimidx<-(max(isimidx)+1):(max(isimidx)+nrow(psim)) #theta phase
rsimidx<-(max(psimidx)+1):(max(psimidx)+nrow(rsim)) #running speed
timesimidx<-(max(rsimidx)+1):(max(rsimidx)+nrow(timesim)) #time by session

#define indices for strain difference plot terms
disimidx<-1:(nrow(isim)/2) #intercepts
dpsimidx<-(max(disimidx)+1):(max(disimidx)+nrow(psim)/2) #theta phase
drsimidx<-(max(dpsimidx)+1):(max(dpsimidx)+nrow(rsim)/2) #running speed
dtimesimidx<-(max(drsimidx)+1):(max(drsimidx)+nrow(timesim)/2) #time by session

# #get simultaneous CI's for main terms
# sims<-rbind(isim,psim)
# rm(isim,psim)
# gc()
# sims<-rbind(sims,rsim)
# rm(rsim)
# gc()
# sims<-rbind(sims,timesim)
# rm(timesim)
# gc()
# 
# ci<-matrix(0,nrow = nrow(sims),ncol = 3) #matrix to hold results
# ci[,2]<-rowMeans(sims)
# sims_sd<-rep(0,times=nrow(sims))
# for(i in 1:nrow(sims)){
#   sims_sd[i]<-sd(sims[i,])
#   sims[i,]<-abs(sims[i,]-ci[i,2])/sims_sd[i]
# }
# gc()
# maxima<-rep(0,times=ns)
# for(i in 1:ns){maxima[i]<-max(sims[,i])}
# m <- quantile(maxima,probs=0.95,type=8)
# ci[,1]<--m*sims_sd
# ci[,3]<-m*sims_sd
# rm(sims)
# gc()


#get simultaneous CI's for strain difference terms
isim <- isim[(nrow(isim)/2+1):nrow(isim),] - isim[1:(nrow(isim)/2),]
psim <- psim[(nrow(psim)/2+1):nrow(psim),] - psim[1:(nrow(psim)/2),]
rsim <- rsim[(nrow(rsim)/2+1):nrow(rsim),] - rsim[1:(nrow(rsim)/2),]
timesim <- timesim[(nrow(timesim)/2+1):nrow(timesim),] - timesim[1:(nrow(timesim)/2),]

sims<-rbind(isim,psim)
rm(isim,psim)
gc()
sims<-rbind(sims,rsim)
rm(rsim)
gc()
sims<-rbind(sims,timesim)
rm(timesim)
gc()

dci<-matrix(0,nrow = nrow(sims),ncol = 3) #matrix to hold results
dci[,2]<-rowMeans(sims)
sims_sd<-rep(0,times=nrow(sims))
for(i in 1:nrow(sims)){
  sims_sd[i]<-sd(sims[i,])
  sims[i,]<-abs(sims[i,]-dci[i,2])/sims_sd[i]
}
gc()
maxima<-rep(0,times=ns)
for(i in 1:ns){maxima[i]<-max(sims[,i])}
m <- quantile(maxima,probs=0.95,type=8)
dci[,1]<--m*sims_sd
dci[,3]<-m*sims_sd
rm(sims)
gc()

################## Begin Plotting ##################################

#r-squared plot
rsqdf<-data.frame(Frequency=freq,
                  Rsquared=as.numeric(Rsq))
p<-ggplot(rsqdf,aes(x=Frequency,y=Rsquared)) 
p + geom_line()+geom_path(size=1)+
    scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",limits=c(2,100),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,0.5,by=0.05),name="R^2",limits=c(0,0.5),expand=c(0,0))+
    theme(axis.line.x = element_line(color="black", size = 1),
          axis.line.y = element_line(color="black", size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank())

#rho plot
rhodf<-data.frame(Frequency=rep(freq,times=2),
                  rho=c(acf[,1],acfrho[,1]),
                  model=rep(c("Normal","AR1"),each=nf))
p<-ggplot(rhodf,aes(x=Frequency,y=rho,group=model))
p+geom_line(size=1,aes(colour=model))+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,1,by=0.25),labels=seq(0,1,by=0.25),limits=c(-0.05,1),expand=c(0,0))+
  # geom_hline(aes(yintercept=0.05))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.ticks = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA))

#acf plot
acfdf<-as.data.frame(c(as.vector(acf),as.vector(acfrho)))
colnames(acfdf)<-c("acf")
acfdf$lags<-rep(rep(1:lags,each=nf),2)
acfdf$frequency<-rep(freq,times=2*lags)
acfdf$rho<-rep(c("normal","AR1"),each=nf*lags)

#interpolate
lbins<-1:30
fbins<-seq(2,100,length=250)
idxnorm<-acfdf$rho=="normal"
idxrho<-acfdf$rho=="AR1"
acfinterp<-interp(x=acfdf$lags[idxnorm],y=acfdf$frequency[idxnorm],acfdf$acf[idxnorm],xo=lbins,yo=fbins)
acfrhointerp<-interp(x=acfdf$lags[idxrho],y=acfdf$frequency[idxrho],acfdf$acf[idxrho],xo=lbins,yo=fbins)
acfnormdf <- data.frame(expand.grid(lags = acfinterp$x, frequency = acfinterp$y), acf = c(acfinterp$z))
acfrhodf <- data.frame(expand.grid(lags = acfrhointerp$x, frequency = acfrhointerp$y), acf = c(acfrhointerp$z))
acfnormdf$rho<-rep(factor("normal"),times=length(acfnormdf[,1]))
acfrhodf$rho<-rep(factor("AR1"),times=length(acfrhodf[,1]))
acfabsrhodf<-acfrhodf
acfabsrhodf$acf<-abs(acfnormdf$acf)-abs(acfrhodf$acf)
acfabsrhodf$rho<-"|Normal|-|AR1|"
acfdf<-rbind(acfnormdf,acfrhodf,acfabsrhodf)


#display plot
myPalette <- colorRampPalette(rev(brewer.pal(10, "PRGn")), space="Lab")
cl<-1
ggplot(acfdf,aes(x=lags,y=frequency,fill=acf,z=acf))+ 
  geom_raster()+scale_x_continuous(breaks=c(1,10,20,30),labels=c(1,10,20,30),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",alpha=0.3,binwidth=0.05,size=0.05,show.legend = TRUE)+
  facet_wrap(~rho,nrow=1)+theme(aspect.ratio = 1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA,size=1))

#intercepts by session
#run for "lpmatrix"
ifit<-as.data.frame(ci[isimidx,])
colnames(ifit)<-c("low","fit","high")
ifit$low<-ifit$low+ifit$fit
ifit$high<-ifit$high+ifit$fit
ifit$session<-rep(rep(c("1","2","3"),each=nf),times=2)
ifit$frequency<-rep(freq,times=6) 
ifit$strain<-rep(c("hyb","alz"),each=nf*3)
ifit$fit <- exp(ifit$fit)
ifit$low <- exp(ifit$low)
ifit$high <- exp(ifit$high)

p<-ggplot(ifit,aes(x=frequency,y=fit,fill=strain,group=strain))
p+geom_ribbon(aes(ymin=low,ymax=high),alpha=0.3)+
  geom_line(size=1,aes(color=strain))+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  facet_wrap(~session,nrow=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.ticks = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA))

ggsave("ampsession_Hyb.eps",device="eps",width=14,height=14,units="in")

#intercepts differences by session
#run for "lpmatrix"
ifit<-as.data.frame(dci[disimidx,])
colnames(ifit)<-c("low","fit","high")
ifit$low<-ifit$low+ifit$fit
ifit$high<-ifit$high+ifit$fit
ifit$session<-rep(c("1","2","3"),each=nf)
ifit$frequency<-rep(freq,times=3) 
ifit$fit <- exp(ifit$fit)
ifit$low <- exp(ifit$low)
ifit$high <- exp(ifit$high)

p<-ggplot(ifit,aes(x=frequency,y=fit))
p+geom_ribbon(aes(ymin=low,ymax=high),alpha=0.3)+
  geom_line(size=1)+
  geom_hline(yintercept = 1)+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  facet_wrap(~session,nrow=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.ticks = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA))

ggsave("ampsession_Hyb.eps",device="eps",width=14,height=14,units="in")

#theta phase plot
tpbins<-seq(-pi,pi,length=250)
fbins<-seq(min(freq),max(freq),length=250)
fithidx <- 1:(nrow(tpfit)/2)
fitaidx <- (nrow(tpfit)/2+1):nrow(tpfit)
cihidx <- psimidx[1:(length(psimidx)/2)]
ciaidx <- psimidx[(length(psimidx)/2+1):length(psimidx)]
fitinterp_h<-akima::interp(x=tpfit[fithidx,3],y=tpfit[fithidx,2],z=ci[cihidx,2],xo=tpbins,yo=fbins)
highinterp_h<-akima::interp(x=tpfit[fithidx,3],y=tpfit[fithidx,2],z=ci[cihidx,2]+ci[cihidx,3],xo=tpbins,yo=fbins)
lowinterp_h<-akima::interp(x=tpfit[fithidx,3],y=tpfit[fithidx,2],z=ci[cihidx,2]+ci[cihidx,1],xo=tpbins,yo=fbins)
fitinterp_a<-akima::interp(x=tpfit[fitaidx,3],y=tpfit[fitaidx,2],z=ci[ciaidx,2],xo=tpbins,yo=fbins)
highinterp_a<-akima::interp(x=tpfit[fitaidx,3],y=tpfit[fitaidx,2],z=ci[ciaidx,2]+ci[ciaidx,3],xo=tpbins,yo=fbins)
lowinterp_a<-akima::interp(x=tpfit[fitaidx,3],y=tpfit[fitaidx,2],z=ci[ciaidx,2]+ci[ciaidx,1],xo=tpbins,yo=fbins)
pfit <- data.frame(rbind(expand.grid(thetaphase = fitinterp_h$x, frequency = fitinterp_h$y), 
                         expand.grid(thetaphase = fitinterp_h$x, frequency = fitinterp_h$y)),
                   fit = c(fitinterp_h$z,fitinterp_a$z), 
                   high = c(highinterp_h$z,highinterp_a$z), 
                   low = c(lowinterp_h$z,lowinterp_a$z))
pfit$fit <- exp(pfit$fit)
pfit$low <- exp(pfit$low)
pfit$high <- exp(pfit$high)
pfit$strain <- rep(c("hyb","alz"),each=nrow(pfit)/2)

#display plot
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
cl<-0.5
ggplot(pfit,aes(x=thetaphase,y=frequency,fill=fit,z=fit))+
  geom_raster()+scale_x_continuous(breaks=seq(-pi,pi,by=pi/2),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(1-cl,1+cl),oob=squish)+
  stat_contour(mapping=aes(x=thetaphase,y=frequency,z=high),colour="blue",breaks=c(1))+
  stat_contour(mapping=aes(x=thetaphase,y=frequency,z=low),colour="red",breaks=c(1))+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  facet_wrap(~strain,nrow=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("speed_Alz.eps",device="eps",width=5,height=5,units="in")

#theta phase difference plot
tpfit <- tpfit[1:(nrow(tpfit)/2),]
tpbins<-seq(-pi,pi,length=250)
fbins<-seq(min(freq),max(freq),length=250)
fitinterp<-akima::interp(x=tpfit[,3],y=tpfit[,2],z=dci[dpsimidx,2],xo=tpbins,yo=fbins)
highinterp<-akima::interp(x=tpfit[,3],y=tpfit[,2],z=dci[dpsimidx,2]+dci[dpsimidx,3],xo=tpbins,yo=fbins)
lowinterp<-akima::interp(x=tpfit[,3],y=tpfit[,2],z=dci[dpsimidx,2]+dci[dpsimidx,1],xo=tpbins,yo=fbins)

pfit <- data.frame(expand.grid(thetaphase = fitinterp$x, frequency = fitinterp$y),
                   fit = c(fitinterp$z), 
                   high = c(highinterp$z), 
                   low = c(lowinterp$z))
pfit$fit <- exp(pfit$fit)
pfit$low <- exp(pfit$low)
pfit$high <- exp(pfit$high)

#display plot
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
cl<-0.5
ggplot(pfit,aes(x=thetaphase,y=frequency,fill=fit,z=fit))+
  geom_raster()+scale_x_continuous(breaks=seq(-pi,pi,by=pi/2),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(1-cl,1+cl),oob=squish)+
  stat_contour(mapping=aes(x=thetaphase,y=frequency,z=high),colour="blue",breaks=c(1))+
  stat_contour(mapping=aes(x=thetaphase,y=frequency,z=low),colour="red",breaks=c(1))+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

#running speed plot
rsfit[,3] <- (rsfit[,3]^2)/(2*pi)*pi*96.5
minrs <- min(rsfit[,3])
maxrs <- max(rsfit[,3])
rbins<-seq(minrs,maxrs,length=250)
fbins<-seq(min(freq),max(freq),length=250)
fithidx <- 1:(nrow(rsfit)/2)
fitaidx <- (nrow(rsfit)/2+1):nrow(rsfit)
cihidx <- rsimidx[1:(length(rsimidx)/2)]
ciaidx <- rsimidx[(length(rsimidx)/2+1):length(rsimidx)]
fitinterp_h<-akima::interp(x=rsfit[fithidx,3],y=rsfit[fithidx,2],z=ci[cihidx,2],xo=rbins,yo=fbins)
highinterp_h<-akima::interp(x=rsfit[fithidx,3],y=rsfit[fithidx,2],z=ci[cihidx,2]+ci[cihidx,3],xo=rbins,yo=fbins)
lowinterp_h<-akima::interp(x=rsfit[fithidx,3],y=rsfit[fithidx,2],z=ci[cihidx,2]+ci[cihidx,1],xo=rbins,yo=fbins)
fitinterp_a<-akima::interp(x=rsfit[fitaidx,3],y=rsfit[fitaidx,2],z=ci[ciaidx,2],xo=rbins,yo=fbins)
highinterp_a<-akima::interp(x=rsfit[fitaidx,3],y=rsfit[fitaidx,2],z=ci[ciaidx,2]+ci[ciaidx,3],xo=rbins,yo=fbins)
lowinterp_a<-akima::interp(x=rsfit[fitaidx,3],y=rsfit[fitaidx,2],z=ci[ciaidx,2]+ci[ciaidx,1],xo=rbins,yo=fbins)
rfit <- data.frame(rbind(expand.grid(runningspeed = fitinterp_h$x, frequency = fitinterp_h$y), 
                         expand.grid(runningspeed = fitinterp_h$x, frequency = fitinterp_h$y)),
                   fit = c(fitinterp_h$z,fitinterp_a$z), 
                   high = c(highinterp_h$z,highinterp_a$z), 
                   low = c(lowinterp_h$z,lowinterp_a$z))
rfit$fit <- exp(rfit$fit)
rfit$low <- exp(rfit$low)
rfit$high <- exp(rfit$high)
rfit$strain <- rep(c("hyb","alz"),each=nrow(rfit)/2)

#display plot
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
cl<-0.5
ggplot(rfit,aes(x=runningspeed,y=frequency,fill=fit,z=fit))+
  geom_raster()+scale_x_continuous(breaks=seq(round(minrs),maxrs,by=10),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(1-cl,1+cl),oob=squish)+
  stat_contour(mapping=aes(x=runningspeed,y=frequency,z=high),colour="blue",breaks=c(1))+
  stat_contour(mapping=aes(x=runningspeed,y=frequency,z=low),colour="red",breaks=c(1))+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  facet_wrap(~strain,nrow=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("speed_Alz.eps",device="eps",width=5,height=5,units="in")

#running speed difference plot
rsfit[,3] <- (rsfit[,3]^2)/(2*pi)*pi*96.5
rsfit <- rsfit[1:(nrow(rsfit)/2),]
minrs <- min(rsfit[,3])
maxrs <- max(rsfit[,3])
rbins<-seq(minrs,maxrs,length=250)
fbins<-seq(min(freq),max(freq),length=250)
fitinterp<-akima::interp(x=rsfit[,3],y=rsfit[,2],z=dci[drsimidx,2],xo=rbins,yo=fbins)
highinterp<-akima::interp(x=rsfit[,3],y=rsfit[,2],z=dci[drsimidx,2]+dci[drsimidx,3],xo=rbins,yo=fbins)
lowinterp<-akima::interp(x=rsfit[,3],y=rsfit[,2],z=dci[drsimidx,2]+dci[drsimidx,1],xo=rbins,yo=fbins)

rfit <- data.frame(expand.grid(runningspeed = fitinterp$x, frequency = fitinterp$y),
                   fit = c(fitinterp$z), 
                   high = c(highinterp$z), 
                   low = c(lowinterp$z))
rfit$fit <- exp(rfit$fit)
rfit$low <- exp(rfit$low)
rfit$high <- exp(rfit$high)

#display plot
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
cl<-0.5
ggplot(rfit,aes(x=runningspeed,y=frequency,fill=fit,z=fit))+
  geom_raster()+scale_x_continuous(breaks=seq(round(minrs),maxrs,by=10),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(1-cl,1+cl),oob=squish)+
  stat_contour(mapping=aes(x=runningspeed,y=frequency,z=high),colour="blue",breaks=c(1))+
  stat_contour(mapping=aes(x=runningspeed,y=frequency,z=low),colour="red",breaks=c(1))+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("speed_Alz.eps",device="eps",width=5,height=5,units="in")


#time by session
#run for "lpmatrix"
tfit<-as.data.frame(ci[timesimidx,])
colnames(tfit)<-c("low","fit","high")
tfit$low<-tfit$low+tfit$fit
tfit$high<-tfit$high+tfit$fit
tfit$session<-rep(rep(c("1","2","3"),each=nb*nf),times=2)
tfit$frequency<-rep(rep(freq,each=nb),times=6) 
tfit$time<-rep(seq(mint,maxt,length.out=nb),times=6*nf)
tfit$strain<-rep(c("hyb","alz"),each=nrow(tfit)/2)
tfit$SessionStrain <- with(tfit,interaction(session,strain))

#interpolate
tibins<-seq(mint,maxt,length=250)
fbins<-seq(min(freq),max(freq),length=250)
idxoneh<-tfit$SessionStrain == "1.hyb"
idxonea<-tfit$SessionStrain == "1.alz"
idxtwoh<-tfit$SessionStrain == "2.hyb"
idxtwoa<-tfit$SessionStrain == "2.alz"
idxthreeh<-tfit$SessionStrain == "3.hyb"
idxthreea<-tfit$SessionStrain == "3.alz"

fitsoneh<-akima::interp(x=tfit$time[idxoneh],y=tfit$frequency[idxoneh],tfit$fit[idxoneh],xo=tibins,yo=fbins)
fitstwoh<-akima::interp(x=tfit$time[idxtwoh],y=tfit$frequency[idxtwoh],tfit$fit[idxtwoh],xo=tibins,yo=fbins)
fitsthreeh<-akima::interp(x=tfit$time[idxthreeh],y=tfit$frequency[idxthreeh],tfit$fit[idxthreeh],xo=tibins,yo=fbins)
highoneh<-akima::interp(x=tfit$time[idxoneh],y=tfit$frequency[idxoneh],tfit$high[idxoneh],xo=tibins,yo=fbins)
hightwoh<-akima::interp(x=tfit$time[idxtwoh],y=tfit$frequency[idxtwoh],tfit$high[idxtwoh],xo=tibins,yo=fbins)
highthreeh<-akima::interp(x=tfit$time[idxthreeh],y=tfit$frequency[idxthreeh],tfit$high[idxthreeh],xo=tibins,yo=fbins)
lowoneh<-akima::interp(x=tfit$time[idxoneh],y=tfit$frequency[idxoneh],tfit$low[idxoneh],xo=tibins,yo=fbins)
lowtwoh<-akima::interp(x=tfit$time[idxtwoh],y=tfit$frequency[idxtwoh],tfit$low[idxtwoh],xo=tibins,yo=fbins)
lowthreeh<-akima::interp(x=tfit$time[idxthreeh],y=tfit$frequency[idxthreeh],tfit$low[idxthreeh],xo=tibins,yo=fbins)

fitsonea<-akima::interp(x=tfit$time[idxonea],y=tfit$frequency[idxonea],tfit$fit[idxonea],xo=tibins,yo=fbins)
fitstwoa<-akima::interp(x=tfit$time[idxtwoa],y=tfit$frequency[idxtwoa],tfit$fit[idxtwoa],xo=tibins,yo=fbins)
fitsthreea<-akima::interp(x=tfit$time[idxthreea],y=tfit$frequency[idxthreea],tfit$fit[idxthreea],xo=tibins,yo=fbins)
highonea<-akima::interp(x=tfit$time[idxonea],y=tfit$frequency[idxonea],tfit$high[idxonea],xo=tibins,yo=fbins)
hightwoa<-akima::interp(x=tfit$time[idxtwoa],y=tfit$frequency[idxtwoa],tfit$high[idxtwoa],xo=tibins,yo=fbins)
highthreea<-akima::interp(x=tfit$time[idxthreea],y=tfit$frequency[idxthreea],tfit$high[idxthreea],xo=tibins,yo=fbins)
lowonea<-akima::interp(x=tfit$time[idxonea],y=tfit$frequency[idxonea],tfit$low[idxonea],xo=tibins,yo=fbins)
lowtwoa<-akima::interp(x=tfit$time[idxtwoa],y=tfit$frequency[idxtwoa],tfit$low[idxtwoa],xo=tibins,yo=fbins)
lowthreea<-akima::interp(x=tfit$time[idxthreea],y=tfit$frequency[idxthreea],tfit$low[idxthreea],xo=tibins,yo=fbins)

tfitoneh <- data.frame(expand.grid(time = fitsoneh$x, frequency = fitsoneh$y), 
                      fit = c(fitsoneh$z), high = c(highoneh$z), low = c(lowoneh$z))
tfittwoh <- data.frame(expand.grid(time = fitstwoh$x, frequency = fitstwoh$y), 
                      fit = c(fitstwoh$z), high = c(hightwoh$z), low = c(lowtwoh$z))
tfitthreeh <- data.frame(expand.grid(time = fitsthreeh$x, frequency = fitsthreeh$y), 
                        fit = c(fitsthreeh$z), high = c(highthreeh$z), low = c(lowthreeh$z))

tfitonea <- data.frame(expand.grid(time = fitsonea$x, frequency = fitsonea$y), 
                       fit = c(fitsonea$z), high = c(highonea$z), low = c(lowonea$z))
tfittwoa <- data.frame(expand.grid(time = fitstwoa$x, frequency = fitstwoa$y), 
                       fit = c(fitstwoa$z), high = c(hightwoa$z), low = c(lowtwoa$z))
tfitthreea <- data.frame(expand.grid(time = fitsthreea$x, frequency = fitsthreea$y), 
                         fit = c(fitsthreea$z), high = c(highthreea$z), low = c(lowthreea$z))

tfitoneh$session<-rep(factor("1"),times=length(tfitoneh[,1]))
tfittwoh$session<-rep(factor("2"),times=length(tfittwoh[,1]))
tfitthreeh$session<-rep(factor("3"),times=length(tfitthreeh[,1]))

tfitonea$session<-rep(factor("1"),times=length(tfitonea[,1]))
tfittwoa$session<-rep(factor("2"),times=length(tfittwoa[,1]))
tfitthreea$session<-rep(factor("3"),times=length(tfitthreea[,1]))

tfit<-rbind(tfitoneh,tfittwoh,tfitthreeh,tfitonea,tfittwoa,tfitthreea)
tfit$fit <- exp(tfit$fit)
tfit$low <- exp(tfit$low)
tfit$high <- exp(tfit$high)
tfit$strain<-rep(c("hyb","alz"),each=nrow(tfit)/2)
tfit$SessionStrain <- with(tfit,interaction(session,strain))

#display plot
cl<-0.5
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
ggplot(tfit,aes(x=time,y=frequency,fill=fit,z=fit))+ 
  geom_raster()+scale_x_continuous(breaks=seq(0,600,by=60),labels=seq(0,10,by=1),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(1-cl,1+cl),oob=squish)+
  #stat_contour(colour="black",alpha=0.3,binwidth=0.05,size=0.05,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=frequency,z=high),colour="blue",breaks=c(1),alpha=1)+
  stat_contour(mapping=aes(x=time,y=frequency,z=low),colour="red",breaks=c(1),alpha=1)+
  facet_wrap(~SessionStrain,nrow=2)+theme(aspect.ratio = 1)+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("timesession_Hyb.eps",device="eps",width=14,height=14,units="in")


#strain difference time by session
tfit<-as.data.frame(dci[dtimesimidx,])
colnames(tfit)<-c("low","fit","high")
tfit$low<-tfit$low+tfit$fit
tfit$high<-tfit$high+tfit$fit
tfit$session<-rep(c("1","2","3"),each=nb*nf)
tfit$frequency<-rep(rep(freq,each=nb),times=3) 
tfit$time<-rep(seq(mint,maxt,length.out=nb),times=3*nf)

#interpolate
tibins<-seq(mint,maxt,length=250)
fbins<-seq(min(freq),max(freq),length=250)
idxone<-tfit$session == "1"
idxtwo<-tfit$session == "2"
idxthree<-tfit$session == "3"

fitsone<-akima::interp(x=tfit$time[idxone],y=tfit$frequency[idxone],tfit$fit[idxone],xo=tibins,yo=fbins)
fitstwo<-akima::interp(x=tfit$time[idxtwo],y=tfit$frequency[idxtwo],tfit$fit[idxtwo],xo=tibins,yo=fbins)
fitsthree<-akima::interp(x=tfit$time[idxthree],y=tfit$frequency[idxthree],tfit$fit[idxthree],xo=tibins,yo=fbins)
highone<-akima::interp(x=tfit$time[idxone],y=tfit$frequency[idxone],tfit$high[idxone],xo=tibins,yo=fbins)
hightwo<-akima::interp(x=tfit$time[idxtwo],y=tfit$frequency[idxtwo],tfit$high[idxtwo],xo=tibins,yo=fbins)
highthree<-akima::interp(x=tfit$time[idxthree],y=tfit$frequency[idxthree],tfit$high[idxthree],xo=tibins,yo=fbins)
lowone<-akima::interp(x=tfit$time[idxone],y=tfit$frequency[idxone],tfit$low[idxone],xo=tibins,yo=fbins)
lowtwo<-akima::interp(x=tfit$time[idxtwo],y=tfit$frequency[idxtwo],tfit$low[idxtwo],xo=tibins,yo=fbins)
lowthree<-akima::interp(x=tfit$time[idxthree],y=tfit$frequency[idxthree],tfit$low[idxthree],xo=tibins,yo=fbins)

tfitone <- data.frame(expand.grid(time = fitsone$x, frequency = fitsone$y), 
                       fit = c(fitsone$z), high = c(highone$z), low = c(lowone$z))
tfittwo <- data.frame(expand.grid(time = fitstwo$x, frequency = fitstwo$y), 
                       fit = c(fitstwo$z), high = c(hightwo$z), low = c(lowtwo$z))
tfitthree <- data.frame(expand.grid(time = fitsthree$x, frequency = fitsthree$y), 
                         fit = c(fitsthree$z), high = c(highthree$z), low = c(lowthree$z))

tfitone$session<-rep(factor("1"),times=length(tfitone[,1]))
tfittwo$session<-rep(factor("2"),times=length(tfittwo[,1]))
tfitthree$session<-rep(factor("3"),times=length(tfitthree[,1]))

tfit<-rbind(tfitone,tfittwo,tfitthree)
tfit$fit <- exp(tfit$fit)
tfit$low <- exp(tfit$low)
tfit$high <- exp(tfit$high)

#display plot
cl<-0.25
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
ggplot(tfit,aes(x=time,y=frequency,fill=fit,z=fit))+ 
  geom_raster()+scale_x_continuous(breaks=seq(0,600,by=60),labels=seq(0,10,by=1),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(1-cl,1+cl),oob=squish)+
  #stat_contour(colour="black",alpha=0.3,binwidth=0.05,size=0.05,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=frequency,z=high),colour="blue",breaks=c(1),alpha=1)+
  stat_contour(mapping=aes(x=time,y=frequency,z=low),colour="red",breaks=c(1),alpha=1)+
  facet_wrap(~session,nrow=1)+theme(aspect.ratio = 1)+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("timesession_Hyb.eps",device="eps",width=14,height=14,units="in")

