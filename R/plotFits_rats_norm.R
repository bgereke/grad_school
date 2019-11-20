library(mgcv)
library(ggplot2)
library(MASS)
library(itsadug)
library(RColorBrewer)
library(wesanderson)
library(directlabels)
library(scales)
library(akima)

setwd("G:/Rats/GAMS_full_CCS")

#load my data
cwd<-getwd()
MYDATA<-readRDS("MYDATA_full_CCS.rds")

#Set model prediction vars
numdays<-length(unique(MYDATA$Day))
days<-seq(1,numdays)
nb<-45
nf<-50
freq<-10^seq(log10(2),log10(100),length.out=nf)
mint<-0
maxt<-600
tbins<-seq(mint,maxt,length.out=nb)
minrs<-min(MYDATA$Running_Speed)
maxrs<-quantile(MYDATA$Running_Speed,0.999) #max(MYDATA$Running_Speed)
rsbins<-seq(minrs,maxrs,length.out=nb)

#create new data
dgrid<-expand.grid(rsbins,tbins,unique(MYDATA$Session))
colnames(dgrid)<-c("runningspeed","time","session")
newdat<-with(MYDATA,data.frame(Running_Speed=dgrid$runningspeed,Time=dgrid$time,
                               Session=dgrid$session))

#define indices and remove unused predictions
pi1idx<-which(!duplicated(newdat[,c('Session')]))[1]
pi2idx<-which(!duplicated(newdat[,c('Session')]))[2]
pi3idx<-which(!duplicated(newdat[,c('Session')]))[3]
pridx<-which(!duplicated(newdat[,c('Running_Speed')]))
ptidx<-which(!duplicated(newdat[,c('Time','Session')]))
prtintidx<-which(!duplicated(newdat[,c('Running_Speed','Time','Session')]))
used<-c(pi1idx,pi2idx,pi3idx,pridx,ptidx,prtintidx)
# used<-c(pridx,ptidx)
newdat<-newdat[sort(unique(used)),]
newdat$Day<-MYDATA$Day[1]
newdat$Trial<-MYDATA$Trial[1]
newdat$Mouse<-MYDATA$Mouse[1]

#define new indices
pi1idx<-which(!duplicated(newdat[,c('Session')]))[1]
pi2idx<-which(!duplicated(newdat[,c('Session')]))[2]
pi3idx<-which(!duplicated(newdat[,c('Session')]))[3]
pridx<-which(!duplicated(newdat[,c('Running_Speed')]))
ptidx<-which(!duplicated(newdat[,c('Time','Session')]))
prtintidx<-which(!duplicated(newdat[,c('Running_Speed','Time','Session')]))
pt1idx<-ptidx[newdat$Session[ptidx]=="1"]
pt2idx<-ptidx[newdat$Session[ptidx]=="2"]
pt3idx<-ptidx[newdat$Session[ptidx]=="3"]
prtint1idx<-prtintidx[newdat$Session[prtintidx]=="1"]
prtint2idx<-prtintidx[newdat$Session[prtintidx]=="2"]
prtint3idx<-prtintidx[newdat$Session[prtintidx]=="3"]


#preallocate plot data for iterms
lags<-30
acf<-matrix(0,nrow=nf,ncol=lags)
acfrho<-acf
rhofit<-matrix(0,nrow=nf,ncol=1)
# ifit<-matrix(0,nrow=nf,ncol=3)
rsfit<-matrix(0,nrow=nb*nf,ncol=3)
timefit<-matrix(0,nrow=nb*nf*3,ncol=4)
rtintfit<-matrix(0,nrow=nb*nb*nf*3,ncol=5)
fRsq<-matrix(0,nrow=1,ncol=nf)
# rsse<-matrix(0,nrow=nb*nf,ncol=1)
# timese<-matrix(0,nrow=nb*nf*3,ncol=1)
# intse<-matrix(0,nrow=nb*nb*nf*3,ncol=1)

#set known values
rsfit[,2]<-rep(freq,each=nb) #frequency
rsfit[,3]<-rep(rsbins,times=nf) #runnings speed
timefit[,2]<-rep(unique(MYDATA$Session),each=nb*nf) #session
timefit[,3]<-rep(rep(freq,each=nb),times=3) #frequency
timefit[,4]<-rep(tbins,times=3*nf) #time
rtintfit[,2]<-rep(unique(MYDATA$Session),each=nb*nb*nf) #session
rtintfit[,3]<-rep(rep(freq,each=nb*nb),times=3) #frequency
rtintfit[,4]<-rep(rep(tbins,each=nb),times=nf*3) #time
rtintfit[,5]<-rep(rsbins,times=nb*nf*3) #running speed

# setwd("C:/Data/GAMS_full_rho")

#prediction using "lpmatrix"
ns<-1500
isim<-matrix(0,nrow = 3*nf,ncol = ns)
rsim<-matrix(0,nrow = nf*nb,ncol = ns)
timesim<-matrix(0,nrow = 3*nf*nb,ncol = ns)
rtintsim<-matrix(0,nrow = 3*nf*nb*nb,ncol = ns)
fnum<-1:nf

# r<-matrix(0,nrow = nf,ncol = 1)
# for(f in 1:nf){r[f]<-ar(MYDATA[,f],order.max=1)$ar}

for(f in 1:nf)
{
  # setwd("C:/Data/GAMS_full")
  filename<-paste0("CCS_f",fnum[f],".rds")
  fmod<-readRDS(filename)
  # setwd("C:/Data/GAMS_full_rho")
  # filename<-paste0("Base_RT_TP_f",fnum[f],"_full_rho_PP.rds")
  # fmodrho<-readRDS(filename)
  
  #rho's
  # res<-resid_gam(fmod,incl_na = TRUE,return_all = TRUE)
  # acf[f,]<-acf_resid(fmod,split_pred = list(MYDATA$Trial),plot=FALSE,max_lag=lags+1)[-1]
  # acfrho[f,]<-acf_resid(fmodrho,split_pred = c("Trial"),plot=FALSE,max_lag=lags+1)[-1]
  # rhofit[f]<-fmodrho$AR1.rho
  fRsq[f] <- 1 - fmod$deviance/fmod$null.deviance
  #
  #get mean coefficient values and variance-covariance matrix
  coefs<-coef(fmod)
  vc<-vcov(fmod)
  #draw random coefficient values from the model
  # set.seed(f)

  sim<-mvrnorm(ns,mu=coefs,Sigma=vc)

  #get predictions
  pred<-predict.bam(fmod, newdat, type = "lpmatrix", newdata.guaranteed = TRUE)

  rm(coefs,vc)
  gc()

  #define indices
  ridx<-((f-1)*nb+1):(f*nb)
  t1idx<-((f-1)*nb+1):(f*nb)
  t2idx<-((f-1)*nb+1+nb*nf):(f*nb+nb*nf)
  t3idx<-((f-1)*nb+1+2*nb*nf):(f*nb+2*nb*nf)
  # rtint1idx<-((f-1)*nb*nb+1):(f*nb*nb)
  # rtint2idx<-((f-1)*nb*nb+1+nb*nb*nf):(f*nb*nb+nf*nb*nb)
  # rtint3idx<-((f-1)*nb*nb+1+2*nb*nb*nf):(f*nb*nb+2*nf*nb*nb)

  #intercepts
  isim[f,] <- pred[pi1idx, 1] %*% t(sim[, 1]) #session 1
  isim[nf+f,] <- pred[pi2idx, 2] %*% t(sim[, 2]) + isim[f,] #session 2
  isim[2*nf+f,] <- pred[pi3idx, 3] %*% t(sim[, 3]) + isim[f,] #session 3

  #runningspeed
  want<-grep("s\\(Running\\_Speed\\)",names(fmod$coefficients))
  rsim[ridx,]<-pred[pridx,want] %*% t(sim[,want])


  #time by session
  want1<-grep("s\\(Time\\)\\:Session1",names(fmod$coefficients))
  want2<-grep("s\\(Time\\)\\:Session2",names(fmod$coefficients))
  want3<-grep("s\\(Time\\)\\:Session3",names(fmod$coefficients))
  timesim[t1idx,] <- pred[pt1idx, want1] %*% t(sim[, want1]) #session 1
  timesim[t2idx,] <- pred[pt2idx, want2] %*% t(sim[, want2]) #session 2
  timesim[t3idx,] <- pred[pt3idx, want3] %*% t(sim[, want3]) #session 3

  # #time/running speed interactions
  # want1<-grep("ti\\(Time\\,Running\\_Speed\\)\\:Session1",names(fmod$coefficients))
  # want2<-grep("ti\\(Time\\,Running\\_Speed\\)\\:Session2",names(fmod$coefficients))
  # want3<-grep("ti\\(Time\\,Running\\_Speed\\)\\:Session3",names(fmod$coefficients))
  # rtintsim[rtint1idx,]<-pred[prtint1idx,want1] %*% t(sim[,want1]) #session 1
  # rtintsim[rtint2idx,]<-pred[prtint2idx,want2] %*% t(sim[,want2]) #session 2
  # rtintsim[rtint3idx,]<-pred[prtint3idx,want3] %*% t(sim[,want3]) #session 3

  rm(pred,sim)
  gc()
  
  print(paste("Completed: ",f," of ",nf),quote=FALSE)
  
}
rm(newdat)
gc()

#define indices for plot terms
isimidx<-1:nrow(isim) #intercepts
rsimidx<-(max(isimidx)+1):(max(isimidx)+nrow(rsim)) #running speed
# rsimidx<-1:nrow(rsim) #running speed
timesimidx<-(max(rsimidx)+1):(max(rsimidx)+nrow(timesim)) #time by session
# rtintsimidx<-(max(timesimidx)+1):(max(timesimidx)+nrow(rtintsim)) #time/running speed interaction by session

#get simultaneous CI's for main terms
sims<-rbind(isim,rsim)
rm(isim,rsim)
# sims<-rbind(rsim)
# rm(rsim)
gc()
sims<-rbind(sims,timesim)
rm(timesim)
gc()
# sims<-rbind(sims,rtintsim)
# rm(rtintsim)
# gc()

# source("C:/Users/Brian/Documents/R/simCIz.R")
# ci<-simCIz(sims)
ci<-matrix(0,nrow = nrow(sims),ncol = 3) #matrix to hold results
ci[,2]<-rowMeans(sims)
sims_sd<-rep(0,times=nrow(sims))
for(i in 1:nrow(sims)){
  sims_sd[i]<-sd(sims[i,])
  sims[i,]<-abs(sims[i,]-ci[i,2])/sims_sd[i]
}
gc()
maxima<-rep(0,times=ns)
for(i in 1:ns){maxima[i]<-max(sims[,i])}
m <- quantile(maxima,probs=0.95,type=8)
ci[,1]<--m*sims_sd
ci[,3]<-m*sims_sd

rm(sims,maxima)
gc()

# #get difference terms
# disim<-matrix(0,nrow = nrow(isim),ncol = ncol(isim))
# disim[1:nf,]<-isim[1:nf,]-isim[(1+nf):(2*nf),]
# disim[(1+nf):(2*nf),]<-isim[1:nf,]-isim[(1+2*nf):(3*nf),]
# disim[(1+2*nf):(3*nf),]<-isim[(1+nf):(2*nf),]-isim[(1+2*nf):(3*nf),]
# rm(isim)
# gc()
# dtimesim<-matrix(0,nrow = nrow(timesim),ncol = ncol(timesim))
# dtimesim[1:(nf*nb),]<-timesim[1:(nf*nb),]-timesim[(1+nb*nf):(2*nb*nf),]
# dtimesim[(1+nb*nf):(2*nb*nf),]<-timesim[1:(nf*nb),]-timesim[(1+2*nb*nf):(3*nb*nf),]
# dtimesim[(1+2*nb*nf):(3*nb*nf),]<-timesim[(1+nb*nf):(2*nb*nf),]-timesim[(1+2*nb*nf):(3*nb*nf),]
# rm(timesim)
# gc()
# # drtintsim<-matrix(0,nrow = nrow(rtintsim),ncol(rtintsim))
# # drtintsim[1:(nf*nb*nb),]<-rtintsim[1:(nf*nb*nb),]-rtintsim[(1+nb*nb*nf):(2*nb*nb*nf),]
# # drtintsim[(1+nb*nb*nf):(2*nb*nb*nf),]<-rtintsim[1:(nf*nb*nb),]-rtintsim[(1+2*nb*nb*nf):(3*nb*nb*nf),]
# # drtintsim[(1+2*nb*nb*nf):(3*nb*nb*nf),]<-rtintsim[(1+nb*nb*nf):(2*nb*nb*nf),]-rtintsim[(1+2*nb*nb*nf):(3*nb*nb*nf),]
# # rm(rtintsim)
# # gc()
# 
# #get simultaneous CI's for difference terms
# sims<-rbind(disim,rsim)
# rm(disim,rsim)
# gc()
# # sims<-rbind(rsim)
# # rm(rsim)
# # gc()
# sims<-rbind(sims,dtimesim)
# rm(dtimesim)
# gc()
# # sims<-rbind(sims,drtintsim)
# # rm(drtintsim)
# # gc()
# 
# # source("C:/Users/Brian/Documents/R/simCIz.R")
# # ci<-simCIz(sims)
# dci<-matrix(0,nrow = nrow(sims),ncol = 3) #matrix to hold results
# dci[,2]<-rowMeans(sims)
# sims_sd<-rep(0,times=nrow(sims))
# for(i in 1:nrow(sims)){
#   sims_sd[i]<-sd(sims[i,])
#   sims[i,]<-abs(sims[i,]-dci[i,2])/sims_sd[i]
# }
# gc()
# maxima<-rep(0,times=ns)
# for(i in 1:ns){maxima[i]<-max(sims[,i])}
# m <- quantile(maxima,probs=0.95,type=8)
# dci[,1]<--m*sims_sd
# dci[,3]<-m*sims_sd
# 
# rm(sims,maxima)
# gc()

#r-squared plot
rsqdf<-data.frame(Frequency=freq,Rsquared=as.numeric(fRsq))
p<-ggplot(rsqdf,aes(x=Frequency,y=Rsquared))
p+geom_line()+geom_path(size=1)+
  geom_hline(aes(yintercept=0),linetype="dashed")+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,0.5,by=0.025),name="R^2")+
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
ifit$session<-rep(rep(c("1","2","3")),each=nf)
ifit$frequency<-rep(freq,times=3) 
# ifit$fit <- exp(ifit$fit)
# ifit$low <- exp(ifit$low)
# ifit$high <- exp(ifit$high)

p<-ggplot(ifit,aes(x=frequency,y=fit))
p+geom_ribbon(aes(ymin=low,ymax=high),color="grey95")+
  geom_line(size=1)+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  geom_hline(aes(yintercept=0),linetype="dashed")+
  facet_wrap(~session,nrow=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.ticks = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA))

ggsave("ampsession_CA1.eps",device="eps",width=14,height=14,units="in")

# n1 <- sum(MYDATA$Session=="1")
# n2 <- sum(MYDATA$Session=="2")
# n3 <- sum(MYDATA$Session=="3")
# fits <- (n1*ifit$fit[ifit$session=="1"]+n2*ifit$fit[ifit$session=="2"]+n3*ifit$fit[ifit$session=="3"])/(n1+n2+n3)

#intercept differences by session
#run for "lpmatrix"
ifit<-as.data.frame(dci[isimidx,])
colnames(ifit)<-c("low","fit","high")
ifit$low<-ifit$low+ifit$fit
ifit$high<-ifit$high+ifit$fit
ifit$session<-rep(rep(c("1-2","1-3","2-3")),each=nf)
ifit$frequency<-rep(freq,times=3) 
# ifit$fit <- exp(ifit$fit)
# ifit$low <- exp(ifit$low)
# ifit$high <- exp(ifit$high)

p<-ggplot(ifit,aes(x=frequency,y=fit))
p+geom_ribbon(aes(ymin=low,ymax=high))+
  geom_line(size=1)+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  geom_hline(aes(yintercept=0),linetype="dashed")+
  facet_wrap(~session,nrow=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.ticks = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA))

ggsave("dampsession_CA1.eps",device="eps",width=14,height=14,units="in")

#running speed plot
#run for "lpmatrix"
rbins<-seq(rsbins[1],rsbins[length(rsbins)],length=250)
fbins<-seq(min(freq),max(freq),length=250)
fitinterp<-akima::interp(x=rsfit[,3],y=rsfit[,2],ci[rsimidx,2],xo=rbins,yo=fbins)
highinterp<-akima::interp(x=rsfit[,3],y=rsfit[,2],ci[rsimidx,2]+ci[rsimidx,3],xo=rbins,yo=fbins)
lowinterp<-akima::interp(x=rsfit[,3],y=rsfit[,2],ci[rsimidx,2]+ci[rsimidx,1],xo=rbins,yo=fbins)
rfit <- data.frame(expand.grid(runningspeed = fitinterp$x, frequency = fitinterp$y), 
                   fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z))
# rfit$fit <- exp(rfit$fit)
# rfit$low <- exp(rfit$low)
# rfit$high <- exp(rfit$high)

#display plot
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
cl<-0.5
ggplot(rfit,aes(x=runningspeed,y=frequency,fill=fit,z=fit))+
  geom_raster()+scale_x_continuous(breaks=seq(round(minrs),maxrs,by=10),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-0.25,0.25),oob=squish)+
  stat_contour(mapping=aes(x=runningspeed,y=frequency,z=high),colour="blue",breaks=c(0))+
  stat_contour(mapping=aes(x=runningspeed,y=frequency,z=low),colour="red",breaks=c(0))+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("speed_CA1.eps",device="eps",width=5,height=5,units="in")

#theta phase plot
#run for "lpmatrix"
ipbins<-seq(pbins[1],pbins[length(pbins)],length=250)
fbins<-seq(min(freq),max(freq),length=250)
fitinterp<-interp(x=phasefit[,3],y=phasefit[,2],ci[psimidx,2],xo=ipbins,yo=fbins)
highinterp<-interp(x=phasefit[,3],y=phasefit[,2],ci[psimidx,2]+ci[psimidx,3],xo=ipbins,yo=fbins)
lowinterp<-interp(x=phasefit[,3],y=phasefit[,2],ci[psimidx,2]+ci[psimidx,1],xo=ipbins,yo=fbins)
pfit <- data.frame(expand.grid(phase = fitinterp$x, frequency = fitinterp$y), 
                   fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z))

# pfit$phase<-pfit$phase+pi
# pfit$phase[pfit$phase>pi]<-pfit$phase[pfit$phase>pi]-2*pi-median(diff(ipbins))

#display plot
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
cl<-max(abs(pfit$fit))
ggplot(pfit,aes(x=phase,y=frequency,fill=fit,z=fit))+
  geom_raster()+scale_x_continuous(breaks=seq(-pi,pi,by=pi/2),
                                   labels=c(expression(paste("-",pi)),expression(paste("-",pi,"/2")),
                                            "0",expression(paste(pi,"/2")),expression(paste(pi))),
                                   expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  stat_contour(mapping=aes(x=phase,y=frequency,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=phase,y=frequency,z=low),colour="red",breaks=c(0),alpha=1)+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("phase_Alz.eps",device="eps",width=5,height=5,units="in")


#time by session
#run for "lpmatrix"
tfit<-as.data.frame(ci[timesimidx,])
colnames(tfit)<-c("low","fit","high")
tfit$low<-tfit$low+tfit$fit
tfit$high<-tfit$high+tfit$fit
tfit$session<-rep(rep(c("1","2","3")),each=nb*nf)
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
# tfit$fit <- exp(tfit$fit)
# tfit$low <- exp(tfit$low)
# tfit$high <- exp(tfit$high)

#display plot
cl<-0.5#max(abs(tfit$fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
ggplot(tfit,aes(x=time,y=frequency,fill=fit,z=fit))+ 
  geom_raster()+scale_x_continuous(breaks=seq(0,600,by=60),labels=seq(0,10,by=1),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-0.15,0.15),oob=squish)+
  #stat_contour(colour="black",alpha=0.3,binwidth=0.05,size=0.05,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=frequency,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=frequency,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~session,nrow=1)+theme(aspect.ratio = 1)+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("timesession_CA1.eps",device="eps",width=14,height=14,units="in")

#differences time by session

#run for "lpmatrix"
tfit<-as.data.frame(dci[timesimidx,])
colnames(tfit)<-c("low","fit","high")
tfit$low<-tfit$low+tfit$fit
tfit$high<-tfit$high+tfit$fit
tfit$session<-rep(rep(c("1","2","3")),each=nb*nf)
tfit$frequency<-rep(rep(freq,each=nb),times=3) 
tfit$time<-rep(seq(mint,maxt,length.out=nb),times=3*nf)

#interpolate for "lpmatrix"
tibins<-seq(tbins[1],tbins[length(tbins)],length=250)
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
# tfit$fit <- exp(tfit$fit)
# tfit$low <- exp(tfit$low)
# tfit$high <- exp(tfit$high)

#display plot
tfit$dsessions<-as.factor(c(rep("1-2",times=250*250),rep("1-3",times=250*250),rep("2-3",times=250*250)))
cl<-0.5#max(abs(tfit$fit))
p<-ggplot(tfit,aes(x=time,y=frequency,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster()+scale_x_continuous(breaks=seq(0,600,by=60),labels=seq(0,10,by=1),expand=c(0,0))+
  scale_y_continuous(breaks=seq(10,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",alpha=0.3,binwidth=0.05,size=0.05,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=frequency,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=frequency,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~dsessions,nrow=1)+theme(aspect.ratio = 1)+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dtimesession_CA1.eps",device="eps",width=14,height=14,units="in")  

#time/running speed interactions by frequency session 1
#run for lpmatrix
trsfit<-as.data.frame(ci[rtintsimidx,])
colnames(trsfit)<-c("low","fit","high")
trsfit$low<-trsfit$low+trsfit$fit
trsfit$high<-trsfit$high+trsfit$fit
trsfit$session<-rtintfit[,2]
trsfit$frequency<-rtintfit[,3]
trsfit$time<-rtintfit[,4]
trsfit$runningspeed<-rtintfit[,5]

#set NA's
th = 1;
b = 1;
MYDATA$Time<-MYDATA$Time-min(MYDATA$Time)
MYDATA$Time<-MYDATA$Time/max(MYDATA$Time)
trsfit$time<-trsfit$time-min(trsfit$time)
trsfit$time<-trsfit$time/max(trsfit$time)
require(gplots)
h1<-hist2d(MYDATA$Time[MYDATA$Session=="1"],MYDATA$Running_Speed[MYDATA$Session=="1"],nbins=nb/b,show=FALSE)
h2<-hist2d(MYDATA$Time[MYDATA$Session=="2"],MYDATA$Running_Speed[MYDATA$Session=="2"],nbins=nb/b,show=FALSE)
h3<-hist2d(MYDATA$Time[MYDATA$Session=="3"],MYDATA$Running_Speed[MYDATA$Session=="3"],nbins=nb/b,show=FALSE)
counts1<-matrix(0,nb,nb)
counts2<-matrix(0,nb,nb)
counts3<-matrix(0,nb,nb)
for(i in 1:30){
  counts1[2*i,]<-rep(h1$counts[i,],each=b)
  counts1[2*i-1,]<-rep(h1$counts[i,],each=b)
  counts2[2*i,]<-rep(h2$counts[i,],each=b)
  counts2[2*i-1,]<-rep(h2$counts[i,],each=b)
  counts3[2*i,]<-rep(h3$counts[i,],each=b)
  counts3[2*i-1,]<-rep(h3$counts[i,],each=b)
}
counts1<-rep(as.vector(t(counts1)),times=nf)
counts2<-rep(as.vector(t(counts2)),times=nf)
counts3<-rep(as.vector(t(counts3)),times=nf)
counts<-c(counts1,counts2,counts3)
trsfit$fit[counts<th]<-NA
trsfit$low[counts<th]<-NA
trsfit$high[counts<th]<-NA
trsfit$fit <-exp(trsfit$fit)
trsfit$low <- exp(trsfit$low)
trsfit$high <- exp(trsfit$high)

#display plot
trsfit$frequency<-round(trsfit$frequency,digits=1)
trsfit$frequency<-as.ordered(trsfit$frequency)
for(f in 1:nf){
  levels(trsfit$frequency)[f]<-paste0(levels(trsfit$frequency)[f]," Hz")
}
cl<-0.5#0.7*max(abs(trsfit$fit),na.rm=TRUE)
p<-ggplot(trsfit[trsfit$session=="1",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(1-cl,1+cl),oob=squish,na.value="white")+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="blue",breaks=c(1),alpha=1)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="red",breaks=c(1),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("speedtimesession1_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time/running speed interactions by frequency session 2
p<-ggplot(trsfit[trsfit$session=="2",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(1-cl,1+cl),oob=squish,na.value="white")+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="blue",breaks=c(1),alpha=1)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="red",breaks=c(1),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("speedtimesession2_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time/running speed interactions by frequency session 3
p<-ggplot(trsfit[trsfit$session=="3",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(1-cl,1+cl),oob=squish,na.value="white")+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="blue",breaks=c(1),alpha=1)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="red",breaks=c(1),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("speedtimesession3_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time/running speed interactions by frequency session 1-2
trsfit<-as.data.frame(dci[rtintsimidx,])
colnames(trsfit)<-c("low","fit","high")
trsfit$low<-trsfit$low+trsfit$fit
trsfit$high<-trsfit$high+trsfit$fit
trsfit$session<-rep(c("1","2","3"),each=nb*nb*nf)
trsfit$frequency<-rep(rep(freq,each=nb*nb),times=3)
trsfit$time<-rep(rep(seq(mint,maxt,length.out=nb),each=nb),times=nf*3)
trsfit$runningspeed<-rep(rsbins,times=nb*nf*3)

#set NA's
th = 1;
b = 1;
MYDATA$Time<-MYDATA$Time-min(MYDATA$Time)
MYDATA$Time<-MYDATA$Time/max(MYDATA$Time)
trsfit$time<-trsfit$time-min(trsfit$time)
trsfit$time<-trsfit$time/max(trsfit$time)
require(gplots)
h1<-hist2d(MYDATA$Time[MYDATA$Session=="1"],MYDATA$Running_Speed[MYDATA$Session=="1"],nbins=nb/b,show=FALSE)
h2<-hist2d(MYDATA$Time[MYDATA$Session=="2"],MYDATA$Running_Speed[MYDATA$Session=="2"],nbins=nb/b,show=FALSE)
h3<-hist2d(MYDATA$Time[MYDATA$Session=="3"],MYDATA$Running_Speed[MYDATA$Session=="3"],nbins=nb/b,show=FALSE)
counts1<-rep(as.vector(t(h1)),times=nf)
counts2<-rep(as.vector(t(h2)),times=nf)
counts3<-rep(as.vector(t(h3)),times=nf)
nas1<-rep(0,times=length(counts1))
nas2<-nas1
nas3<-nas1
nas1[counts1<th | counts2<th]<-NA
nas2[counts1<th | counts3<th]<-NA
nas3[counts2<th | counts3<th]<-NA
nas<-c(nas1,nas2,nas3)
trsfit$fit[is.na(nas)]<-NA
trsfit$low[is.na(nas)]<-NA
trsfit$high[is.na(nas)]<-NA

trsfit$frequency<-round(trsfit$frequency,digits=1)
trsfit$frequency<-as.ordered(trsfit$frequency)
for(f in 1:nf){
  levels(trsfit$frequency)[f]<-paste0(levels(trsfit$frequency)[f]," Hz")
}
cl<-0.8*max(abs(trsfit$fit),na.rm=TRUE)
p<-ggplot(trsfit[trsfit$session=="1",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish,na.value="white")+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="blue",breaks=c(0))+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="red",breaks=c(0))+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dspeedtimesession1_Hyb.eps",device="eps",width=5.5,height=11,units="in")

#time/running speed interactions by frequency session 1-3
p<-ggplot(trsfit[trsfit$session=="2",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish,na.value="white")+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="blue",breaks=c(0))+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="red",breaks=c(0))+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dspeedtimesession2_Hyb.eps",device="eps",width=5.5,height=11,units="in")

#time/running speed interactions by frequency session 2-3
p<-ggplot(trsfit[trsfit$session=="3",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish,na.value="white")+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="blue",breaks=c(0))+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="red",breaks=c(0))+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dspeedtimesession3_Hyb.eps",device="eps",width=5.5,height=11,units="in")

#running speed theta phase interaction tiled by frequency
#run for lpmatrix
rspfit<-as.data.frame(ci[rpintsimidx,])
colnames(rspfit)<-c("low","fit","high")
rspfit$low<-rspfit$low+rspfit$fit
rspfit$high<-rspfit$high+rspfit$fit
rspfit$frequency<-rpintfit[,2]
rspfit$phase<-rpintfit[,3]
rspfit$runningspeed<-rpintfit[,4]

# rspfit$phase<-rspfit$phase+pi
# rspfit$phase[rspfit$phase>pi]<-rspfit$phase[rspfit$phase>pi]-2*pi-median(diff(pbins))

#display plot
rspfit$frequency<-round(rspfit$frequency,digits=1)
rspfit$frequency<-as.ordered(rspfit$frequency)
for(f in 1:nf){
  levels(rspfit$frequency)[f]<-paste0(levels(rspfit$frequency)[f]," Hz")
}
cl<-max(abs(rspfit$fit))
p<-ggplot(rspfit,aes(x=runningspeed,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,1,by=0.5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi),labels=c(expression(paste("-",pi)),"0",expression(paste(pi))),
                     expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=runningspeed,y=phase,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=runningspeed,y=phase,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=5) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("speedphase_Alz.eps",device="eps",width=14,height=14,units="in")

#running speed theta phase interaction tiled by running speed
#run for lpmatrix
rspfit<-as.data.frame(ci[rpintsimidx,])
colnames(rspfit)<-c("low","fit","high")
rspfit$fit<-rspfit$fit + rep(ci[psimidx,2],each=nb)
rspfit$low<-rspfit$low+rspfit$fit
rspfit$high<-rspfit$high+rspfit$fit
rspfit$frequency<-rpintfit[,2]
rspfit$phase<-rpintfit[,3]
rspfit$runningspeed<-rpintfit[,4]
keep<-(round((1:nb/nb)/0.2,digits=3) %% 1) == 0
keep[1]<-TRUE
keep<-rep(keep,nb*nf)
rspfit<-rspfit[keep,]

#interpolate for each running speed
#run for "lpmatrix"
pbins<-seq(pbins[1],pbins[length(pbins)],length=250)
fbins<-seq(min(freq),max(freq),length=250)
rs<-unique(rspfit$runningspeed)
fitinterp<-interp(x=rspfit$phase[rspfit$runningspeed==rs[1]],
                  y=rspfit$frequency[rspfit$runningspeed==rs[1]],
                  rspfit$fit[rspfit$runningspeed==rs[1]],
                  xo=pbins,yo=fbins)
highinterp<-interp(x=rspfit$phase[rspfit$runningspeed==rs[1]],
                   y=rspfit$frequency[rspfit$runningspeed==rs[1]],
                   rspfit$high[rspfit$runningspeed==rs[1]],
                   xo=pbins,yo=fbins)
lowinterp<-interp(x=rspfit$phase[rspfit$runningspeed==rs[1]],
                  y=rspfit$frequency[rspfit$runningspeed==rs[1]],
                  rspfit$low[rspfit$runningspeed==rs[1]],
                  xo=pbins,yo=fbins)
rfit <- data.frame(expand.grid(phase = fitinterp$x, frequency = fitinterp$y), 
                   fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z))
for(r in 2:length(rs)){
  fitinterp<-interp(x=rspfit$phase[rspfit$runningspeed==rs[r]],
                    y=rspfit$frequency[rspfit$runningspeed==rs[r]],
                    rspfit$fit[rspfit$runningspeed==rs[r]],
                    xo=pbins,yo=fbins)
  highinterp<-interp(x=rspfit$phase[rspfit$runningspeed==rs[r]],
                     y=rspfit$frequency[rspfit$runningspeed==rs[r]],
                     rspfit$high[rspfit$runningspeed==rs[r]],
                     xo=pbins,yo=fbins)
  lowinterp<-interp(x=rspfit$phase[rspfit$runningspeed==rs[r]],
                    y=rspfit$frequency[rspfit$runningspeed==rs[r]],
                    rspfit$low[rspfit$runningspeed==rs[r]],
                    xo=pbins,yo=fbins)
  rfit <- rbind(rfit,data.frame(expand.grid(phase = fitinterp$x, frequency = fitinterp$y), 
                     fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z)))
  }
rfit$runningspeed<-round(rep(rs,each=250^2),digits = 2)

#display plot
cl<-max(abs(rfit$fit))
p<-ggplot(rfit,aes(x=phase,y=frequency,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster()+scale_x_continuous(breaks=seq(-pi,pi,pi/2),labels=c("-pi","-pi/2","0","pi/2","pi"),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,250,by=25),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  # stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  # stat_contour(mapping=aes(x=phase,y=frequency,z=high),colour="purple",breaks=c(0),alpha=1)+
  # stat_contour(mapping=aes(x=phase,y=frequency,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~runningspeed)


#time theta phase interaction session1
#run for lpmatrix
tpfit<-as.data.frame(ci[tpintsimidx,])
colnames(tpfit)<-c("low","fit","high")
tpfit$low<-tpfit$low+tpfit$fit
tpfit$high<-tpfit$high+tpfit$fit
tpfit$session<-as.factor(tpintfit[,2])
tpfit$frequency<-tpintfit[,3]
tpfit$phase<-tpintfit[,4]
tpfit$time<-tpintfit[,5]

# tpfit$diff<-tpfit$fit
# tpfit$diff[tpfit$session=="1"]<-tpfit$fit[tpfit$session=="1"]-tpfit$fit[tpfit$session=="2"]
# tpfit$diff[tpfit$session=="2"]<-tpfit$fit[tpfit$session=="1"]-tpfit$fit[tpfit$session=="3"]
# tpfit$diff[tpfit$session=="3"]<-tpfit$fit[tpfit$session=="2"]-tpfit$fit[tpfit$session=="3"]
# tpfit$fit<-tpfit$diff

# tpfit$phase<-tpfit$phase+pi
# tpfit$phase[tpfit$phase>pi]<-tpfit$phase[tpfit$phase>pi]-2*pi-median(diff(pbins))

#display plot
tpfit$frequency<-round(tpfit$frequency,digits=1)
tpfit$frequency<-as.ordered(tpfit$frequency)
for(f in 1:nf){
  levels(tpfit$frequency)[f]<-paste0(levels(tpfit$frequency)[f]," Hz")
}
cl<-max(abs(tpfit$fit))
p<-ggplot(tpfit[tpfit$session=="1",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi),
                     labels=c(expression(paste("-",pi)),"0",expression(paste(pi))),
                     expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("phasetimesession1_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time theta phase interaction session2
p<-ggplot(tpfit[tpfit$session=="2",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi),
                     labels=c(expression(paste("-",pi)),"0",expression(paste(pi))),
                     expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("phasetimesession2_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time theta phase interaction session3
p<-ggplot(tpfit[tpfit$session=="3",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi),
                     labels=c(expression(paste("-",pi)),"0",expression(paste(pi))),
                     expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("phasetimesession3_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time/theta phase interactions by frequency session 1-2
tpfit<-as.data.frame(dci[tpintsimidx,])
colnames(tpfit)<-c("low","fit","high")
tpfit$low<-tpfit$low+tpfit$fit
tpfit$high<-tpfit$high+tpfit$fit
tpfit$session<-as.factor(tpintfit[,2])
tpfit$frequency<-tpintfit[,3]
tpfit$phase<-tpintfit[,4]
tpfit$time<-tpintfit[,5]

# tpfit$phase<-tpfit$phase+pi
# tpfit$phase[tpfit$phase>pi]<-tpfit$phase[tpfit$phase>pi]-2*pi-median(diff(pbins))

tpfit$frequency<-round(tpfit$frequency,digits=1)
tpfit$frequency<-as.ordered(tpfit$frequency)
for(f in 1:nf){
  levels(tpfit$frequency)[f]<-paste0(levels(tpfit$frequency)[f]," Hz")
}
cl<-max(abs(tpfit$fit))
p<-ggplot(tpfit[tpfit$session=="1",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi),
                     labels=c(expression(paste("-",pi)),"0",expression(paste(pi))),
                     expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dphasetimesession1_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time/theta phase interactions by frequency session 1-3
p<-ggplot(tpfit[tpfit$session=="2",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi),
                     labels=c(expression(paste("-",pi)),"0",expression(paste(pi))),
                     expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dphasetimesession2_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time/theta phase interactions by frequency session 2-3
p<-ggplot(tpfit[tpfit$session=="3",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi),
                     labels=c(expression(paste("-",pi)),"0",expression(paste(pi))),
                     expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dphasetimesession3_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time theta phase interaction tiled by time
#run for lpmatrix
tpfit<-as.data.frame(ci[tpintsimidx,])
colnames(tpfit)<-c("low","fit","high")
tpfit$low<-tpfit$low+tpfit$fit
tpfit$high<-tpfit$high+tpfit$fit
tpfit$session<-as.factor(tpintfit[,2])
tpfit$frequency<-tpintfit[,3]
tpfit$phase<-tpintfit[,4]
tpfit$time<-tpintfit[,5]
tpfit<-tpfit[tpfit$session=="1",]
tpfit$fit<-tpfit$fit + rep(ci[psimidx,2],each=nb)
keep<-(round((1:nb/nb)/0.2,digits=3) %% 1) == 0
keep[1]<-TRUE
keep<-rep(keep,nb*nf)
tpfit<-tpfit[keep,]

#interpolate for each time 
#run for "lpmatrix"
pbins<-seq(pbins[1],pbins[length(pbins)],length=250)
fbins<-seq(min(freq),max(freq),length=250)
t<-unique(tpfit$time)
fitinterp<-interp(x=tpfit$phase[tpfit$time==t[1]],
                  y=tpfit$frequency[tpfit$time==t[1]],
                  tpfit$fit[tpfit$time==t[1]],
                  xo=pbins,yo=fbins)
highinterp<-interp(x=tpfit$phase[tpfit$time==t[1]],
                   y=tpfit$frequency[tpfit$time==t[1]],
                   tpfit$high[tpfit$time==t[1]],
                   xo=pbins,yo=fbins)
lowinterp<-interp(x=tpfit$phase[tpfit$time==t[1]],
                  y=tpfit$frequency[tpfit$time==t[1]],
                  tpfit$low[tpfit$time==t[1]],
                  xo=pbins,yo=fbins)
tfit <- data.frame(expand.grid(phase = fitinterp$x, frequency = fitinterp$y), 
                   fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z))
for(r in 2:length(t)){
  fitinterp<-interp(x=tpfit$phase[tpfit$time==t[r]],
                    y=tpfit$frequency[tpfit$time==t[r]],
                    tpfit$fit[tpfit$time==t[r]],
                    xo=pbins,yo=fbins)
  highinterp<-interp(x=tpfit$phase[tpfit$time==t[r]],
                     y=tpfit$frequency[tpfit$time==t[r]],
                     tpfit$high[tpfit$time==t[r]],
                     xo=pbins,yo=fbins)
  lowinterp<-interp(x=tpfit$phase[tpfit$time==t[r]],
                    y=tpfit$frequency[tpfit$time==t[r]],
                    tpfit$low[tpfit$time==t[r]],
                    xo=pbins,yo=fbins)
  tfit <- rbind(tfit,data.frame(expand.grid(phase = fitinterp$x, frequency = fitinterp$y), 
                                fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z)))
}
tfit$time<-round(rep(t,each=250^2),digits = 2)

#display plot
cl<-max(abs(tfit$fit))
p<-ggplot(tfit,aes(x=phase,y=frequency,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster()+scale_x_continuous(breaks=seq(-pi,pi,pi/2),labels=c("-pi","-pi/2","0","pi/2","pi"),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,250,by=25),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  # stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  # stat_contour(mapping=aes(x=phase,y=frequency,z=high),colour="purple",breaks=c(0),alpha=1)+
  # stat_contour(mapping=aes(x=phase,y=frequency,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~time)