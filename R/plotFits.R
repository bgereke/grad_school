library(mgcv)
library(RevoUtilsMath)
library(ggplot2)
library(MASS)
library(itsadug)
library(RColorBrewer)
library(directlabels)
library(scales)
library(akima)

# if(require("RevoUtilsMath")){
#   setMKLthreads(3)
#   print("threads set")
# }

setwd("C:/Data/GAMS_CV")

#load my data
cwd<-getwd()
# setwd(file.path(cwd, "GAMS"))
# MYDATA<-read.table("MYDATA.csv", header = TRUE, sep = ",", row.names = 1)
MYDATA<-readRDS("MYDATAPP_full.rds ")
#MYDATA<-MYDATA[c(101:150,(ncol(MYDATA)-10):ncol(MYDATA))] #CICS only
#MYDATA<-MYDATA[-c(51:(ncol(MYDATA)-11))] #PICS only
# MYDATA<-MYDATA[c(201:250,(ncol(MYDATA)-10):ncol(MYDATA))] #PP only

#Set model prediction vars
numdays<-length(unique(MYDATA$Day))
days<-seq(1,numdays)
nb<-60
nf<-50
mint<-0
maxt<-600
freq<-10^seq(log10(2),log10(100),length.out=nf)
tbins<-seq(mint,maxt,length.out=nb)
minrs<-min(MYDATA$Running_Speed)
maxrs<-max(MYDATA$Running_Speed)
rsbins<-seq(minrs,maxrs,length.out=nb)
minp<- -pi
maxp<- pi
pbins<-seq(minp,maxp,length.out=nb)


#create new data
dgrid<-expand.grid(rsbins,tbins,unique(MYDATA$Session))
colnames(dgrid)<-c("runningspeed","time","session")
newdat<-with(MYDATA,data.frame(Running_Speed=dgrid$runningspeed,Time=dgrid$time,Session=dgrid$session,Trial=rep("1",length(dgrid[,1]))))

#preallocate plot data for iterms
lags<-30
acf<-matrix(0,nrow=nf,ncol=lags)
acfrho<-acf
rhofit<-matrix(0,nrow=nf,ncol=1)
ifit<-matrix(0,nrow=nf,ncol=3)
rsfit<-matrix(0,nrow=nb*nf,ncol=3)
timefit<-matrix(0,nrow=nb*nf*3,ncol=4)
intfit<-matrix(0,nrow=nb*nb*nf*3,ncol=5)
fRsq<-matrix(0,nrow=1,ncol=nf)
rsse<-matrix(0,nrow=nb*nf,ncol=1)
timese<-matrix(0,nrow=nb*nf*3,ncol=1)
intse<-matrix(0,nrow=nb*nb*nf*3,ncol=1)

#set known values
rsfit[,2]<-rep(freq,each=nb) #frequency
rsfit[,3]<-rep(rsbins,times=nf) #runnings speed
timefit[,2]<-rep(unique(MYDATA$Session),each=nb*nf) #session
timefit[,3]<-rep(rep(freq,each=nb),times=3) #frequency
timefit[,4]<-rep(tbins,times=3*nf) #time
intfit[,2]<-rep(unique(MYDATA$Session),each=nb*nb*nf) #session
intfit[,3]<-rep(rep(freq,each=nb*nb),times=3) #frequency
intfit[,4]<-rep(rep(tbins,each=nb),times=nf*3) #time
intfit[,5]<-rep(rsbins,times=nb*nf*3) #running speed

setwd("C:/Data/GAMS_full_rho")

#prediction using "iterms"
for(f in 1:nf)
{
  filename<-paste0("VTbSIbSrho_f",f,"_full_PP.rds")
  fmod<-readRDS(filename)
  fRsq[f]<-summary(fmod)$r.sq
  # pred<-predict.bam(fmod, newdat, type = "iterms",se.fit = TRUE,
  #                   newdata.guaranteed = TRUE)
  # 
  # #define indices
  # pi1idx<-1
  # pi2idx<-nb*nb+1
  # pi3idx<-2*nb*nb+1
  # pridx<-1:nb
  # pt1idx<-seq(1,nb*nb,by=nb)
  # pt2idx<-seq(nb*nb+1,2*nb*nb,by=nb)
  # pt3idx<-seq(2*nb*nb+1,3*nb*nb,by=nb)
  # pint1idx<-1:(nb*nb)
  # pint2idx<-(nb*nb+1):(2*nb*nb)
  # pint3idx<-(2*nb*nb+1):(3*nb*nb)
  # ridx<-((f-1)*nb+1):(f*nb)
  # t1idx<-((f-1)*nb+1):(f*nb)
  # t2idx<-((f-1)*nb+1+nb*nf):(f*nb+nb*nf)
  # t3idx<-((f-1)*nb+1+2*nb*nf):(f*nb+2*nb*nf)
  # int1idx<-((f-1)*nb*nb+1):(f*nb*nb)
  # int2idx<-((f-1)*nb*nb+1+nb*nb*nf):(f*nb*nb+nf*nb*nb)
  # int3idx<-((f-1)*nb*nb+1+2*nb*nb*nf):(f*nb*nb+2*nf*nb*nb)
  # 
  # #rho's and intercepts
  # rhofit[f]<-fmod$AR1.rho
  # ifit[f,1]<-pred$fit[pi1idx,1] #session 1
  # ifit[f,2]<-pred$fit[pi2idx,1] #session 2
  # ifit[f,3]<-pred$fit[pi3idx,1] #session 3
  # 
  # #running speed 
  # rsfit[ridx,1]<-pred$fit[pridx,2]
  # rsse[ridx]<-pred$se.fit[pridx,2]
  # 
  # #time by session 
  # timefit[t1idx,1]<-pred$fit[pt1idx,3] #session 1
  # timefit[t2idx,1]<-pred$fit[pt2idx,4] #session 2
  # timefit[t3idx,1]<-pred$fit[pt3idx,5] #session 3
  # timese[t1idx]<-pred$se.fit[pt1idx,3] #session 1
  # timese[t2idx]<-pred$se.fit[pt2idx,4] #session 2
  # timese[t3idx]<-pred$se.fit[pt3idx,5] #session 3
  # 
  # #time/running speed interactions by session
  # intfit[int1idx,1]<-pred$fit[pint1idx,6] #session 1
  # intfit[int2idx,1]<-pred$fit[pint2idx,7] #session 2
  # intfit[int3idx,1]<-pred$fit[pint3idx,8] #session 3
  # intse[int1idx]<-pred$se.fit[pint1idx,6] #session 1
  # intse[int2idx]<-pred$se.fit[pint2idx,7] #session 2
  # intse[int3idx]<-pred$se.fit[pint3idx,8] #session 3
  
  print(paste("Completed: ",f," of ",nf),quote=FALSE)
}

#prediction using "lpmatrix"
ns<-1000
isim<-matrix(0,nrow = 3*nf,ncol = ns)
rsim<-matrix(0,nrow = nf*nb,ncol = ns)
timesim<-matrix(0,nrow = 3*nf*nb,ncol = ns)
intsim<-matrix(0,nrow = 3*nf*nb*nb,ncol = ns)
for(f in 1:nf)
{
  filename<-paste0("VTbSIbSrho_f",f,"_full_PP.rds")
  fmod<-readRDS(filename)
  # filename<-paste0("VTbSIbSrho_f",f,"_full_PICS.rds")
  # fmodrho<-readRDS(filename)
  pred<-predict.bam(fmod, newdat, type = "lpmatrix", newdata.guaranteed = TRUE)
  
  #rho's
  # res<-resid_gam(fmod,incl_na = TRUE,return_all = TRUE)
  # acf[f,]<-acf_resid(fmod,split_pred = c("Trial"),plot=FALSE,max_lag=lags+1)[-1]
  # acfrho[f,]<-acf_resid(fmodrho,split_pred = c("Trial"),plot=FALSE,max_lag=lags+1)[-1]
  # rhofit[f]<-fmodrho$AR1.rho
  
  #get mean coefficient values and var-ance-covariance matric
  coefs<-coef(fmod)
  vc<-vcov(fmod)
  #draw random coefficient values from the model
  set.seed(35)
  sim<-mvrnorm(ns,mu=coefs,Sigma=vc)

  rm(coefs,vc,fmod)
  gc()

  #define indices
  pi1idx<-1
  pi2idx<-nb*nb+1
  pi3idx<-2*nb*nb+1
  pridx<-1:nb
  pt1idx<-seq(1,nb*nb,by=nb)
  pt2idx<-seq(nb*nb+1,2*nb*nb,by=nb)
  pt3idx<-seq(2*nb*nb+1,3*nb*nb,by=nb)
  pint1idx<-1:(nb*nb)
  pint2idx<-(nb*nb+1):(2*nb*nb)
  pint3idx<-(2*nb*nb+1):(3*nb*nb)
  ridx<-((f-1)*nb+1):(f*nb)
  t1idx<-((f-1)*nb+1):(f*nb)
  t2idx<-((f-1)*nb+1+nb*nf):(f*nb+nb*nf)
  t3idx<-((f-1)*nb+1+2*nb*nf):(f*nb+2*nb*nf)
  int1idx<-((f-1)*nb*nb+1):(f*nb*nb)
  int2idx<-((f-1)*nb*nb+1+nb*nb*nf):(f*nb*nb+nf*nb*nb)
  int3idx<-((f-1)*nb*nb+1+2*nb*nb*nf):(f*nb*nb+2*nf*nb*nb)

  #intercepts
  isim[f,] <- pred[pi2idx, 1] %*% t(sim[, 1]) #session 1
  isim[nf+f,] <- pred[pi2idx, 2] %*% t(sim[, 2]) #session 2
  isim[2*nf+f,] <- pred[pi3idx, 3] %*% t(sim[, 3]) #session 3

  #runningspeed
  want<-grep("s\\(Running\\_Speed\\)",colnames(pred))
  rsim[ridx,]<-pred[pridx,want] %*% t(sim[,want])

  #time by session
  want1<-grep("s\\(Time\\)\\:Session1",colnames(pred))
  want2<-grep("s\\(Time\\)\\:Session2",colnames(pred))
  want3<-grep("s\\(Time\\)\\:Session3",colnames(pred))
  timesim[t1idx,] <- pred[pt1idx, want1] %*% t(sim[, want1]) #session 1
  timesim[t2idx,] <- pred[pt2idx, want2] %*% t(sim[, want2]) #session 2
  timesim[t3idx,] <- pred[pt3idx, want3] %*% t(sim[, want3]) #session 3

  #time/running speed interactions
  want1<-grep("ti\\(Time\\,Running\\_Speed\\)\\:Session1",colnames(pred))
  want2<-grep("ti\\(Time\\,Running\\_Speed\\)\\:Session2",colnames(pred))
  want3<-grep("ti\\(Time\\,Running\\_Speed\\)\\:Session3",colnames(pred))
  intsim[int1idx,]<-pred[pint1idx,want1] %*% t(sim[,want1]) #session 1
  intsim[int2idx,]<-pred[pint2idx,want2] %*% t(sim[,want2]) #session 2
  intsim[int3idx,]<-pred[pint3idx,want3] %*% t(sim[,want3]) #session 3

  rm(pred,sim)
  gc()
  
  print(paste("Completed: ",f," of ",nf),quote=FALSE)
  
}

#run for "lpmatrix" only
#get simultaneous CI's for main terms
source("C:/Users/Brian/Documents/R/simCIz.R")
ci<-simCIz(rbind(isim,rsim,timesim,intsim))
gc()
#get simultaneous CI's for difference terms
dtimesim<-matrix(0,nrow = nrow(timesim),ncol = ncol(timesim))
dtimesim[1:(nf*nb),]<-timesim[1:(nf*nb),]-timesim[(1+nb*nf):(2*nb*nf),]
dtimesim[(1+nb*nf):(2*nb*nf),]<-timesim[1:(nf*nb),]-timesim[(1+2*nb*nf):(3*nb*nf),]
dtimesim[(1+2*nb*nf):(3*nb*nf),]<-timesim[(1+nb*nf):(2*nb*nf),]-timesim[(1+2*nb*nf):(3*nb*nf),]
dintsim<-matrix(0,nrow = nrow(intsim),ncol(intsim))
dintsim[1:(nf*nb*nb),]<-intsim[1:(nf*nb*nb),]-intsim[(1+nb*nb*nf):(2*nb*nb*nf),]
dintsim[(1+nb*nb*nf):(2*nb*nb*nf),]<-intsim[1:(nf*nb*nb),]-intsim[(1+2*nb*nb*nf):(3*nb*nb*nf),]
dintsim[(1+2*nb*nb*nf):(3*nb*nb*nf),]<-intsim[(1+nb*nb*nf):(2*nb*nb*nf),]-intsim[(1+2*nb*nb*nf):(3*nb*nb*nf),]
gc()
dci<-simCIz(rbind(isim,rsim,dtimesim,dintsim))
gc()

#define indices for plot terms
isimidx<-1:nrow(isim) #intercepts
rsimidx<-(nrow(isim)+1):(nrow(isim)+nrow(rsim)) #running speed
timesimidx<-(nrow(isim)+nrow(rsim)+1):(nrow(isim)+nrow(rsim)+nrow(timesim)) #time by session
intsimidx<-(nrow(isim)+nrow(rsim)+nrow(timesim)+1):(nrow(isim)+nrow(rsim)+nrow(timesim)+nrow(intsim)) #time/running speed interaction by session

# rm(isim,timesim,intsim,disim,dtimesim,dintsim)
# gc()

#r-squared plot
rsqdf<-data.frame(Frequency=freq,Rsquared=as.numeric(fRsq))
p<-ggplot(rsqdf,aes(x=Frequency,y=Rsquared))
p+geom_line()+geom_path(size=1)+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,0.2,by=0.025),name="R^2")+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())


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
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
cl<-1
ggplot(acfdf,aes(x=lags,y=frequency,fill=acf,z=acf))+ 
  geom_tile()+scale_x_continuous(breaks=c(1,10,20,30),labels=c(1,10,20,30),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",alpha=0.3,binwidth=0.05,size=0.05,show.legend = TRUE)+
  facet_wrap(~rho,nrow=1)+theme(aspect.ratio = 1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())


#running speed plot
#run for iterms
rbins<-seq(rsbins[1],rsbins[length(rsbins)],length=250)
fbins<-seq(2,100,length=250)
fitinterp<-interp(x=rsfit[,3],y=rsfit[,2],rsfit[,1],xo=rbins,yo=fbins)
highinterp<-interp(x=rsfit[,3],y=rsfit[,2],rsfit[,1]+2*rsse[,1],xo=rbins,yo=fbins)
lowinterp<-interp(x=rsfit[,3],y=rsfit[,2],rsfit[,1]-2*rsse[,1],xo=rbins,yo=fbins)
rfit <- data.frame(expand.grid(runningspeed = fitinterp$x, frequency = fitinterp$y), 
                   fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z))

#run for "lpmatrix"
rbins<-seq(rsbins[1],rsbins[length(rsbins)],length=250)
fbins<-seq(2,100,length=250)
fitinterp<-interp(x=rsfit[,3],y=rsfit[,2],ci[rsimidx,2],xo=rbins,yo=fbins)
highinterp<-interp(x=rsfit[,3],y=rsfit[,2],ci[rsimidx,2]+ci[rsimidx,3],xo=rbins,yo=fbins)
lowinterp<-interp(x=rsfit[,3],y=rsfit[,2],ci[rsimidx,2]+ci[rsimidx,1],xo=rbins,yo=fbins)
rfit <- data.frame(expand.grid(runningspeed = fitinterp$x, frequency = fitinterp$y), 
                   fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z))

#display plot
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
cl<-max(abs(rfit$fit))
ggplot(rfit,aes(x=runningspeed,y=frequency,fill=fit,z=fit))+
  geom_tile()+scale_x_continuous(breaks=seq(round(minrs),round(maxrs),by=0.2),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  stat_contour(mapping=aes(x=runningspeed,y=frequency,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=runningspeed,y=frequency,z=low),colour="orange",breaks=c(0),alpha=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1)


#time by session
#run for iterms
tfit<-as.data.frame(timefit)
colnames(tfit)<-c("fit","session","frequency","time")
tfit$fit<-tfit$fit#+rep(c(ifit[,1],ifit[,2],ifit[,3]),each=nb)
tfit$low<-tfit$fit - 2*timese[,1] 
tfit$high<-tfit$fit + 2*timese[,1]

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
fbins<-seq(2,100,length=250)
idxone<-tfit$session == "1"
idxtwo<-tfit$session == "2"
idxthree<-tfit$session == "3"
fitsone<-interp(x=tfit$time[idxone],y=tfit$frequency[idxone],tfit$fit[idxone],xo=tibins,yo=fbins)
fitstwo<-interp(x=tfit$time[idxtwo],y=tfit$frequency[idxtwo],tfit$fit[idxtwo],xo=tibins,yo=fbins)
fitsthree<-interp(x=tfit$time[idxthree],y=tfit$frequency[idxthree],tfit$fit[idxthree],xo=tibins,yo=fbins)
highone<-interp(x=tfit$time[idxone],y=tfit$frequency[idxone],tfit$high[idxone],xo=tibins,yo=fbins)
hightwo<-interp(x=tfit$time[idxtwo],y=tfit$frequency[idxtwo],tfit$high[idxtwo],xo=tibins,yo=fbins)
highthree<-interp(x=tfit$time[idxthree],y=tfit$frequency[idxthree],tfit$high[idxthree],xo=tibins,yo=fbins)
lowone<-interp(x=tfit$time[idxone],y=tfit$frequency[idxone],tfit$low[idxone],xo=tibins,yo=fbins)
lowtwo<-interp(x=tfit$time[idxtwo],y=tfit$frequency[idxtwo],tfit$low[idxtwo],xo=tibins,yo=fbins)
lowthree<-interp(x=tfit$time[idxthree],y=tfit$frequency[idxthree],tfit$low[idxthree],xo=tibins,yo=fbins)
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

#display plot
cl<-max(abs(tfit$fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
ggplot(tfit,aes(x=time,y=frequency,fill=fit,z=fit))+ 
  geom_tile()+scale_x_continuous(breaks=seq(0,600,by=60),labels=seq(0,10,by=1),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",alpha=0.3,binwidth=0.05,size=0.05,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=frequency,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=frequency,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~session,nrow=1)+theme(aspect.ratio = 1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())

#differences time by session
#run for iterms (must run time by session plot first and can't get CI's)
tfit$diffs<-tfit$fit
tfit$diffs[tfit$session=="1"]<-tfit$fit[tfit$session=="1"] - tfit$fit[tfit$session=="2"]
tfit$diffs[tfit$session=="2"]<-tfit$fit[tfit$session=="1"] - tfit$fit[tfit$session=="3"]
tfit$diffs[tfit$session=="3"]<-tfit$fit[tfit$session=="2"] - tfit$fit[tfit$session=="3"]

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
fbins<-seq(2,100,length=250)
idxone<-tfit$session == "1"
idxtwo<-tfit$session == "2"
idxthree<-tfit$session == "3"
fitsone<-interp(x=tfit$time[idxone],y=tfit$frequency[idxone],tfit$fit[idxone],xo=tibins,yo=fbins)
fitstwo<-interp(x=tfit$time[idxtwo],y=tfit$frequency[idxtwo],tfit$fit[idxtwo],xo=tibins,yo=fbins)
fitsthree<-interp(x=tfit$time[idxthree],y=tfit$frequency[idxthree],tfit$fit[idxthree],xo=tibins,yo=fbins)
highone<-interp(x=tfit$time[idxone],y=tfit$frequency[idxone],tfit$high[idxone],xo=tibins,yo=fbins)
hightwo<-interp(x=tfit$time[idxtwo],y=tfit$frequency[idxtwo],tfit$high[idxtwo],xo=tibins,yo=fbins)
highthree<-interp(x=tfit$time[idxthree],y=tfit$frequency[idxthree],tfit$high[idxthree],xo=tibins,yo=fbins)
lowone<-interp(x=tfit$time[idxone],y=tfit$frequency[idxone],tfit$low[idxone],xo=tibins,yo=fbins)
lowtwo<-interp(x=tfit$time[idxtwo],y=tfit$frequency[idxtwo],tfit$low[idxtwo],xo=tibins,yo=fbins)
lowthree<-interp(x=tfit$time[idxthree],y=tfit$frequency[idxthree],tfit$low[idxthree],xo=tibins,yo=fbins)
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

#display plot
tfit$dsessions<-as.factor(c(rep("1-2",times=250*250),rep("1-3",times=250*250),rep("2-3",times=250*250)))
cl<-max(abs(tfit$diffs))
p<-ggplot(tfit,aes(x=time,y=frequency,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_tile()+scale_x_continuous(breaks=seq(0,600,by=60),labels=seq(0,10,by=1),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",alpha=0.3,binwidth=0.05,size=0.05,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=frequency,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=frequency,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~dsessions,nrow=1)+theme(aspect.ratio = 1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())  

#time/running speed interactions by frequency session 1
trsfit<-as.data.frame(intfit)
colnames(trsfit)<-c("fit","session","frequency","time","runningspeed")
trsfit$session<-as.factor(trsfit$session)
trsfit$high<-trsfit$fit+2*intse[1,]
trsfit$low<-trsfit$fit-2*intse[1,]

#run for lpmatrix
trsfit<-as.data.frame(ci[intsimidx,])
colnames(trsfit)<-c("low","fit","high")
trsfit$low<-trsfit$low+trsfit$fit
trsfit$high<-trsfit$high+trsfit$fit
trsfit$session<-intfit[,2]
trsfit$frequency<-intfit[,3]
trsfit$time<-intfit[,4]
trsfit$runningspeed<-intfit[,5]

#display plot
trsfit$frequency<-round(trsfit$frequency,digits=1)
trsfit$frequency<-as.ordered(trsfit$frequency)
for(f in 1:nf){
  levels(trsfit$frequency)[f]<-paste0(levels(trsfit$frequency)[f]," Hz")
}
trsfit<-trsfit[trsfit$runningspeed<=0.75,]
cl<-0.7*max(abs(trsfit$fit))
p<-ggplot(trsfit[trsfit$session=="1",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=5) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))

#time/running speed interactions by frequency session 2
p<-ggplot(trsfit[trsfit$session=="2",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=5) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))

#time/running speed interactions by frequency session 3
p<-ggplot(trsfit[trsfit$session=="3",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=5) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))

#time/running speed interactions by frequency session 1-2
trsfit<-as.data.frame(dci[intsimidx,])
colnames(trsfit)<-c("low","fit","high")
trsfit$low<-trsfit$low+trsfit$fit
trsfit$high<-trsfit$high+trsfit$fit
trsfit$session<-rep(c("1","2","3"),each=nb*nb*nf)
trsfit$frequency<-rep(rep(freq,each=nb*nb),times=3)
trsfit$time<-rep(rep(seq(mint,maxt,length.out=nb),each=nb),times=nf*3)
trsfit$runningspeed<-rep(rsbins,times=nb*nf*3)

trsfit$frequency<-round(trsfit$frequency,digits=1)
trsfit$frequency<-as.ordered(trsfit$frequency)
for(f in 1:nf){
  levels(trsfit$frequency)[f]<-paste0(levels(trsfit$frequency)[f]," Hz")
}
trsfit<-trsfit[trsfit$runningspeed<=0.75,]
cl<-0.7*max(abs(trsfit$fit))
p<-ggplot(trsfit[trsfit$session=="1",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=5) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))

#time/running speed interactions by frequency session 1-3
cl<-0.7*max(abs(trsfit$fit))
p<-ggplot(trsfit[trsfit$session=="2",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=5) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))

#time/running speed interactions by frequency session 2-3
cl<-0.7*max(abs(trsfit$fit))
p<-ggplot(trsfit[trsfit$session=="3",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=5) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))