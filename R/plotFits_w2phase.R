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

setwd("C:/Data/GAMS_full")

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
freq<-freq[30:nf]
nf<-21
tbins<-seq(mint,maxt,length.out=nb)
minrs<-min(MYDATA$Running_Speed)
maxrs<-max(MYDATA$Running_Speed)
rsbins<-seq(minrs,maxrs,length.out=nb)
minp<- -pi
maxp<- pi
pbins<-seq(minp,maxp,length.out=nb)


#create new data
dgrid<-expand.grid(rsbins,tbins,pbins,pbins,unique(MYDATA$Session))
colnames(dgrid)<-c("runningspeed","time","stphase","ftphase","session")
newdat<-with(MYDATA,data.frame(Running_Speed=dgrid$runningspeed,Time=dgrid$time,STphase=dgrid$stphase,FTphase=dgrid$ftphase,
                               Session=dgrid$session,Trial=rep("1",length(dgrid[,1])),Day=rep("1",length(dgrid[,1]))))

#define indices and remove unused predictions
pi1idx<-which(!duplicated(newdat[,c('Session')]))[1]
pi2idx<-which(!duplicated(newdat[,c('Session')]))[2]
pi3idx<-which(!duplicated(newdat[,c('Session')]))[3]
pridx<-which(!duplicated(newdat[,c('Running_Speed')]))
ptidx<-which(!duplicated(newdat[,c('Time','Session')]))
pstidx<-which(!duplicated(newdat[,c('STphase')]))
pftidx<-which(!duplicated(newdat[,c('FTphase')]))
prtintidx<-which(!duplicated(newdat[,c('Running_Speed','Time','Session')])) 
prstintidx<-which(!duplicated(newdat[,c('Running_Speed','STphase')]))
ptstintidx<-which(!duplicated(newdat[,c('Time','STphase','Session')])) 
prftintidx<-which(!duplicated(newdat[,c('Running_Speed','FTphase')]))
ptftintidx<-which(!duplicated(newdat[,c('Time','FTphase','Session')])) 
used<-c(pi1idx,pi2idx,pi3idx,pridx,ptidx,pstidx,pftidx,prtintidx,
        prstintidx,prftintidx,ptstintidx,ptftintidx)
newdat<-newdat[sort(unique(used)),]

#define new indices
pi1idx<-which(!duplicated(newdat[,c('Session')]))[1]
pi2idx<-which(!duplicated(newdat[,c('Session')]))[2]
pi3idx<-which(!duplicated(newdat[,c('Session')]))[3]
pridx<-which(!duplicated(newdat[,c('Running_Speed')]))
ptidx<-which(!duplicated(newdat[,c('Time','Session')]))
pstidx<-which(!duplicated(newdat[,c('STphase')]))
pftidx<-which(!duplicated(newdat[,c('FTphase')]))
prtintidx<-which(!duplicated(newdat[,c('Running_Speed','Time','Session')])) 
prstintidx<-which(!duplicated(newdat[,c('Running_Speed','STphase')]))
ptstintidx<-which(!duplicated(newdat[,c('Time','STphase','Session')])) 
prftintidx<-which(!duplicated(newdat[,c('Running_Speed','FTphase')]))
ptftintidx<-which(!duplicated(newdat[,c('Time','FTphase','Session')]))  


#preallocate plot data for iterms
lags<-30
acf<-matrix(0,nrow=nf,ncol=lags)
acfrho<-acf
rhofit<-matrix(0,nrow=nf,ncol=1)
ifit<-matrix(0,nrow=nf,ncol=3)
rsfit<-matrix(0,nrow=nb*nf,ncol=3)
stfit<-matrix(0,nrow=nb*nf,ncol=3)
ftfit<-matrix(0,nrow=nb*nf,ncol=3)
timefit<-matrix(0,nrow=nb*nf*3,ncol=4)
rtintfit<-matrix(0,nrow=nb*nb*nf*3,ncol=5)
rstintfit<-matrix(0,nrow=nb*nb*nf,ncol=4)
rftintfit<-matrix(0,nrow=nb*nb*nf,ncol=4)
tstintfit<-matrix(0,nrow=nb*nb*nf*3,ncol=5)
tftintfit<-matrix(0,nrow=nb*nb*nf*3,ncol=5)
fRsq<-matrix(0,nrow=1,ncol=nf)

#set known values
rsfit[,2]<-rep(freq,each=nb) #frequency
rsfit[,3]<-rep(rsbins,times=nf) #runnings speed
stfit[,2]<-rep(freq,each=nb) #frequency
stfit[,3]<-rep(pbins,times=nf) #stphase
ftfit[,2]<-rep(freq,each=nb) #frequency
ftfit[,3]<-rep(pbins,times=nf) #ftphase
timefit[,2]<-rep(unique(MYDATA$Session),each=nb*nf) #session
timefit[,3]<-rep(rep(freq,each=nb),times=3) #frequency
timefit[,4]<-rep(tbins,times=3*nf) #time
rtintfit[,2]<-rep(unique(MYDATA$Session),each=nb*nb*nf) #session
rtintfit[,3]<-rep(rep(freq,each=nb*nb),times=3) #frequency
rtintfit[,4]<-rep(rep(tbins,each=nb),times=nf*3) #time
rtintfit[,5]<-rep(rsbins,times=nb*nf*3) #running speed
rstintfit[,2]<-rep(freq,each=nb*nb) #frequency
rstintfit[,3]<-rep(rep(pbins,each=nb),times=nf) #stphase
rstintfit[,4]<-rep(rsbins,times=nb*nf) #running speed
rftintfit[,2]<-rep(freq,each=nb*nb) #frequency
rftintfit[,3]<-rep(rep(pbins,each=nb),times=nf) #ftphase
rftintfit[,4]<-rep(rsbins,times=nb*nf) #running speed
tstintfit[,2]<-rep(unique(MYDATA$Session),each=nb*nb*nf) #session
tstintfit[,3]<-rep(rep(freq,each=nb*nb),times=3) #frequency
tstintfit[,4]<-rep(rep(pbins,each=nb),times=nf*3) #stphase
tstintfit[,5]<-rep(tbins,times=nb*nf*3) #time
tftintfit[,2]<-rep(unique(MYDATA$Session),each=nb*nb*nf) #session
tftintfit[,3]<-rep(rep(freq,each=nb*nb),times=3) #frequency
tftintfit[,4]<-rep(rep(pbins,each=nb),times=nf*3) #ftphase
tftintfit[,5]<-rep(tbins,times=nb*nf*3) #time

setwd("C:/Data/GAMS_stft")

#prediction using "lpmatrix"
ns<-1000
isim<-matrix(0,nrow = 3*nf,ncol = ns)
rsim<-matrix(0,nrow = nf*nb,ncol = ns)
stsim<-matrix(0,nrow = nf*nb,ncol = ns)
ftsim<-matrix(0,nrow = nf*nb,ncol = ns)
timesim<-matrix(0,nrow = 3*nf*nb,ncol = ns)
rtintsim<-matrix(0,nrow = 3*nf*nb*nb,ncol = ns)
rstintsim<-matrix(0,nrow = nf*nb*nb,ncol = ns)
tstintsim<-matrix(0,nrow = 3*nf*nb*nb,ncol = ns)
rftintsim<-matrix(0,nrow = nf*nb*nb,ncol = ns)
tftintsim<-matrix(0,nrow = 3*nf*nb*nb,ncol = ns)
fnum<-30:50
for(f in 1:nf)
{
  filename<-paste0("VTP_bS_rho_f",fnum[f],"_full_PP.rds")
  fmod<-readRDS(filename)
  # filename<-paste0("VTbSIbSrho_f",f,"_full_PICS.rds")
  # fmodrho<-readRDS(filename)
  
  #rho's
  # res<-resid_gam(fmod,incl_na = TRUE,return_all = TRUE)
  # acf[f,]<-acf_resid(fmod,split_pred = c("Trial"),plot=FALSE,max_lag=lags+1)[-1]
  # acfrho[f,]<-acf_resid(fmodrho,split_pred = c("Trial"),plot=FALSE,max_lag=lags+1)[-1]
  rhofit[f]<-fmod$AR1.rho
  # fRsq[f]<-summary(fmod)$r.sq
  
  #get mean coefficient values and var-ance-covariance matric
  coefs<-coef(fmod)
  vc<-vcov(fmod)
  #draw random coefficient values from the model
  set.seed(35)
  sim<-mvrnorm(ns,mu=coefs,Sigma=vc)

  #get predictions
  pred<-predict.bam(fmod, newdat[,], type = "lpmatrix", newdata.guaranteed = TRUE)

  rm(coefs,vc,fmod)
  gc()

  #define indices
  pt1idx<-ptidx[newdat$Session[ptidx]=="1"]
  pt2idx<-ptidx[newdat$Session[ptidx]=="2"]
  pt3idx<-ptidx[newdat$Session[ptidx]=="3"]
  prtint1idx<-prtintidx[newdat$Session[prtintidx]=="1"]
  prtint2idx<-prtintidx[newdat$Session[prtintidx]=="2"]
  prtint3idx<-prtintidx[newdat$Session[prtintidx]=="3"]
  ptstint1idx<-ptstintidx[newdat$Session[ptstintidx]=="1"]
  ptstint2idx<-ptstintidx[newdat$Session[ptstintidx]=="2"]
  ptstint3idx<-ptstintidx[newdat$Session[ptstintidx]=="3"]
  ptftint1idx<-ptftintidx[newdat$Session[ptftintidx]=="1"]
  ptftint2idx<-ptftintidx[newdat$Session[ptftintidx]=="2"]
  ptftint3idx<-ptftintidx[newdat$Session[ptftintidx]=="3"]

  ridx<-((f-1)*nb+1):(f*nb)
  stidx<-((f-1)*nb+1):(f*nb)
  ftidx<-((f-1)*nb+1):(f*nb)
  t1idx<-((f-1)*nb+1):(f*nb)
  t2idx<-((f-1)*nb+1+nb*nf):(f*nb+nb*nf)
  t3idx<-((f-1)*nb+1+2*nb*nf):(f*nb+2*nb*nf)
  rtint1idx<-((f-1)*nb*nb+1):(f*nb*nb)
  rtint2idx<-((f-1)*nb*nb+1+nb*nb*nf):(f*nb*nb+nf*nb*nb)
  rtint3idx<-((f-1)*nb*nb+1+2*nb*nb*nf):(f*nb*nb+2*nf*nb*nb)
  rstintidx<-((f-1)*nb*nb+1):(f*nb*nb)
  tstint1idx<-((f-1)*nb*nb+1):(f*nb*nb)
  tstint2idx<-((f-1)*nb*nb+1+nb*nb*nf):(f*nb*nb+nf*nb*nb)
  tstint3idx<-((f-1)*nb*nb+1+2*nb*nb*nf):(f*nb*nb+2*nf*nb*nb)
  rftintidx<-((f-1)*nb*nb+1):(f*nb*nb)
  tftint1idx<-((f-1)*nb*nb+1):(f*nb*nb)
  tftint2idx<-((f-1)*nb*nb+1+nb*nb*nf):(f*nb*nb+nf*nb*nb)
  tftint3idx<-((f-1)*nb*nb+1+2*nb*nb*nf):(f*nb*nb+2*nf*nb*nb)

  #intercepts
  isim[f,] <- pred[pi1idx, 1] %*% t(sim[, 1]) #session 1
  isim[nf+f,] <- pred[pi2idx, 2] %*% t(sim[, 2]) #session 2
  isim[2*nf+f,] <- pred[pi3idx, 3] %*% t(sim[, 3]) #session 3

  #runningspeed
  want<-grep("s\\(Running\\_Speed\\)",colnames(pred))
  rsim[ridx,]<-pred[pridx,want] %*% t(sim[,want])

  #stphase
  want<-grep("s\\(STphase\\)",colnames(pred))
  stsim[stidx,]<-pred[pstidx,want] %*% t(sim[,want])
  
  #ftphase
  want<-grep("s\\(FTphase\\)",colnames(pred))
  ftsim[ftidx,]<-pred[pftidx,want] %*% t(sim[,want])

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
  rtintsim[rtint1idx,]<-pred[prtint1idx,want1] %*% t(sim[,want1]) #session 1
  rtintsim[rtint2idx,]<-pred[prtint2idx,want2] %*% t(sim[,want2]) #session 2
  rtintsim[rtint3idx,]<-pred[prtint3idx,want3] %*% t(sim[,want3]) #session 3

  #running speed/stphase interactions
  want1<-grep("ti\\(Running\\_Speed\\,STphase\\)",colnames(pred))
  rstintsim[rstintidx,]<-pred[prstintidx,want1] %*% t(sim[,want1])
  
  #running speed/ftphase interactions
  want1<-grep("ti\\(Running\\_Speed\\,FTphase\\)",colnames(pred))
  rftintsim[rftintidx,]<-pred[prftintidx,want1] %*% t(sim[,want1])

  #time/stphase interactions
  want1<-grep("ti\\(Time\\,STphase\\)\\:Session1",colnames(pred))
  want2<-grep("ti\\(Time\\,STphase\\)\\:Session2",colnames(pred))
  want3<-grep("ti\\(Time\\,STphase\\)\\:Session3",colnames(pred))
  tstintsim[tstint1idx,]<-pred[ptstint1idx,want1] %*% t(sim[,want1]) #session 1
  tstintsim[tstint2idx,]<-pred[ptstint2idx,want2] %*% t(sim[,want2]) #session 2
  tstintsim[tstint3idx,]<-pred[ptstint3idx,want3] %*% t(sim[,want3]) #session 3
  
  #time/ftphase interactions
  want1<-grep("ti\\(Time\\,FTphase\\)\\:Session1",colnames(pred))
  want2<-grep("ti\\(Time\\,FTphase\\)\\:Session2",colnames(pred))
  want3<-grep("ti\\(Time\\,FTphase\\)\\:Session3",colnames(pred))
  tftintsim[tftint1idx,]<-pred[ptftint1idx,want1] %*% t(sim[,want1]) #session 1
  tftintsim[tftint2idx,]<-pred[ptftint2idx,want2] %*% t(sim[,want2]) #session 2
  tftintsim[tftint3idx,]<-pred[ptftint3idx,want3] %*% t(sim[,want3]) #session 3

  rm(pred,sim)
  gc()
  
  print(paste("Completed: ",f," of ",nf),quote=FALSE)
  
}
rm(newdat)
gc()

#run for "lpmatrix" only
#get simultaneous CI's for main terms
source("C:/Users/Brian/Documents/R/simCIz.R")
ci<-simCIz(rbind(isim,rsim,stsim,ftsim,timesim,rtintsim,rstintsim,rftintsim,tstintsim,tftintsim))
gc()
#get simultaneous CI's for difference terms
dtimesim<-matrix(0,nrow = nrow(timesim),ncol = ncol(timesim))
dtimesim[1:(nf*nb),]<-timesim[1:(nf*nb),]-timesim[(1+nb*nf):(2*nb*nf),]
dtimesim[(1+nb*nf):(2*nb*nf),]<-timesim[1:(nf*nb),]-timesim[(1+2*nb*nf):(3*nb*nf),]
dtimesim[(1+2*nb*nf):(3*nb*nf),]<-timesim[(1+nb*nf):(2*nb*nf),]-timesim[(1+2*nb*nf):(3*nb*nf),]
drtintsim<-matrix(0,nrow = nrow(rtintsim),ncol(rtintsim))
drtintsim[1:(nf*nb*nb),]<-rtintsim[1:(nf*nb*nb),]-rtintsim[(1+nb*nb*nf):(2*nb*nb*nf),]
drtintsim[(1+nb*nb*nf):(2*nb*nb*nf),]<-rtintsim[1:(nf*nb*nb),]-rtintsim[(1+2*nb*nb*nf):(3*nb*nb*nf),]
drtintsim[(1+2*nb*nb*nf):(3*nb*nb*nf),]<-rtintsim[(1+nb*nb*nf):(2*nb*nb*nf),]-rtintsim[(1+2*nb*nb*nf):(3*nb*nb*nf),]
dtpintsim<-matrix(0,nrow = nrow(tpintsim),ncol(tpintsim))
dtpintsim[1:(nf*nb*nb),]<-tpintsim[1:(nf*nb*nb),]-tpintsim[(1+nb*nb*nf):(2*nb*nb*nf),]
dtpintsim[(1+nb*nb*nf):(2*nb*nb*nf),]<-tpintsim[1:(nf*nb*nb),]-tpintsim[(1+2*nb*nb*nf):(3*nb*nb*nf),]
dtpintsim[(1+2*nb*nb*nf):(3*nb*nb*nf),]<-tpintsim[(1+nb*nb*nf):(2*nb*nb*nf),]-tpintsim[(1+2*nb*nb*nf):(3*nb*nb*nf),]
gc()
dci<-simCIz(rbind(isim,rsim,psim,dtimesim,drtintsim,rpintsim,dtpintsim))
gc()

#define indices for plot terms
isimidx<-1:nrow(isim) #intercepts
rsimidx<-(max(isimidx)+1):(max(isimidx)+nrow(rsim)) #running speed
stsimidx<-(max(rsimidx)+1):(max(rsimidx)+nrow(stsim)) #stphase
ftsimidx<-(max(stsimidx)+1):(max(stsimidx)+nrow(ftsim)) #ftphase
timesimidx<-(max(ftsimidx)+1):(max(ftsimidx)+nrow(timesim)) #time by session
rtintsimidx<-(max(timesimidx)+1):(max(timesimidx)+nrow(rtintsim)) #time/running speed interaction by session
rstintsimidx<-(max(rtintsimidx)+1):(max(rtintsimidx)+nrow(rstintsim)) #runningspeed/stphase interaction
rftintsimidx<-(max(rstintsimidx)+1):(max(rstintsimidx)+nrow(rftintsim)) #runningspeed/ftphase interaction
tstintsimidx<-(max(rftintsimidx)+1):(max(rftintsimidx)+nrow(tstintsim)) #time/stphase interaction
tftintsimidx<-(max(tstintsimidx)+1):(max(tstintsimidx)+nrow(tftintsim)) #time/ftphase interaction

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
fbins<-seq(min(freq),max(freq),length=250)
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

#theta phase plot
#run for "lpmatrix"
pbins<-seq(pbins[1],pbins[length(pbins)],length=250)
fbins<-seq(min(freq),max(freq),length=250)
fitinterp<-interp(x=stfit[,3],y=stfit[,2],ci[stsimidx,2],xo=pbins,yo=fbins)
highinterp<-interp(x=stfit[,3],y=stfit[,2],ci[stsimidx,2]+ci[stsimidx,3],xo=pbins,yo=fbins)
lowinterp<-interp(x=stfit[,3],y=stfit[,2],ci[stsimidx,2]+ci[stsimidx,1],xo=pbins,yo=fbins)
pfit <- data.frame(expand.grid(phase = fitinterp$x, frequency = fitinterp$y), 
                   fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z))

#display plot
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
cl<-max(abs(pfit$fit))
ggplot(pfit,aes(x=phase,y=frequency,fill=fit,z=fit))+
  geom_tile()+scale_x_continuous(breaks=seq(-pi,pi,by=pi/2),labels=c("-pi","-pi/2","0","pi/2","pi"),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  stat_contour(mapping=aes(x=phase,y=frequency,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=phase,y=frequency,z=low),colour="orange",breaks=c(0),alpha=1)+
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
fbins<-seq(min(freq),max(freq),length=250)
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
fbins<-seq(min(freq),max(freq),length=250)
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
  scale_y_continuous(breaks=seq(20,100,by=10),expand=c(0,0))+
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
trsfit<-as.data.frame(ci[rtintsimidx,])
colnames(trsfit)<-c("low","fit","high")
trsfit$low<-trsfit$low+trsfit$fit
trsfit$high<-trsfit$high+trsfit$fit
trsfit$session<-rtintfit[,2]
trsfit$frequency<-rtintfit[,3]
trsfit$time<-rtintfit[,4]
trsfit$runningspeed<-rtintfit[,5]

#display plot
trsfit$frequency<-round(trsfit$frequency,digits=1)
trsfit$frequency<-as.ordered(trsfit$frequency)
for(f in 1:nf){
  levels(trsfit$frequency)[f]<-paste0(levels(trsfit$frequency)[f]," Hz")
}
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
trsfit<-as.data.frame(dci[rtintsimidx,])
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
cl<-max(abs(trsfit$fit))
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
cl<-max(abs(trsfit$fit))
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
cl<-max(abs(trsfit$fit))
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

#running speed theta phase interaction tiled by frequency
#run for lpmatrix
rspfit<-as.data.frame(ci[rstintsimidx,])
colnames(rspfit)<-c("low","fit","high")
rspfit$low<-rspfit$low+rspfit$fit
rspfit$high<-rspfit$high+rspfit$fit
rspfit$frequency<-rstintfit[,2]
rspfit$phase<-rstintfit[,3]
rspfit$runningspeed<-rstintfit[,4]

#display plot
rspfit$frequency<-round(rspfit$frequency,digits=1)
rspfit$frequency<-as.ordered(rspfit$frequency)
for(f in 1:nf){
  levels(trsfit$frequency)[f]<-paste0(levels(trsfit$frequency)[f]," Hz")
}
cl<-max(abs(rspfit$fit))
p<-ggplot(rspfit,aes(x=runningspeed,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,1,by=0.5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi/2),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=runningspeed,y=phase,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=runningspeed,y=phase,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=5) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))

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
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  # stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  # stat_contour(mapping=aes(x=phase,y=frequency,z=high),colour="purple",breaks=c(0),alpha=1)+
  # stat_contour(mapping=aes(x=phase,y=frequency,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~runningspeed)

#time theta phase interaction session1
#run for lpmatrix
tpfit<-as.data.frame(ci[tstintsimidx,])
colnames(tpfit)<-c("low","fit","high")
tpfit$low<-tpfit$low+tpfit$fit
tpfit$high<-tpfit$high+tpfit$fit
tpfit$session<-as.factor(tstintfit[,2])
tpfit$frequency<-tstintfit[,3]
tpfit$phase<-tstintfit[,4]
tpfit$time<-tstintfit[,5]

#display plot
tpfit$frequency<-round(tpfit$frequency,digits=1)
tpfit$frequency<-as.ordered(tpfit$frequency)
for(f in 1:nf){
  levels(trsfit$frequency)[f]<-paste0(levels(trsfit$frequency)[f]," Hz")
}
cl<-max(abs(tpfit$fit))
p<-ggplot(tpfit[tpfit$session=="1",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi/2),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=5) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))

#time theta phase interaction session2
p<-ggplot(tpfit[tpfit$session=="2",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi/2),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=5) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))

#time theta phase interaction session3
p<-ggplot(tpfit[tpfit$session=="3",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi/2),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=5) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))

#time/theta phase interactions by frequency session 1-2
tpfit<-as.data.frame(dci[tpintsimidx,])
colnames(tpfit)<-c("low","fit","high")
tpfit$low<-tpfit$low+tpfit$fit
tpfit$high<-tpfit$high+tpfit$fit
tpfit$session<-rep(c("1","2","3"),each=nb*nb*nf)
tpfit$frequency<-rep(rep(freq,each=nb*nb),times=3)
tpfit$time<-rep(rep(seq(mint,maxt,length.out=nb),each=nb),times=nf*3)
tpfit$phase<-rep(seq(min(pbins),max(pbins),length.out=nb),times=nb*nf*3)

tpfit$frequency<-round(tpfit$frequency,digits=1)
tpfit$frequency<-as.ordered(tpfit$frequency)
for(f in 1:nf){
  levels(tpfit$frequency)[f]<-paste0(levels(tpfit$frequency)[f]," Hz")
}
cl<-max(abs(tpfit$fit))
p<-ggplot(tpfit[tpfit$session=="1",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi/2),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=5) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))

#time/theta phase interactions by frequency session 1-3
p<-ggplot(tpfit[tpfit$session=="2",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi/2),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=5) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))

#time/theta phase interactions by frequency session 2-3
p<-ggplot(tpfit[tpfit$session=="3",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi/2),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="purple",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=5) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))

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
tpfit<-tpfit[tpfit$session=="3",]
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
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  # stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  # stat_contour(mapping=aes(x=phase,y=frequency,z=high),colour="purple",breaks=c(0),alpha=1)+
  # stat_contour(mapping=aes(x=phase,y=frequency,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~time)