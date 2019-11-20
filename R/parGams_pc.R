#runGAMS_all_mice

# read MYDATA if not running previous code
setwd("/work/02592/bgereke/Data/GAMS_full_wint_rsnorm")
MYDATA<-readRDS("MYDATA.rds")

#normalize running speed
for(d in 1:max(as.numeric(MYDATA$Day)))
{
  MYDATA$Running_Speed[MYDATA$Day==toString(d)]<-MYDATA$Running_Speed[MYDATA$Day==toString(d)]/max(MYDATA$Running_Speed[MYDATA$Day==toString(d)])
}

#Set model prediciton vars
numdays<-length(unique(MYDATA$Day))
nb<-200
nf<-50
mint<-0
maxt<-600
day<-toString(1)
minrs<-min(MYDATA$Running_Speed[MYDATA$Day==day])
maxrs<-max(MYDATA$Running_Speed[MYDATA$Day==day])
rsbins<-seq(minrs,maxrs,length.out=nb)
tbins<-seq(0,600,length.out=nb)

#make new data to predict responses from
newdat<-with(MYDATA,data.frame(Running_Speed=c(rep(mean(Running_Speed),3*nb)),
                               Theta_Phase=c(rep(mean(Theta_Phase),3*nb)),
                               Time=c(rep(tbins,3)),
                               Session=c(rep("1",nb),rep("2",nb),rep("3",nb)),
                               Day=rep(day,3*nb)))
#add grids for interaction terms
inb = 200
irsbins<-seq(minrs,maxrs,length.out=inb)
ithbins<-seq(-pi,pi,length.out=inb)
itbins<-seq(0,600,length.out=inb)
trsgrid<-expand.grid(itbins,irsbins,unique(MYDATA$Session))
colnames(trsgrid)<-c("t","rs","ses")
tthgrid<-expand.grid(itbins,ithbins,unique(MYDATA$Session))
colnames(trsgrid)<-c("t","thp","ses")
numpts<-length(trsgrid$rs)
intdat<-with(MYDATA,data.frame(Running_Speed=c(trsgrid$rs,rep(mean(Running_Speed),numpts)),
                               Theta_Phase=c(rep(mean(Theta_Phase),numpts),tthgrid$thp),
                               Time=c(trsgrid$t,tthgrid$t),
                               Session=c(trsgrid$ses,tthgrid$ses),
                               Day=rep(day,2*numpts)))
newdat<-rbind(newdat,intdat)

#preallocate predictions (these are the matrices I want to save for each permutation)
sessiononeFR<-matrix(nrow=nf,ncol=nb)
sessiontwoFR<-matrix(nrow=nf,ncol=nb)
sessionthreeFR<-matrix(nrow=nf,ncol=nb)
timerunningspeedoneFR<-matrix(nrow=nf,ncol=numpts)
timerunningspeedtwoFR<-matrix(nrow=nf,ncol=numpts)
timerunningspeedthreeFR<-matrix(nrow=nf,ncol=numpts)
timethetaphaseoneFR<-matrix(nrow=nf,ncol=numpts)
timethetaphasetwoFR<-matrix(nrow=nf,ncol=numpts)
timethetaphasethreeFR<-matrix(nrow=nf,ncol=numpts)

#load snow libraries 
library(Rmpi)
library(snow)

cl<-getMPIcluster()

parLapply(cl,1:10000,function(i)
{
  #Shuffle session id's (should be the same for each f, but unique for each i)
  tempSess <- MYDATA$Session
  for(d in 1:max(as.numeric(MYDATA$Day)))
  {
    samp <- sample(3) #make sure this is truly unique across nodes 
    tempSess[as.numeric(MYDATA$Day) == d & as.numeric(MYDATA$Session) == 1] <- rep(toString(samp[1]),sum(as.numeric(MYDATA$Day) == d & as.numeric(MYDATA$Session) == 1))
    tempSess[as.numeric(MYDATA$Day) == d & as.numeric(MYDATA$Session) == 2] <- rep(toString(samp[2]),sum(as.numeric(MYDATA$Day) == d & as.numeric(MYDATA$Session) == 2))
    tempSess[as.numeric(MYDATA$Day) == d & as.numeric(MYDATA$Session) == 3] <- rep(toString(samp[3]),sum(as.numeric(MYDATA$Day) == d & as.numeric(MYDATA$Session) == 3))
  }
  MYDATA$Session <- tempSess
  rm(tempSess)
  
  #run gam on all frequencies
  #load other libraries (Is this fine or do these need to be exported to cluster nodes?)
  library(mgcv)
  library(parallel) 
  library(MASS)
  
  mb<-"cr"
  #for(f in 1:nf)
  #{
  mclapply(1:50,function(f)
  {
    fmod <- bam(MYDATA[,f] ~ s(Running_Speed,bs=mb,k=8,fx=TRUE,by=Day)+s(Theta_Phase,bs="cc",k=6,fx=TRUE,by=Day)+s(Time,bs=mb,k=9,fx=TRUE,by=Session)
                +ti(Running_Speed,Theta_Phase,bs=c(mb,"cc"),fx=TRUE,k=c(4,4),by=Day)+ti(Time,Running_Speed,bs=c(mb,mb),k=c(6,3),fx=TRUE,by=Session)
                +ti(Time,Theta_Phase,bs=c(mb,"cc"),k=c(6,3),fx=TRUE,by=Session),method="fREML",family=gaussian(),data = MYDATA)
    
    # assign(names(MYDATA)[f],fmod)
    # filename<-paste(names(MYDATA)[f],".rds",sep="")
    # saveRDS(fmod,filename)
    
    #lpmatrix for simultaneous confidence intervals
    newpred<-predict(fmod, newdat, type = "lpmatrix")
    coefs<-coef(fmod)
    vc<-vcov(fmod)
    sim<-mvrnorm(numsims,mu=coefs,Sigma=vc)
    
    #get Session One fits
    want<-grep("Time",colnames(newpred))
    dwant<-diff(want)
    didx<-which(dwant>1)
    want<-want[1:didx[1]]
    lwant<-length(want)
    numcols<-lwant/3
    want1<-want[1:numcols]
    fits<-newpred[,want1] %*% t(sim[,want1])
    sessiononeFR[f,]<-apply(fits[(1):(nb),],1,quantile,probs=c(0.5))
    #get Session Two fits
    want2<-want[(numcols+1):(2*numcols)]
    fits<-newpred[,want2] %*% t(sim[,want2])
    sessiontwoFR[f,]<-apply(fits[(nb+1):(2*nb),],1,quantile,probs=c(0.5))
    #get Session Three fits
    want3<-want[(2*numcols+1):(3*numcols)]
    fits<-newpred[,want3] %*% t(sim[,want3])
    sessionthreeFR[f,]<-apply(fits[(2*nb+1):(3*nb),],1,quantile,probs=c(0.5))
    
    #get interactions
    
    #get Time x Running_Speed Session 1
    want<-grep("Time,Running_Speed",colnames(newpred))
    lwant<-length(want)
    numcols<-lwant/3
    want1<-want[1:numcols]
    fits<-newpred[,want1] %*% t(sim[,want1])
    timerunningspeedoneFR[f,]<-apply(fits[(3*nb+1+numpts):(3*nb+2*numpts),],1,quantile,probs=c(0.5))
    #get Time x Running_Speed Session 2
    want2<-want[(numcols+1):(2*numcols)]
    fits<-newpred[,want2] %*% t(sim[,want2])
    timerunningspeedtwoFR[f,]<-apply(fits[(3*nb+1+2*numpts):(3*nb+3*numpts),],1,quantile,probs=c(0.5))
    #get Time x Running_Speed Session 3
    want3<-want[(2*numcols+1):(3*numcols)]
    fits<-newpred[,want3] %*% t(sim[,want3])
    timerunningspeedthreeFR[f,]<-apply(fits[(3*nb+1+3*numpts):(3*nb+4*numpts),],1,quantile,probs=c(0.5))
    #get Time x Theta_Phase Session 1
    want<-grep("Time,Theta_Phase",colnames(newpred))
    lwant<-length(want)
    numcols<-lwant/3
    want1<-want[1:numcols]
    fits<-newpred[,want1] %*% t(sim[,want1])
    timerunningspeedoneFR[f,]<-apply(fits[(3*nb+1+4*numpts):(3*nb+5*numpts),],1,quantile,probs=c(0.5))
    #get Time x Theta_Phase Session 2
    want2<-want[(numcols+1):(2*numcols)]
    fits<-newpred[,want2] %*% t(sim[,want2])
    timerunningspeedtwoFR[f,]<-apply(fits[(3*nb+1+5*numpts):(3*nb+6*numpts),],1,quantile,probs=c(0.5))
    #get Time x Theta_Phase Session 3
    want3<-want[(2*numcols+1):(3*numcols)]
    fits<-newpred[,want3] %*% t(sim[,want3])
    timerunningspeedthreeFR[f,]<-apply(fits[(3*nb+1+6*numpts):(3*nb+7*numpts),],1,quantile,probs=c(0.5))
    
    #print(paste("Completed: ",2*f," Hz"),quote=FALSE)
  },mc.cores=10)
  
  
  #save FR matrices as sessiononeFR_i.rds or something like that into some folder
  filename<-paste("sessiononeFR_",i,".rds",sep="")
  saveRDS(sessiononeFR,filename)
  filename<-paste("sessiontwoFR_",i,".rds",sep="")
  saveRDS(sessiontwoFR,filename)
  filename<-paste("sessionthreeFR_",i,".rds",sep="")
  saveRDS(sessionthreeFR,filename)
  
  filename<-paste("timerunningspeedoneFR_",i,".rds",sep="")
  saveRDS(timerunningspeedoneFR,filename)
  filename<-paste("timerunningspeedtwoFR_",i,".rds",sep="")
  saveRDS(stimerunningspeedtwoFR,filename)
  filename<-paste("timerunningspeedtwoFR_",i,".rds",sep="")
  saveRDS(timerunningspeedtwoFR,filename)
  
})

