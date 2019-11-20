#runGAMS_all_mice

jobid <- Sys.getenv('SLURM_JOB_ID')

#load packages
library(mgcv)
library(Rmpi)
library(snow)
library(rlecuyer)

#load data
setwd("/work/02592/bgereke/Data/GAMS_full")
MYDATA<-readRDS("MYDATAPOW_full.rds")

nboots <- 1000
nmice <- length(unique(MYDATA$Mouse))
ndays <- length(unique(MYDATA$Day))
bootidx <- matrix(0,nrow = 2*ndays,ncol = nboots)

#get bootstrap indices
for(b in 1:nboots){
  msamp <- sample(nmice,nmice,replace = TRUE)
  tmp <- 1
  for(m in 1:nmice){
    days_m <- unique(MYDATA$Day[MYDATA$Mouse == toString(msamp[m])])
    ndays_m <- length(days_m)
    bootidx[tmp:(tmp+ndays_m-1),b] <- sample(days_m,ndays_m,replace = TRUE)
    tmp <- tmp + ndays_m
  }
}

#set some vars
nb <- 60
mint<-0
maxt<-600
tbins<-seq(mint,maxt,length.out=nb)
minrs<-0
maxrs<-1
rsbins<-seq(minrs,maxrs,length.out=nb)
minp<- -pi
maxp<- pi
pbins<-seq(minp,maxp,length.out=nb)

#create cluster
cl<-getMPIcluster()

#load stuff onto each worker
clusterExport(cl,list("jobid","MYDATA","bootidx","nb","tbins","rsbins","pbins"))
ignore<-clusterEvalQ(cl,{library(mgcv);NULL})
clusterSetupRNG(cl,type="RNGstream")

#define parallel framework for task chunking with load balancing
parLapplyLB<-function(cl,x,fun,...)
{
  clusterCall(cl,LB.init,fun,...)
  r<-clusterApplyLB(cl,x,LB.worker)
  clusterEvalQ(cl,rm('.LB.fun','.LB.args',pos=globalenv()))
  r
}
LB.init<-function(fun,...)
{
  assign('.LB.fun',fun,pos=globalenv())
  assign('.LB.args',list(...),pos=globalenv())
  NULL
}
LB.worker<-function(x)
{
  do.call('.LB.fun',c(list(x),.LB.args))
}

#do parallel computation 
b<-seq(1,nboots)
parLapplyLB(cl,b,function(b) 
{
  #load packages
  library(mgcv)
  
  #remove unused days to save memory
  unused <- unique(MYDATA$Day)
  used <- unique(bootidx[,b])
  for(u in 1:length(used)){
    idx <- unused==toString(used[u])
    unused <- unused[!idx]
  }
  for(u in 1:length(unused)){
    MYDATA <- MYDATA[MYDATA$Day!=unused[u],]
  }
  
  #duplicate repeatedly used days
  for(u in 1:length(used)){
    n_d <- sum(bootidx[,b]==used[u])
    if(n_d>1){
      MYDATA <- rbind(MYDATA,do.call("rbind",replicate(n_d-1,MYDATA[MYDATA$Day==toString(used[u]),],simplify=FALSE)))
    }
  }
  
  #preallocate predictions
  nf = 50
  intercepts <- matrix(0,nrow=3,ncol=nf)
  runningspeed <- matrix(0,nrow=nb,ncol=nf)
  thetaphase <- matrix(0,nrow=nb,ncol=nf)
  speedphase <- matrix(0,nrow=nb^2,ncol=nf)
  timesession <- matrix(0,nrow=3*nb,ncol=nf)
  speedtimesession <- matrix(0,nrow=3*nb^2,ncol=nf)
  phasetimesession <- matrix(0,nrow=3*nb^2,ncol=nf)
  
  for(f in 1:nf){
    fmod <- bam(MYDATA[,f] ~ Session +
                  s(Running_Speed,bs="cr") +
                  s(Theta_Phase,bs="cc") + 
                  s(Time,bs="cr",by=Session)+
                  ti(Running_Speed,Theta_Phase,k=c(5,5),bs=c("cr","cc")) +
                  ti(Running_Speed,Time,k=c(3,3),bs=c("cr","cr"),by=Session) + 
                  ti(Time,Theta_Phase,k=c(5,5),bs=c("cr","cc"),by=Session)+
                  s(Running_Speed,bs="cr",by=Day)+
                  s(Theta_Phase,bs="cc",by=Day)+
                  s(Time,bs="cr",by=Day),
                method="fREML",
                family=gaussian(),
                data = MYDATA)
    
    #intercepts
    ints <- fmod$coefficients[1:3]
    ints[2:3] <- ints[2:3] + ints[1]
    intercepts[,f] <- ints
    
    #running speed
    newdat<-with(MYDATA,data.frame(Running_Speed=rep(rsbins,times=3),
                                   Theta_Phase=rep(pbins,times=3),
                                   Time=rep(tbins,times=3),
                                   Session=rep(c("1","2","3"),each=nb),
                                   Day=rep("1",times=3*nb),
                                   Trial=rep("1",times=3*nb)))
    pred<-predict.bam(fmod,newdat,type="terms",se.fit=FALSE,newdata.guaranteed = TRUE)
    runningspeed[,f] <- pred[1:nb,2]
    
    #theta phase
    thetaphase[,f] <- pred[1:nb,3]
    
    #time x session
    timesession[,f] <- rowSums(pred[,4:6])
    
    #running speed x theta phase
    dgrid<-expand.grid(rsbins,pbins)
    colnames(dgrid)<-c("runningspeed","thetaphase")
    newdat<-with(MYDATA,data.frame(Running_Speed=dgrid$runningspeed,
                                   Theta_Phase = dgrid$thetaphase,
                                   Time=rep(1,times=nrow(dgrid)),
                                   Session=rep("1",times=nrow(dgrid)),
                                   Day=rep("1",times=nrow(dgrid)),
                                   Trial=rep("1",times=nrow(dgrid))))
    pred<-predict.bam(fmod,newdat,type="terms",se.fit=FALSE,newdata.guaranteed = TRUE)
    speedphase[,f] <- pred[,7]
    
    #running speed x time x session
    dgrid<-expand.grid(rsbins,tbins,c("1","2","3"))
    colnames(dgrid)<-c("runningspeed","time","session")
    newdat<-with(MYDATA,data.frame(Running_Speed=dgrid$runningspeed,
                                   Theta_Phase = rep(1,times=nrow(dgrid)),
                                   Time=dgrid$time,
                                   Session=dgrid$session,
                                   Day=rep("1",times=nrow(dgrid)),
                                   Trial=rep("1",times=nrow(dgrid))))
    pred<-predict.bam(fmod,newdat,type="terms",se.fit=FALSE,newdata.guaranteed = TRUE)
    speedtimesession[,f] <- rowSums(pred[,8:10])
    
    #theta phase x time x session
    dgrid<-expand.grid(pbins,tbins,c("1","2","3"))
    colnames(dgrid)<-c("thetaphase","time","session")
    newdat<-with(MYDATA,data.frame(Running_Speed=rep(1,times=nrow(dgrid)),
                                   Theta_Phase = dgrid$thetaphase,
                                   Time=dgrid$time,
                                   Session=dgrid$session,
                                   Day=rep("1",times=nrow(dgrid)),
                                   Trial=rep("1",times=nrow(dgrid))))
    pred<-predict.bam(fmod,newdat,type="terms",se.fit=FALSE,newdata.guaranteed = TRUE)
    phasetimesession[,f] <- rowSums(pred[,11:13])
    
    print(paste("Completed: ",f," of ",nf),quote=FALSE)
  }
  #concatenate fits
  fits <- rbind(intercepts,runningspeed,thetaphase,speedphase,timesession,speedtimesession,phasetimesession)
  
  #save fits to disk
  filename<-paste0("boot_",b,".rds")
  saveRDS(fits,filename)
  NULL
})

stopCluster(cl)

