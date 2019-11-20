#runGAMS_all_mice

jobid <- Sys.getenv('SLURM_JOB_ID')

#load packages
library(mgcv)
library(Rmpi)
library(snow)
library(rlecuyer)

#load data
setwd("/work/02592/bgereke/Data/GAMS_full_wint_rsnorm")
MYDATA<-readRDS("MYDATA.rds")
sess<-unique(MYDATA$Session)

#create cluster
cl<-getMPIcluster()

#load stuff onto each worker
clusterExport(cl,list("MYDATA","sess","jobid"))
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
numsims<-24
i<-seq(1,numsims)
parLapplyLB(cl,i,function(i) 
{
  
  #load packages
  library(mgcv)
  
  #Shuffle session id's (should be the same for each f, but unique for each i)
  tempSess <- MYDATA$Session
  for(d in 1:max(as.numeric(MYDATA$Day)))
  {
    samp <- sample(3) #make sure this is truly unique across nodes 
    tempSess[as.numeric(MYDATA$Day) == d & as.numeric(MYDATA$Session) == 1] <- factor(rep(toString(samp[1]),sum(as.numeric(MYDATA$Day) == d & as.numeric(MYDATA$Session) == 1)))
    tempSess[as.numeric(MYDATA$Day) == d & as.numeric(MYDATA$Session) == 2] <- factor(rep(toString(samp[2]),sum(as.numeric(MYDATA$Day) == d & as.numeric(MYDATA$Session) == 2)))
    tempSess[as.numeric(MYDATA$Day) == d & as.numeric(MYDATA$Session) == 3] <- factor(rep(toString(samp[3]),sum(as.numeric(MYDATA$Day) == d & as.numeric(MYDATA$Session) == 3)))
  }
  MYDATA$Session <- tempSess
  rm(tempSess)
  
  #Set model prediction vars
  numdays<-length(unique(MYDATA$Day))
  days<-seq(1,numdays)
  nb<-100
  nf<-50
  mint<-0
  maxt<-600
  tbins<-seq(mint,maxt,length.out=nb)
  minrs<-min(MYDATA$Running_Speed)
  maxrs<-max(MYDATA$Running_Speed)
  rsbins<-seq(minrs,maxrs,length.out=nb)
  minp<- -pi
  maxp<- pi
  pbins<-seq(minp,maxp,length.out=nb)
  
  
  #create new data for predictions
  dgrid<-expand.grid(tbins,sess)
  colnames(dgrid)<-c("time","session")
  newdat<-with(MYDATA,data.frame(Time=dgrid$time,Session=dgrid$session))
  
  #preallocate predictions
  nv<-3
  modpred<-matrix(nrow=length(newdat[,1]),ncol=nv*nf)
  fRsq<-matrix(nrow=1,ncol=nf) 
  
  mb<-"cr"
  for(f in 1:nf)
  {
    "before"
    fmod <- bam(MYDATA[,f] ~ s(Time,bs=mb,by=Session)
                ,method="fREML",family=gaussian(),data = MYDATA,nthreads=1)
    "after"
    fRsq[f]<-summary(fmod)$r.sq
    modpred[,(nv*(f-1)+1):(nv*f)]<-predict(fmod, newdat, type = "terms")
  }
  
  #save model predictions to disk
  setwd("/work/02592/bgereke/Data/GAMS_full_wint_rsnorm")
  filename<-paste0("modpred_",jobid,"_",i,".rds")
  saveRDS(modpred,filename)

  NULL
})

stopCluster(cl)
