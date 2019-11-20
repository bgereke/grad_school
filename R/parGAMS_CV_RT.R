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
setwd("/work/02592/bgereke/Data/GAMS_CV")

nf<-50
# r<-matrix(0,nrow = nf,ncol = 1)
# for(f in 1:nf){r[f]<-readRDS(paste0("rhos_RT_f",f,".rds"))}


#process days in serial
for(m in 1:6)
{
  #leave one mouse out
  exmouse <- unique(MYDATA$Mouse)[m]
  MYDATA_m <- MYDATA[MYDATA$Mouse != exmouse,]
  
  #create cluster
  cl<-getMPIcluster()
  
  #load stuff onto each worker
  clusterExport(cl,list("MYDATA_m","m"))
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
  numsims<-50
  f<-seq(1,numsims)
  parLapplyLB(cl,f,function(f) 
  {
    
    #load packages
    library(mgcv)
    
    fmod <- bam(MYDATA_m[,f] ~ Session +
			       s(Running_Speed,bs="cr") +
			       s(Theta_Phase,bs="cc") + 
			       s(Time,bs="cr",by=Session)+
			       ti(Running_Speed,Theta_Phase,k=c(5,5),bs=c("cr","cc")) + 
			       ti(Running_Speed,Time,k=c(3,3),bs=c("cr","cr"),by=Session),
             method="fREML",
			       family=gaussian(),
			       data = MYDATA_m)
    
    #save model to disk
    setwd("/work/02592/bgereke/Data/GAMS_CV")
    filename<-paste0("Base_TS_RT_","f",f,"_m",m,"_POW.rds")
    saveRDS(fmod,filename)
    
    NULL
  })
  
}

stopCluster(cl)