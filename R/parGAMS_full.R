#runGAMS_all_mice

jobid <- Sys.getenv('SLURM_JOB_ID')

#load packages
library(mgcv)
library(Rmpi)
library(snow)
library(rlecuyer)

#load data
setwd("/work/02592/bgereke/Data/GAMS_CV")
MYDATA<-readRDS("MYDATA_full.rds")
MYDATA<-MYDATA[-c(51:(ncol(MYDATA)-11))] #PICS only

#create cluster
cl<-getMPIcluster()

#load stuff onto each worker
clusterExport(cl,list("MYDATA","jobid"))
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
  
  mb<-"cr"
  rb<-"fs"
  fmod <- bam(MYDATA[,f] ~ s(Running_Speed,bs=mb) + s(Running_Speed,Day,bs=rb,xt=mb,m=1)
              ,method="fREML",family=gaussian(),data = MYDATA)
  
  #save model to disk
  setwd("/work/02592/bgereke/Data/GAMS_CV")
  filename<-paste0("V_","f",f,"_full_PICS.rds")
  saveRDS(fmod,filename)
  
  NULL
})

stopCluster(cl)
