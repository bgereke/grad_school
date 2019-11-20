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
MYDATA<-MYDATA[c(101:150,(ncol(MYDATA)-10):ncol(MYDATA))] #CICS only
#MYDATA<-MYDATA[-c(51:(ncol(MYDATA)-11))] #PICS only
#MYDATA<-MYDATA[c(201:250,(ncol(MYDATA)-10):ncol(MYDATA))] #PP only

#process days in serial
for(d in 1:length(unique(MYDATA$Day)))
{
  #leave one day out
  exday <- unique(MYDATA$Day)[d]
  MYDATA_d <- MYDATA[MYDATA$Day != exday,]
  
  #create cluster
  cl<-getMPIcluster()
  
  #load stuff onto each worker
  clusterExport(cl,list("MYDATA_d","jobid","d"))
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
    fmod <- bam(MYDATA_d[,f] ~ s(Running_Speed,bs=mb) + s(Time,bs=mb,by=Session) + s(Running_Speed,Trial,bs=rb,xt=mb,m=1) + s(Time,Trial,bs=rb,xt=mb,m=1)

                ,method="fREML",family=gaussian(),data = MYDATA_d)
    
    #save model to disk
    setwd("/work/02592/bgereke/Data/GAMS_CV")
    filename<-paste0("VTbS_","f",f,"_d",d,"_CICS.rds")
    saveRDS(fmod,filename)
    
    NULL
  })
  
}

stopCluster(cl)
