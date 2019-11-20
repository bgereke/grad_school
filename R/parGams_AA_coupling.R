#runGAMS_all_mice

jobid <- Sys.getenv('SLURM_JOB_ID')

#load packages
library(mgcv)
library(Rmpi)
library(snow)
library(rlecuyer)

#load data
setwd("/work/02592/bgereke/Data/GAMS_AA_coupling")
MYDATA<-readRDS("MYDATA_PA.rds")

#process days in serial
for(d in 1:length(unique(MYDATA$Day)))
{
  exday <- unique(MYDATA$Day)[d]
  MYDATA_d <- MYDATA[MYDATA$Day != exday,]
  #create cluster
  cl<-getMPIcluster()
  
  #load stuff onto each worker
  clusterExport(cl,list("MYDATA_d","jobid"))
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
  r<-seq(1,numsims)
  parLapplyLB(cl,r,function(r) 
  {
    
    #load packages
    library(mgcv)
    
    #Set model prediction vars
    nummmice<-length(unique(MYDATA$Mouse))
    numdays<-length(unique(MYDATA$Day))
    days<-seq(1,numdays)
    nb<-100
    nf<-50
    # mint<-0
    # maxt<-600
    # tbins<-seq(mint,maxt,length.out=nb)
    # minrs<-min(MYDATA$Running_Speed)
    # maxrs<-max(MYDATA$Running_Speed)
    # rsbins<-seq(minrs,maxrs,length.out=nb)
    # minp<- -pi
    # maxp<- pi
    # pbins<-seq(minp,maxp,length.out=nb)
    # mina<- -4
    # maxa<- 4
    # abins<-seq(mina,maxa,length.out=nb)
    # 
    # 
    # #create new data for predictions
    # dgrid<-expand.grid(abins,unique(MYDATA$Day))
    # colnames(dgrid)<-c("amp","day")
    # m<-matrix(length(unique(MYDATA$Day)),ncol=1)
    # for(d in 1:length(unique(MYDATA$Day)))
    # {
    #   idx<-which(MYDATA$Day == toString(d))
    #   m[d]<-MYDATA$Mouse[idx[1]]
    # }
    # newdat<-with(MYDATA,data.frame(Amp=dgrid$amp,Mouse=factor(rep(m,rep(nb,length(m)))),Day=dgrid$day))
    # 
    # #preallocate predictions
    # nv<-3
    # modpred<-matrix(nrow=length(dgrid[,1]),ncol=nv*nf)
    # fRsq<-matrix(nrow=1,ncol=nf)
    
    mb<-"cr"
    rb<-"fs"
    form<-as.formula(paste0("f",r," ~ s(f",f,",bs=mb)+s(f",f,",Mouse,bs=rb,xt=mb)+s(f",f,",Day,bs=rb,xt=mb)"))
    fmod <- bam(form,method="fREML",family=gaussian(),data = MYDATA)
    fRsq[f]<-summary(fmod)$r.sq
    colnames(newdat)[1]<-colnames(MYDATA)[f]
    modpred[,(nv*(f-1)+1):(nv*f)]<-predict(fmod, newdat, type = "terms")
    
    #save model predictions to disk
    setwd("/work/02592/bgereke/Data/GAMS_AA_coupling")
    filename<-paste0("modpred_",jobid,"_",r,".rds")
    saveRDS(modpred,filename)
    
    #save R-squared values to disk
    filename<-paste0("AARsq_",jobid,"_",r,".rds")
    saveRDS(fRsq,filename)
    
    NULL
  })
  
  stopCluster(cl)
  
}
