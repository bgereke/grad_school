#runGAMS_all_mice

#load packages
library(mgcv)
library(stats)
library(parallel)
library(itsadug)
library(pracma)
library(fdasrvf)

#load data
setwd("E:/For Brian/GAMS_full/")
MYDATA<-readRDS("E:/For Brian/GAMS_full/MYDATAPOW_w.rds")
MYDATA$cPhase <- cos(MYDATA$Theta_Phase)
MYDATA$sPhase <- sin(MYDATA$Theta_Phase)
# MYDATA$Strain <- ordered(MYDATA$Strain,levels=c("hyb","3xTG"))
# contrasts(MYDATA$Strain) <- 'contr.treatment'

# setwd("G:/Rats/GAMS_full_CA1/")
# MYDATA<-readRDS("MYDATAPOW_full_CA1.rds")

# setwd("G:/Rats/GAMS_full_CCS/")
# MYDATA<-readRDS("MYDATA_full_CCS.rds")

# MYDATA <- MYDATA[MYDATA$Running_Speed>5,]
# MYDATA$ARStart<-rep(FALSE,length(MYDATA[,1]))
# MYDATA$ARStart[abs(tdiff)>1]<-TRUE
# MYDATA$ARStart[1]<-TRUE
nf<-50
freq<-10^seq(log10(2),log10(100),length.out=nf)
width <- 7
s <- (width+sqrt(2+width^2))/(4*pi*freq)
dt <- median(diff(MYDATA$Time))
c <- Mod(exp(complex(real=0,imaginary=1)*width/s*dt)*exp(-0.5*dt^2/(2*s^2)))


# MYDATA <- MYDATA[MYDATA$Strain == "hyb",]

# #register spectra
# specs <- matrix(0,nrow=nf,ncol=length(unique(MYDATA$Day)))
# for(d in 1:length(unique(MYDATA$Day))){
#   specs[,d] <- colMeans(MYDATA[MYDATA$Day==toString(d),1:nf])
# }
# warps <- multiple_align_functions(specs,log2(freq),rowMeans(specs),lambda=0)
# MYDATA_w <- MYDATA
# for(d in 1:length(unique(MYDATA$Day))){
#   MYDATA_d <- MYDATA[MYDATA$Day==toString(d),1:nf]
#   for(i in 1:nrow(MYDATA_d)){
#     MYDATA_d[i,] <- warp_f_gamma(MYDATA_d[i,],log2(freq),warps$gam[,d])
#   }
#   MYDATA_w[MYDATA$Day==toString(d),1:nf] <- MYDATA_d
#   print(paste("Completed: ",d," of ",length(unique(MYDATA$Day))),quote=FALSE)
# }
# saveRDS(MYDATA_w,file="MYDATAPOW_w.rds")
# 
# rm(MYDATA,MYDATA_d,specs,warps)
# gc()

# #setup cluster
# no_cores <- 1 #detectCores() #- 1
# cl<-makeCluster(no_cores, type='PSOCK')
# clusterExport(cl,list("nf","MYDATA"))
# clusterEvalQ(cl,{library(mgcv);library(stats);NULL})

#do parallel computation
# parLapply(cl,21:nf,function(f)
for(f in 1:nf)
{
  
  if(f > 1){
    
    fmod <- bam(MYDATA[,f] ~ StrainSession + 
                s(Running_Speed,by=Strain,bs="cr") +
                s(Theta_Phase,by=Strain,bs="cc") +
                s(Time,bs="cr",by=StrainSession) +
                s(Theta_Phase,Day,bs="fs",xt="cc",k=5) +
                s(Running_Speed,Day,bs="fs",xt="cr",k=5) +
                s(Time,Trial,bs="fs",xt="cr",k=5),
                coef = coefs,
                method="fREML",
                data = MYDATA,
                family = Gamma(link="log"),
                discrete=TRUE,
                nthreads=4)
    
    coefs <- fmod$coefficients
    
  } else {
      
    fmod <- bam(MYDATA[,f] ~ StrainSession + 
                s(Running_Speed,by=Strain,bs="cr") +
                s(Theta_Phase,by=Strain,bs="cc") +
                s(Time,bs="cr",by=StrainSession) +
                s(Theta_Phase,Day,bs="fs",xt="cc",k=5) +
                s(Running_Speed,Day,bs="fs",xt="cr",k=5) +
                s(Time,Trial,bs="fs",xt="cr",k=5),
                method="fREML",
                data = MYDATA,
                family = Gamma(link="log"),
                discrete=TRUE,
                nthreads=4)
    
    coefs <- fmod$coefficients
    
    }
  
  # r <- 0.9*start_value_rho(fmod)
  # 
  # fmod <- bam(MYDATA[,f] ~ StrainSession + Day +
  #             # s(Running_Speed,by=Strain,bs="cr") +
  #             # s(Theta_Phase,by=Strain,bs="cc") +
  #             # s(Time,bs="cr",by=StrainSession) +
  #             s(Mouse,bs="re"),
  #             coef = coefs,
  #             method="fREML",
  #             AR.start = ARStart,
  #             rho = r,
  #             data = MYDATA,
  #             family = Gamma(link="log"),
  #             discrete=TRUE,
  #             nthreads=4)

  #save model to disk
  filename<-paste0("warped_test_f",f,".rds")
  saveRDS(fmod,filename)
  
  print(paste("Completed: ",f," of ",nf),quote=FALSE)
  
  NULL
}#)
#stopCluster(cl)
