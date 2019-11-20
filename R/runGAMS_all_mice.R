#runGAMS_all_mice

#load my functions
setwd("C:/Users/Brian/Documents/R")

#load libraries  
library(RevoUtilsMath)
library(mgcv)
library(ks)
library(signal)
library(data.table)
library(circular)

# if(require("RevoUtilsMath")){
#   setMKLthreads(2)
#   print("threads set")
# }

#go to data directory
setwd("C:/Data")
miceDirs<-readLines("mice.txt",warn=FALSE,skipNul=TRUE)
numMice<-length(miceDirs)
day<-1
trial<-1
rsth<-5 #running speed threshold
tth<-600 #time threshold
nf<-60 #number of frequencies

for(m in 1:numMice)
  {
  print(c("Reading data for: ",miceDirs[m]),quote=FALSE)
  setwd(miceDirs[m])
  
  # read the dates.txt file
  dates<-readLines("dates.txt",warn=FALSE,skipNul=TRUE)
  numDates<-length(dates)
  
  for(d in 1:numDates)
    {
    print(dates[d])
    setwd(paste(miceDirs[m],"\\",dates[d],sep=""))
    infile<-"infile - Copy.txt"
    #infile<-"infile.txt"
    # infile<-"InFile.txt"
    
    #load and merge data from sessions
    inlines<-readLines(infile,warn=FALSE,skipNul=TRUE)
    
    #session 1 data
    #position data
    s1data<-read.table(paste(inlines[5],"\\CS\\rtable_tracking.csv",sep=""),sep=",",header=TRUE)
    s1data$Session <- factor(rep(toString(1),nrow(s1data)))
    s1data$Trial <- factor(rep(toString(trial),nrow(s1data)))
    s1data<-s1data[s1data$Time<=tth,]
    # #spiking data
    # sp1data<-read.table(paste(inlines[5],"\\CS\\rtable_TS.csv",sep=""),sep=",",header=TRUE)[,1:2]
    # sp1data$CellID<-as.factor(sp1data$CellID)
    # sp1data$Session <- factor(rep(toString(1),nrow(sp1data)))
    # sp1data$Trial <- factor(rep(toString(trial),nrow(sp1data)))
    # sp1data<-sp1data[sp1data$SpTime<=tth,]
    # #spectral data
    # csc1data<-fread(paste(inlines[5],"\\CS\\rtable_CSC.csv",sep=""),sep=",",header=TRUE,drop=c(1:100))
    # csc1data<-as.data.frame(csc1data)
    # csc1data<-csc1data[csc1data$Time>0&csc1data$Time<=tth,]
    # #get spike slow/fast theta phases
    # spkphase1<-matrix(0,nrow=nrow(sp1data),ncol=2)
    # for(i in 1:nrow(spkphase1)){
    #   ttidx<-which.min((sp1data$SpTime[i]-csc1data$Time)^2)
    #   spkphase1[i,1]<-Arg(complex(real=csc1data$RephasePP_f14[ttidx],imaginary=csc1data$ImphasePP_f14[ttidx]))
    #   spkphase1[i,2]<-Arg(complex(real=csc1data$RephasePP_f20[ttidx],imaginary=csc1data$ImphasePP_f20[ttidx]))
    # }
    # #downsample spectral data
    # ttidx<-rep(0,nrow(s1data))
    # for(i in 1:nrow(s1data)){ttidx[i]<-which.min((csc1data$Time-s1data$Time[i])^2)}
    # csc1data<-csc1data[ttidx,]
    gc()
    
    trial<-trial+1
    
    #session 2 data
    #position data
    s2data<-read.table(paste(inlines[6],"\\CS\\rtable_tracking.csv",sep=""),sep=",",header=TRUE)
    s2data$Session <- factor(rep(toString(2),nrow(s2data)))
    s2data$Trial <- factor(rep(toString(trial),nrow(s2data)))
    s2data<-s2data[s2data$Time<=tth,]
    # #spiking data
    # sp2data<-read.table(paste(inlines[6],"\\CS\\rtable_TS.csv",sep=""),sep=",",header=TRUE)[,1:2]
    # sp2data$CellID<-as.factor(sp2data$CellID)
    # sp2data$Session <- factor(rep(toString(2),nrow(sp2data)))
    # sp2data$Trial <- factor(rep(toString(trial),nrow(sp2data)))
    # sp2data<-sp2data[sp2data$SpTime<=tth,]
    # #spectral data
    # csc2data<-fread(paste(inlines[6],"\\CS\\rtable_CSC.csv",sep=""),sep=",",header=TRUE,drop=c(1:100))
    # csc2data<-as.data.frame(csc2data)
    # csc2data<-csc2data[csc2data$Time>0&csc2data$Time<=tth,]
    # #get spike slow/fast theta phases
    # spkphase2<-matrix(0,nrow=nrow(sp2data),ncol=2)
    # for(i in 1:nrow(spkphase2)){
    #   ttidx<-which.min((sp2data$SpTime[i]-csc2data$Time)^2)
    #   spkphase2[i,1]<-Arg(complex(real=csc2data$RephasePP_f14[ttidx],imaginary=csc2data$ImphasePP_f14[ttidx]))
    #   spkphase2[i,2]<-Arg(complex(real=csc2data$RephasePP_f20[ttidx],imaginary=csc2data$ImphasePP_f20[ttidx]))
    # }
    # #downsample spectral data
    # ttidx<-rep(0,nrow(s2data))
    # for(i in 1:nrow(s2data)){ttidx[i]<-which.min((csc2data$Time-s2data$Time[i])^2)}
    # csc2data<-csc2data[ttidx,]
    gc()
    
    trial<-trial+1
    
    #session 3 data
    #position data
    s3data<-read.table(paste(inlines[7],"\\CS\\rtable_tracking.csv",sep=""),sep=",",header=TRUE)
    s3data$Session <- factor(rep(toString(3),nrow(s3data)))
    s3data$Trial <- factor(rep(toString(trial),nrow(s3data)))
    s3data<-s3data[s3data$Time<=tth,]
    # #spiking data
    # sp3data<-read.table(paste(inlines[7],"\\CS\\rtable_TS.csv",sep=""),sep=",",header=TRUE)[,1:2]
    # sp3data$CellID<-as.factor(sp3data$CellID)
    # sp3data$Session <- factor(rep(toString(3),nrow(sp3data)))
    # sp3data$Trial <- factor(rep(toString(trial),nrow(sp3data)))
    # sp3data<-sp3data[sp3data$SpTime<=tth,]
    # #spectral data
    # csc3data<-fread(paste(inlines[7],"\\CS\\rtable_CSC.csv",sep=""),sep=",",header=TRUE,drop=c(1:100))
    # csc3data<-as.data.frame(csc3data)
    # csc3data<-csc3data[csc3data$Time>0&csc3data$Time<=tth,]
    # #get spike slow/fast theta phases
    # spkphase3<-matrix(0,nrow=nrow(sp3data),ncol=2)
    # for(i in 1:nrow(spkphase3)){
    #   ttidx<-which.min((sp3data$SpTime[i]-csc3data$Time)^2)
    #   spkphase3[i,1]<-Arg(complex(real=csc3data$RephasePP_f14[ttidx],imaginary=csc3data$ImphasePP_f14[ttidx]))
    #   spkphase3[i,2]<-Arg(complex(real=csc3data$RephasePP_f20[ttidx],imaginary=csc3data$ImphasePP_f20[ttidx]))
    # }
    # #downsample spectral data
    # ttidx<-rep(0,nrow(s3data))
    # for(i in 1:nrow(s3data)){ttidx[i]<-which.min((csc3data$Time-s3data$Time[i])^2)}
    # csc3data<-csc3data[ttidx,]
    gc()
    
    trial<-trial+1
    
    # #rotate theta phases to max place cell firing
    # #session 1
    # csc1data$STphase<-Arg(complex(real=csc1data$RephasePP_f14,imaginary=csc1data$ImphasePP_f14))
    # csc1data$FTphase<-Arg(complex(real=csc1data$RephasePP_f20,imaginary=csc1data$ImphasePP_f20))
    # phasegrid<-seq(-pi,pi,length.out=512)
    # d1<-density.circular(spkphase1[,1],bw=10,from=circular(-pi),to=circular(pi))
    # d2<-density.circular(spkphase1[,2],bw=10,from=circular(-pi),to=circular(pi))
    # d1max<-phasegrid[which.max(d1$y)]
    # d2max<-phasegrid[which.max(d2$y)]
    # csc1data$STphase<-(unwrap(csc1data$STphase)-d1max-pi) %% (2*pi) - pi
    # csc1data$FTphase<-(unwrap(csc1data$FTphase)-d2max-pi) %% (2*pi) - pi
    # #session 2
    # csc2data$STphase<-Arg(complex(real=csc2data$RephasePP_f14,imaginary=csc2data$ImphasePP_f14))
    # csc2data$FTphase<-Arg(complex(real=csc2data$RephasePP_f20,imaginary=csc2data$ImphasePP_f20))
    # phasegrid<-seq(-pi,pi,length.out=512)
    # d1<-density.circular(spkphase2[,1],bw=10,from=circular(-pi),to=circular(pi))
    # d2<-density.circular(spkphase2[,2],bw=10,from=circular(-pi),to=circular(pi))
    # d1max<-phasegrid[which.max(d1$y)]
    # d2max<-phasegrid[which.max(d2$y)]
    # csc2data$STphase<-(unwrap(csc2data$STphase)-d1max-pi) %% (2*pi) - pi
    # csc2data$FTphase<-(unwrap(csc2data$FTphase)-d2max-pi) %% (2*pi) - pi
    # #session3
    # csc3data$STphase<-Arg(complex(real=csc3data$RephasePP_f14,imaginary=csc3data$ImphasePP_f14))
    # csc3data$FTphase<-Arg(complex(real=csc3data$RephasePP_f20,imaginary=csc3data$ImphasePP_f20))
    # phasegrid<-seq(-pi,pi,length.out=512)
    # d1<-density.circular(spkphase3[,1],bw=10,from=circular(-pi),to=circular(pi))
    # d2<-density.circular(spkphase3[,2],bw=10,from=circular(-pi),to=circular(pi))
    # d1max<-phasegrid[which.max(d1$y)]
    # d2max<-phasegrid[which.max(d2$y)]
    # csc3data$STphase<-(unwrap(csc3data$STphase)-d1max-pi) %% (2*pi) - pi
    # csc3data$FTphase<-(unwrap(csc3data$FTphase)-d2max-pi) %% (2*pi) - pi
    
    #smooth position samples and compute running speed by finite differences
    s1data$Position <- unwrap(s1data$Position)
    s2data$Position <- unwrap(s2data$Position)
    s3data$Position <- unwrap(s3data$Position)
    smPos1 <- bam(s1data$Position ~ s(Time,bs="cr",k=200,fx=TRUE),data = s1data)
    smPos2 <- bam(s2data$Position ~ s(Time,bs="cr",k=200,fx=TRUE),data = s2data)
    smPos3 <- bam(s3data$Position ~ s(Time,bs="cr",k=200,fx=TRUE),data = s3data)
    eps <- 0.0001
    newd <- data.frame(Time=s1data$Time)
    s1data$Position <- predict.bam(smPos1,newd)
    newd <- data.frame(Time=s1data$Time+eps)
    predeps <- predict.bam(smPos1,newd)
    s1data$Running_Speed <- abs(predeps - s1data$Position)/eps
    newd <- data.frame(Time=s2data$Time)
    s2data$Position <- predict.bam(smPos2,newd)
    newd <- data.frame(Time=s2data$Time+eps)
    s2data$Running_Speed <- abs(predict.bam(smPos2,newd) - s2data$Position) / eps
    newd <- data.frame(Time=s3data$Time)
    s3data$Position <- predict.bam(smPos3,newd)
    newd <- data.frame(Time=s3data$Time+eps)
    s3data$Running_Speed <- abs(predict.bam(smPos3,newd) - s3data$Position) / eps
    s1data$Position <- (s1data$Position-pi) %% (2*pi) - pi
    s2data$Position <- (s2data$Position-pi) %% (2*pi) - pi
    s3data$Position <- (s3data$Position-pi) %% (2*pi) - pi
    
    # s1data<-cbind(csc1data[,1:(2*nf)],s1data)
    # s2data<-cbind(csc2data[,1:(2*nf)],s2data)
    # s3data<-cbind(csc3data[,1:(2*nf)],s3data)
    mydata <- rbind(s1data,s2data,s3data)
    # spkphase<-rbind(spkphase1,spkphase2,spkphase3)
    # rm(s1data,s2data,s3data,csc1data,csc2data,csc3data)
    gc()
    mydata$Day<-factor(rep(toString(day),length(mydata[,1])))
    mydata$Mouse<-factor(rep(toString(m),length(mydata[,1])))
    if (m<=3)
    {
      mydata$Strain<-factor(rep("hyb",length(mydata[,1])))
    } else{
      mydata$Strain<-factor(rep("B6",length(mydata[,1])))
    }
    day<-day+1
    # mydata <- mydata[mydata$Running_Speed >= rsth,]  #set running speed threshold
    mydata <- mydata[mydata$Time <= tth,]  #set time threshold 
    # mydata$Running_Speed<-scale(mydata$Running_Speed,center=TRUE,scale=TRUE)
    
    #remove data close to reward location
    vmap <- bam(mydata$Running_Speed ~ s(Position,bs="cc",k=30,fx=TRUE),data = mydata)
    pmap <- seq(-pi,pi,length.out =1000)
    newd <- data.frame(Position=pmap)
    vpred <- predict.bam(vmap,newd)
    rewloc <- pmap[which(vpred == min(vpred))]
    #rewloc<-read.table(paste(getwd(),"\\rewloc\\rewloc.csv",sep=""),sep=",",header=TRUE)
    delta<-abs(Arg(complex(argument=rewloc)*Conj(complex(argument=mydata$Position))))/(2*pi)
    D<-96.5 #track diameter
    rd<-15 #distance from reward loc
    mydata<-mydata[delta>=rd/(pi*D),]
    
    # # #rotate theta phases to max place cell firing
    # # mydata$STphase<-Arg(complex(real=mydata$RephasePP_f14,imaginary=mydata$ImphasePP_f14))
    # # mydata$FTphase<-Arg(complex(real=mydata$RephasePP_f20,imaginary=mydata$ImphasePP_f20))
    # # phasegrid<-seq(-pi,pi,length.out=512)
    # # d1<-density.circular(spkphase[,1],bw=10,from=circular(-pi),to=circular(pi))
    # # d2<-density.circular(spkphase[,2],bw=10,from=circular(-pi),to=circular(pi))
    # # d1max<-phasegrid[which.max(d1$y)]
    # # d2max<-phasegrid[which.max(d2$y)]
    # # mydata$STphase<-(unwrap(mydata$STphase)-d1max-pi) %% (2*pi) - pi
    # # mydata$FTphase<-(unwrap(mydata$FTphase)-d2max-pi) %% (2*pi) - pi
    # 
    # #get PC1 data only
    # mydata[,1:nf]<-matrix(log(Mod(complex(real=as.matrix(mydata[,1:nf]),imaginary=as.matrix(mydata[,(nf+1):(2*nf)])))^2),ncol=nf)
    # mydata<-mydata[,-((nf+1):(2*nf))]
    # #transform response variables
    # for(f in 1:nf)
    # {
    #   # mydata[,f]<-log(mydata[,f])
    #   # mydata[,f+4*nf]<-log(Mod(complex(real=mydata[,f+4*nf],imaginary=mydata[,f+5*nf]))^2)
    # 
    #   mycdf<-kcde(mydata[,f],gridsize=1000,eval.points=mydata[,f],xmin=min(mydata[,f])-sd(mydata[,f]),xmax=max(mydata[,f])+sd(mydata[,f]))
    #   mydata[,f]<-qnorm(mycdf$estimate)
    # 
    #   # mycdf<-kcde(mydata[,f+2*nf],gridsize=1000,eval.points=mydata[,f+2*nf],xmin=min(mydata[,f+2*nf])-sd(mydata[,f+2*nf]),xmax=max(mydata[,f+2*nf])+sd(mydata[,f+2*nf]))
    #   # mydata[,f+2*nf]<-qnorm(mycdf$estimate)
    #   # 
    #   # mycdf<-kcde(mydata[,f+4*nf],gridsize=1000,eval.points=mydata[,f+4*nf],xmin=min(mydata[,f+4*nf])-sd(mydata[,f+4*nf]),xmax=max(mydata[,f+4*nf])+sd(mydata[,f+4*nf]))
    #   # mydata[,f+4*nf]<-qnorm(mycdf$estimate)
    #   colnames(mydata)[f] <- paste0("PP_f",toString(f))
    # }
    # # mydata <- mydata[,-((5*nf+1):(6*nf))]
    

    if(d==1){
      mouseDATA<-mydata
    } else{
      mouseDATA<-rbind(mouseDATA,mydata)
      rm(mydata)
    }
  }
  
  if(m==1){
    MYDATA<-mouseDATA
  } else{
    MYDATA<-rbind(MYDATA,mouseDATA)
    rm(mouseDATA)
  }
}

# Do if reading MYDATA and need to rescale running speed
for(d in 1:max(as.numeric(MYDATA$Day)))
{
  MYDATA$Running_Speed[MYDATA$Day==toString(d)]<-MYDATA$Running_Speed[MYDATA$Day==toString(d)]/max(MYDATA$Running_Speed[MYDATA$Day==toString(d)])
}

tdiff<-diff(MYDATA$Time)
tdiff<-c(0,tdiff)
MYDATA$ARStart<-rep(FALSE,length(MYDATA[,1]))
MYDATA$ARStart[abs(tdiff)>1]<-TRUE
MYDATA$ARStart[1]<-TRUE

MYDATA$MouseSession<-with(MYDATA,interaction(Mouse,Session))
MYDATA$DaySession<-with(MYDATA,interaction(Day,Session))
# MYDATA<-MYDATA[-c(51:200)]

#run and save gam for each frequency
setwd("C:/Data")
cwd<-getwd()
dir.create(file.path(cwd, "GAMS_full"), showWarnings = FALSE)
setwd(file.path(cwd, "GAMS_full"))
#save MYDATA
saveRDS(MYDATA,file="MYDATAPP_full.rds")
# read MYDATA if not running previous code
# MYDATA<-readRDS("MYDATA_full.rds")

# #create new data for predictions
# mina<- -3
# maxa<- 3
# nb=100
# abins<-seq(mina,maxa,length.out=nb)
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

#run gam on all frequencies
mb<-"cr"
rb<-"fs"
library(itsadug)
for(f in 30:nf)
{
  fmod <- bam(MYDATA[,f] ~ Session+
                           s(Running_Speed,bs="cr")+
                           s(Time,bs="cr",by=Session)+
                           s(FTphase,bs="cc")+
                           ti(Time,Running_Speed,bs="cr",by=Session,k=c(3,3))+
                           ti(Time,FTphase,bs=c("cr","cc"),by=Session,k=c(5,5))+
                           ti(Running_Speed,FTphase,bs=c("cr","cc"),k=c(5,5))+
                           s(Running_Speed,Day,bs="fs",xt="cr",m=1)+
                           s(Time,Trial,bs="fs",xt="cr",m=1)+
                           s(FTphase,Day,bs="fs",xt="cc",m=1),
                           method="fREML",
                           family=gaussian(),
                           data = MYDATA)
  
  r1 <- start_value_rho(fmod)
  
  filename<-paste0("VTP_bS_f",f,"_full_PP.rds")
  saveRDS(fmod,filename)
  # 
  # filename<-paste(names(MYDATA)[f],"_test.rds",sep="")
  # saveRDS(fmod,filename)
  # form<-as.formula(paste0("f20 ~ s(f",f,",bs=mb)+s(f",f,",Mouse,bs=rb,xt=mb)+s(f",f,",Day,bs=rb,xt=mb)"))
  # fmod <- bam(form,method="fREML",family=gaussian(),data = MYDATA)
  # fRsq[f]<-summary(fmod)$r.sq
  # colnames(newdat)[1]<-colnames(MYDATA)[f]
  # modpred[,(nv*(f-1)+1):(nv*f)]<-predict(fmod, newdat, type = "terms")
  # fmod<-lmer(MYDATA[,f] ~ 1+(1|Strain)+(1|Mouse)+(1|Day)+(1|DaySession),data = MYDATA)
  
  print(paste("Completed: ",f," of ",nf),quote=FALSE)
}

# #save model predictions to disk
# filename<-paste0("modpred",".rds")
# saveRDS(modpred,filename)
# 
# #save R-squared values to disk
# filename<-paste0("AARsq",".rds")
# saveRDS(fRsq,filename)

a<-matrix(0,
         nrow=1,
         ncol=2)

