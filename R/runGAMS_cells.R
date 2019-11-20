#runGAMS_all_mice

#load my functions
setwd("C:/Users/Brian/Documents/R")

#load libraries  
library(RevoUtilsMath)
library(mgcv)
library(ks)
library(signal)
library(refund)
library(data.table)

setwd("C:/Data")

#load data to get place fields boundaries
data <- read.table("rtable_PlaceCells.csv",sep=",",header=TRUE)

#set appropriate variables to factors
data$CELL <- as.factor(data$CELL)
data$DAY <- as.factor(data$DAY) 
data$MOUSE <- as.factor(data$MOUSE)
data$SESSION <- as.factor(data$SESSION) 

#get mean rate for each cell
cells <- unique(data$CELL)
nc <- length(cells)
data$MURATES <- data$RATES
for(c in 1:nc){
  data$MURATES[data$CELL==cells[c]] <- mean(data$RATES[data$CELL==cells[c]])
}
data <- data[data$MURATES>3,] #remove low firing rate cells
cells <- unique(data$CELL)
nc <- length(cells)

#load mouse list
miceDirs<-readLines("mice.txt",warn=FALSE,skipNul=TRUE)
numMice<-length(miceDirs)
day<-1
trial<-1
rsth<-5    #running speed threshold
tth<-600   #time threshold
nf<-50     #number of frequencies
dt<-0.005  #bin width in seconds

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
    #spiking data
    sp1data<-read.table(paste(inlines[5],"\\CS\\rtable_TS.csv",sep=""),sep=",",header=TRUE)[,1:2]
    sp1data$CellID<-as.factor(sp1data$CellID)
    sp1data$Session <- factor(rep(toString(1),nrow(sp1data)))
    sp1data$Trial <- factor(rep(toString(trial),nrow(sp1data)))
    sp1data<-sp1data[sp1data$SpTime<=tth,]
    #spectral data
    csc1data<-fread(paste(inlines[5],"\\CS\\rtable_CSC.csv",sep=""),sep=",",header=TRUE,drop=c(1:100))
    csc1data<-as.data.frame(csc1data)
    csc1data<-csc1data[csc1data$Time>0&csc1data$Time<=tth,]
    gc()
    
    trial<-trial+1
    
    #session 2 data
    #position data
    s2data<-read.table(paste(inlines[6],"\\CS\\rtable_tracking.csv",sep=""),sep=",",header=TRUE)
    s2data$Session <- factor(rep(toString(2),nrow(s2data)))
    s2data$Trial <- factor(rep(toString(trial),nrow(s2data)))
    s2data<-s2data[s2data$Time<=tth,]
    #spiking data
    sp2data<-read.table(paste(inlines[6],"\\CS\\rtable_TS.csv",sep=""),sep=",",header=TRUE)[,1:2]
    sp2data$CellID<-as.factor(sp2data$CellID)
    sp2data$Session <- factor(rep(toString(2),nrow(sp2data)))
    sp2data$Trial <- factor(rep(toString(trial),nrow(sp2data)))
    sp2data<-sp2data[sp2data$SpTime<=tth,]
    #spectral data
    csc2data<-fread(paste(inlines[6],"\\CS\\rtable_CSC.csv",sep=""),sep=",",header=TRUE,drop=c(1:100))
    csc2data<-as.data.frame(csc2data)
    csc2data<-csc2data[csc2data$Time>0&csc2data$Time<=tth,]
    gc()
    
    trial<-trial+1
    
    #session 3 data
    #position data
    s3data<-read.table(paste(inlines[7],"\\CS\\rtable_tracking.csv",sep=""),sep=",",header=TRUE)
    s3data$Session <- factor(rep(toString(3),nrow(s3data)))
    s3data$Trial <- factor(rep(toString(trial),nrow(s3data)))
    s3data<-s3data[s3data$Time<=tth,]
    #spiking data
    sp3data<-read.table(paste(inlines[7],"\\CS\\rtable_TS.csv",sep=""),sep=",",header=TRUE)[,1:2]
    sp3data$CellID<-as.factor(sp3data$CellID)
    sp3data$Session <- factor(rep(toString(3),nrow(sp3data)))
    sp3data$Trial <- factor(rep(toString(trial),nrow(sp3data)))
    sp3data<-sp3data[sp3data$SpTime<=tth,]
    #spectral data
    csc3data<-fread(paste(inlines[7],"\\CS\\rtable_CSC.csv",sep=""),sep=",",header=TRUE,drop=c(1:100))
    csc3data<-as.data.frame(csc3data)
    csc3data<-csc3data[csc3data$Time>0&csc3data$Time<=tth,]
    gc()
    
    trial<-trial+1
    
    #create master data frame for spike trains
    dfTrain <- as.data.frame(rep(0,tth/dt*3))
    colnames(dfTrain) <- "Train"
    dfTrain$Session <- as.factor(rep(c("1","2","3"),each=tth/dt))
    dfTrain$Time <- rep(seq(dt,600,by=dt),times=3)
    dfTrain$Position <- rep(0,tth/dt*3)
    dfTrain$Running_Speed <- rep(0,tth/dt*3)
    dfTrain<-dfTrain[!(dfTrain$Time<min(csc1data$Time)&dfTrain$Session=="1"),]
    dfTrain<-dfTrain[!(dfTrain$Time<min(csc2data$Time)&dfTrain$Session=="2"),]
    dfTrain<-dfTrain[!(dfTrain$Time<min(csc3data$Time)&dfTrain$Session=="3"),]
    #get corresponding csc indices
    l1<-sum(dfTrain$Session=="1")
    l2<-sum(dfTrain$Session=="2")
    l3<-sum(dfTrain$Session=="3")
    idx1<-seq.int(1,nrow(csc1data),by=2) #by=2 for 1 ms bins
    idx2<-seq.int(1,nrow(csc2data),by=2)
    idx3<-seq.int(1,nrow(csc3data),by=2)
    csc1data<-csc1data[idx1,]
    csc2data<-csc2data[idx2,]
    csc3data<-csc3data[idx3,]
    #normalize cscs
    mod<-matrix(0,nrow=l1+l2+l3,ncol=nf)
    phase<-matrix(0,nrow=l1+l2+l3,ncol=nf)
    for(f in 1:nf){
      mod[,f]<-log(Mod(complex(real=rbind(c(csc1data[,f],csc2data[,f],csc3data[,f])),
                                imaginary=rbind(c(csc1data[,f],csc2data[,f],csc3data[,f]))))^2)
      phase[,f]<-Arg(complex(real=rbind(c(csc1data[,f],csc2data[,f],csc3data[,f])),
                              imaginary=rbind(c(csc1data[,f],csc2data[,f],csc3data[,f]))))
      mycdf<-kcde(mod[,f],gridsize=500,eval.points=mod[,f],xmin=min(mod[,f])-sd(mod[,f]),xmax=max(mod[,f])+sd(mod[,f]))
      mod[,f]<-qnorm(mycdf$estimate)
      print(paste0("Normalized: ",f," of ",nf))
    }
    rm(csc1data,csc2data,csc3data,mycdf)
    gc()
    #append cosine phases
    for(f in 1:nf){
      dfTrain$nv<-cos(phase[,f])
      names(dfTrain)[names(dfTrain)=='nv']<-paste0("cf",f)
    }
    #append sine phases
    for(f in 1:nf){
      dfTrain$nv<-sin(phase[,f])
      names(dfTrain)[names(dfTrain)=='nv']<-paste0("sf",f)
    }
    rm(mod,phase)
    gc()
    
    #smooth position samples and compute running speed by finite differences
    s1data$Position <- unwrap(s1data$Position)
    s2data$Position <- unwrap(s2data$Position)
    s3data$Position <- unwrap(s3data$Position)
    smPos1 <- bam(Position ~ s(Time,bs="cr",k=200,fx=TRUE),data = s1data)
    smPos2 <- bam(Position ~ s(Time,bs="cr",k=200,fx=TRUE),data = s2data)
    smPos3 <- bam(Position ~ s(Time,bs="cr",k=200,fx=TRUE),data = s3data)
    rm(s1data,s2data,s3data)
    eps <- 0.0001
    newd <- data.frame(Time=dfTrain$Time[dfTrain$Session=="1"])
    dfTrain$Position[dfTrain$Session=="1"] <- predict.bam(smPos1,newd)
    newd <- newd+eps
    predeps <- predict.bam(smPos1,newd)
    dfTrain$Running_Speed[dfTrain$Session=="1"] <- abs(predeps - dfTrain$Position[dfTrain$Session=="1"])/eps
    newd <- data.frame(Time=dfTrain$Time[dfTrain$Session=="2"])
    dfTrain$Position[dfTrain$Session=="2"] <- predict.bam(smPos2,newd)
    newd <- newd+eps
    dfTrain$Running_Speed[dfTrain$Session=="2"] <- abs(predict.bam(smPos2,newd) - dfTrain$Position[dfTrain$Session=="2"])/eps
    newd <- data.frame(Time=dfTrain$Time[dfTrain$Session=="3"])
    dfTrain$Position[dfTrain$Session=="3"] <- predict.bam(smPos3,newd)
    newd <- newd+eps
    dfTrain$Running_Speed[dfTrain$Session=="3"] <- abs(predict.bam(smPos3,newd) - dfTrain$Position[dfTrain$Session=="3"])/eps
    dfTrain$Position <- (dfTrain$Position-pi) %% (2*pi) - pi
    
    rm(smPos1,smPos2,smPos3)
    
    # #remove data close to reward location
    # vmap <- bam(dfTrain$Running_Speed ~ s(Position,bs="cc",k=30,fx=TRUE),data = mydata)
    # pmap <- seq(-pi,pi,length.out =1000)
    # newd <- data.frame(Position=pmap)
    # vpred <- predict.bam(vmap,newd)
    # rewloc <- pmap[which(vpred == min(vpred))]
    # delta<-abs(Arg(complex(argument=rewloc)*Conj(complex(argument=dfTrain$Position))))/(2*pi)
    # D<-96.5 #track diameter
    # rd<-15 #distance from reward loc
    # dfTrain<-dfTrain[delta>=rd/(pi*D),]
    
    #run and save gam for each cell
    # setwd("C:/Data")
    # cwd<-getwd()
    # dir.create(file.path(cwd, "GAMS_full"), showWarnings = FALSE)
    # setwd(file.path(cwd, "GAMS_full"))
    #save MYDATA
    # saveRDS(MYDATA,file="MYDATA_full.rds")
    # read MYDATA if not running previous code
    # MYDATA<-readRDS("MYDATA_full.rds")
    
    #run gam on all cells
    nc<-length(unique(sp1data$CellID))
    for(c in 1:nc)
    {
      #bin spike train for each session
      s1times<-sp1data$SpTime[sp1data$CellID==toString(c)]
      s2times<-sp2data$SpTime[sp2data$CellID==toString(c)]
      s3times<-sp3data$SpTime[sp3data$CellID==toString(c)]
      dfTrain$Train<-0
      for(s in 1:length(s1times)){
        idx<-which.min((s1times[s]-dfTrain$Time[dfTrain$Session=="1"])^2)
        dfTrain$Train[idx]<-dfTrain$Train[idx]+1
      }
      for(s in 1:length(s2times)){
        idx<-which.min((s2times[s]-dfTrain$Time[dfTrain$Session=="2"])^2)
        s1end<-max(which(dfTrain$Session=="1"))
        dfTrain$Train[s1end+idx]<-dfTrain$Train[s1end+idx]+1
      }
      for(s in 1:length(s3times)){
        idx<-which.min((s3times[s]-dfTrain$Time[dfTrain$Session=="3"])^2)
        s2end<-max(which(dfTrain$Session=="2"))
        dfTrain$Train[s2end+idx]<-dfTrain$Train[s2end+idx]+1
      }
      
      #create history terms for each session
      wmax<-50

      lagone <- embed(dfTrain$Train[dfTrain$Session=="1"],wmax+1)
      lagtwo <- embed(dfTrain$Train[dfTrain$Session=="2"],wmax+1)
      lagthree <- embed(dfTrain$Train[dfTrain$Session=="3"],wmax+1)
      lagone<-lagone[,2:ncol(lagone)]
      lagtwo<-lagtwo[,2:ncol(lagtwo)]
      lagthree<-lagthree[,2:ncol(lagthree)]
      History<-as.matrix(rbind(lagone,lagtwo,lagthree))
      Hbins<-matrix(log((1:wmax)+5),nrow(History),ncol(History),byrow = TRUE)

      rm(lagone,lagtwo,lagthree)
      gc()
      
      start2<-min(which(dfTrain$Session=="2"))
      start3<-min(which(dfTrain$Session=="3"))
      tempTrain <- dfTrain[-c(1:wmax,start2:(start2+wmax-1),start3:(start3+wmax-1)),]
      ridge<-diag(100)
      Fbins<-matrix(seq(log(2),log(100),length.out = nf),nrow(tempTrain),nf,byrow = TRUE)
      C<-as.matrix(tempTrain[, (5 + 1):(5 + nf)])
      S<-as.matrix(tempTrain[, (5 + nf + 1):(5 + 2 * nf)])
      
      fmod <- bam(Train ~ s(Fbins,by=C,bs="cr",k=10,fx=TRUE)+
                          s(Fbins,by=S,bs="cr",k=10,fx=TRUE),
                          family=poisson(link="log"),
                          data = tempTrain)
      gc()
      
      #for external terms
      nb<-200
      rsbins<-seq(min(tempTrain$Running_Speed),max(tempTrain$Running_Speed),length.out=nb)
      pbins<-seq(-pi,pi,length.out = nb)
      tbins<-seq(0,600,length.out=nb)
      newdat<-with(tempTrain,data.frame(Running_Speed=rep(rsbins,times=3),Position=rep(pbins,times=3),
                                        Session=rep(unique(tempTrain$Session),each=nb),Time=rep(tbins,times=3)))
      newdat$Hbins<-matrix(log((1:wmax)+5),nrow(newdat),ncol(Hbins),byrow = TRUE)
      newdat$History<-matrix(1,nrow(newdat),ncol(Hbins))
      
      #for history term
      newdat<-with(tempTrain,data.frame(Running_Speed=rep(1,times=wmax),Position=rep(1,times=wmax),
                                        Session=rep(unique(tempTrain$Session)[1],each=wmax),Time=rep(1,times=wmax)))
      newdat$Hbins<-log((1:wmax)+5)
      newdat$History<-rep(1,times=wmax)
      newdat$Fbins<-seq(log(2),log(100),length.out = nf)
      newdat$C<-rep(1,times=nf)
      newdat$S<-rep(1,times=nf)
      
      
      pred<-predict.bam(fmod, newdat, type = "terms")
      
      # filename<-paste0("f",f,".rds")
      # saveRDS(fmod,filename)
      
      print(paste("Completed: ",c," of ",nc),quote=FALSE)
      
      rm(fmod)
      gc()
    }
  }
}

