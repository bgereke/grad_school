#create grouped data frame 
cyclep<-NULL
#load my functions
setwd("C:/Users/Brian/Documents/R")

#load libraries  
library(mgcv)
library(ks)
library(signal)
library(data.table)
library(circular)
library(lokern)
library(tsModel)

# if(require("RevoUtilsMath")){
#   setMKLthreads(2)
#   print("threads set")
# }

#go to data directory
setwd("E:/For Brian")
# setwd("C:/Data")
miceDirs<-readLines("mice.txt",warn=FALSE,skipNul=TRUE)
numMice<-length(miceDirs)
numSessions<-3
day<-1
trial<-1
rsth<-5 #running speed threshold
tth<-600 #time threshold
nf<-50 #number of frequencies

for(m in 1:numMice)
{
  print(c("Reading data for: ",miceDirs[m]),quote=FALSE)
  setwd(miceDirs[m])
  
  # read the dates.txt file
  dates<-readLines("dates.txt",warn=FALSE,skipNul=TRUE)
  numDates<-length(dates)
  
  for(d in 1:numDates)
  {
    for(s in 1:numSessions){
      print(dates[d])
      setwd(paste(miceDirs[m],"\\",dates[d],sep=""))
      # infile<-"infile - Copy.txt"
      #infile<-"infile.txt"
      infile<-"InFile.txt"
      
      #load and merge data from sessions
      inlines<-readLines(infile,warn=FALSE,skipNul=TRUE)
      
      #session 1 data
      #position data
      sdata<-fread(paste(inlines[s+4],"\\CS\\rtable_tracking.csv",sep=""),sep=",",header=TRUE)
      sdata<-as.data.frame(sdata)
      sdata$Session <- factor(rep(toString(s),nrow(sdata)))
      sdata$Trial <- factor(rep(toString(trial),nrow(sdata)))
      # #spiking data
      # spdata<-read.table(paste(inlines[s+4],"\\CS\\rtable_TS.csv",sep=""),sep=",",header=TRUE)
      # spdata$CellID<-as.factor(spdata$CellID)
      # spdata$Session <- factor(rep(toString(s),nrow(spdata)))
      # spdata$Trial <- factor(rep(toString(trial),nrow(spdata)))
      #spectral data
      cscdata<-fread(paste(inlines[s+4],"\\CS\\rtable_CSC.csv",sep=""),sep=",",header=TRUE)
      cscdata<-as.data.frame(cscdata)
      cscdata <- cscdata[,colnames(cscdata)!="Wave_Phase"]
      # #assymetry indexes
      # asymmetry<-read.table(paste(inlines[s+4],"\\CS\\sym.csv",sep=""),sep=",",header=TRUE)
      # asymmetry$Session <- factor(rep(toString(s),nrow(asymmetry)))
      # asymmetry$Trial <- factor(rep(toString(trial),nrow(asymmetry)))
      # #lag data
      # maxdata<-fread(paste(inlines[s+4],"\\CS\\rtable_maxima.csv",sep=""),sep=",",header=TRUE)
      # maxdata<-as.data.frame(maxdata)
      # lagdata <- Lag(maxdata[,53],1:80)                                                 
      
      gc()
      
      trial<-trial+1
      
      #smooth position samples and compute running speed/acceleration by finite differences
      sdata$Position <- unwrap(sdata$Position)
      # sp <- 1/max(sdata$Time)
      # smPos <- loess(Position ~ Time,span = sp,data=sdata,family="symmetric",degree=1)$fitted
      # rs <- abs(diff(smPos))/diff(sdata$Time)
      # rs <- sqrt(c(rs[1],rs))
      # smrs <- abs(loess(rs ~ Time,span = sp,data=sdata)$fitted)
      # acc <- diff(smrs)/diff(sdata$Time)
      # acc <- c(acc[1],acc)
      
      smrs <- glkerns(sdata$Time,sdata$Position,deriv=1,bandwidth=0.5,x.out=sdata$Time)
      smrs <- sqrt(abs(smrs$est))
      
      # sdata$Position <- (smPos-pi) %% (2*pi) - pi
      # sdata$Position <- smPos
      sdata$Running_Speed <- smrs
      # sdata$Acceleration <- acc
      
      sdata <- sdata[,colnames(sdata)!="Time"]
      sdata<-cbind(cscdata,sdata)
      
      if(s==1){
        # spkphase <- spdata
        mydata <- sdata
        # Lagdata <- lagdata
        # Maxdata <- maxdata
        # Asymm <- asymmetry
      } else {
        # spkphase <- rbind(spkphase,spdata)
        mydata <- rbind(mydata,sdata)
        # Lagdata <- rbind(Lagdata,lagdata)
        # Maxdata <- rbind(Maxdata,maxdata)
        # Asymm <- rbind(Asymm,asymmetry)
      }
      
      rm(sdata,cscdata)
      # rm(spdata)
      gc()
      
    } #end sessions
    
    # #rotate spike phases to max place cell firing for day
    # phasegrid<-seq(-pi,pi,length.out=512)
    # den<-density.circular(spkphase$Theta_Phase,bw=10,from=circular(-pi),to=circular(pi))
    # dmax<-phasegrid[which.max(den$y)]
    # mydata$Theta_Phase <- (unwrap(mydata$Theta_Phase)-dmax-pi) %% (2*pi) - pi
    # # wden<-density.circular(spkphase$Wave_Phase,bw=10,from=circular(-pi),to=circular(pi))
    # # wdmax<-phasegrid[which.max(wden$y)]
    # # mydata$dtwphase <- (unwrap(mydata$Wave_Phase)-wdmax-pi) %% (2*pi) - pi
    # cyclep<-c(cyclep den)
    
    #add factors
    mydata$Day<-factor(rep(toString(day),length(mydata[,1])))
    mydata$Mouse<-factor(rep(toString(m),length(mydata[,1])))
    # Asymm$Day<-factor(rep(toString(day),length(Asymm[,1])))
    # Asymm$Mouse<-factor(rep(toString(m),length(Asymm[,1])))
    if (m<=3)
    {
      mydata$Strain<-factor(rep("hyb",length(mydata[,1])))
      # Asymm$Strain<-factor(rep("hyb",length(Asymm[,1])))
    } else{
      mydata$Strain<-factor(rep("alz",length(mydata[,1])))
      # Asymm$Strain<-factor(rep("B6",length(Asymm[,1])))
    }
    day<-day+1
    mydata <- mydata[mydata$Time <= tth,]  #set time threshold
    # Asymm <- Asymm[Asymm$Time <= tth,]
    
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
    
    #transform/normalize response variables
    # for(f in 1:nf)
    # {
    #   # mydata[,f]<-log(mydata[,f])
    #   # mycdf<-kcde(mydata[,f],gridsize=500,eval.points=mydata[,f],xmin=min(mydata[,f])-sd(mydata[,f]),xmax=max(mydata[,f])+sd(mydata[,f]))
    #   # mydata[,f]<-qnorm(mycdf$estimate)
    #   # mydata[,f] <- scale(mydata[,f])
    #   mydata[,f] <- mydata[,f]/mean(mydata[mydata$Running_Speed<0.05,f])
    # }
    
    if(d==1){
      mouseDATA<-mydata
      # mouseAsymm<-Asymm
    } else{
      mouseDATA<-rbind(mouseDATA,mydata)
      # mouseAsymm<-rbind(mouseAsymm,Asymm)
      rm(mydata)#,Asymm)
    }
  }#end days
  
  if(m==1){
    MYDATA<-mouseDATA
    # ASYMM<-mouseAsymm
  } else{
    MYDATA<-rbind(MYDATA,mouseDATA)
    # ASYMM<-rbind(ASYMM,mouseAsymm)
    rm(mouseDATA)#,mouseAsymm)
  }
}#end mice

# # Do if reading MYDATA and need to rescale running speed
# MYDATA$Norm_Speed <- MYDATA$Running_Speed
# MYDATA$Scale_Speed <- MYDATA$Running_Speed
# for(d in 1:max(as.numeric(MYDATA$Day)))
# {
#   MYDATA$Norm_Speed[MYDATA$Day==toString(d)]<-MYDATA$Running_Speed[MYDATA$Day==toString(d)]/max(MYDATA$Running_Speed[MYDATA$Day==toString(d)])
#   MYDATA$Scale_Speed[MYDATA$Day==toString(d)]<-scale(MYDATA$Running_Speed[MYDATA$Day==toString(d)])
#   MYDATA$Acceleration[MYDATA$Day==toString(d)]<-MYDATA$Acceleration[MYDATA$Day==toString(d)]/max(abs(MYDATA$Acceleration[MYDATA$Day==toString(d)]))
# }

tdiff<-diff(MYDATA$Time)
tdiff<-c(0,tdiff)
MYDATA$ARStart<-rep(FALSE,length(MYDATA[,1]))
MYDATA$ARStart[abs(tdiff)>1]<-TRUE
MYDATA$ARStart[1]<-TRUE

MYDATA$MouseSession<-with(MYDATA,interaction(Mouse,Session))
MYDATA$DaySession<-with(MYDATA,interaction(Day,Session))
MYDATA$StrainSession<-with(MYDATA,interaction(Strain,Session))

#save data frame
setwd("E:/For Brian")
cwd<-getwd()
dir.create(file.path(cwd, "GAMS_full"), showWarnings = FALSE)
setwd(file.path(cwd, "GAMS_full"))
#save MYDATA
saveRDS(MYDATA,file="MYDATAPOW_full.rds")