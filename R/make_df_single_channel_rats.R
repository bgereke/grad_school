#create grouped data frame 
cyclep<-NULL
#load my functions
setwd("C:/Users/Brian/Documents/R")

#load libraries  
# library(RevoUtilsMath)
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
setwd("G:/Rats")
cwd<-getwd()
miceDirs<-readLines("rats_CA3.txt",warn=FALSE,skipNul=TRUE)
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
  dates<-readLines("dates_CA3.txt",warn=FALSE,skipNul=TRUE)
  numDates<-length(dates)
  
  first_day <- TRUE
  
  for(d in 1:numDates)
  {
    
    print(dates[d])
    setwd(paste(miceDirs[m],"\\",dates[d],sep=""))
    # infile<-"infile - Copy.txt"
    infile<-"infile.txt"
    # infile<-"InFile.txt"
    
    #load and merge data from sessions
    inlines<-readLines(infile,warn=FALSE,skipNul=TRUE)
    inlines <- sub("C:\\\\Data","G:",inlines)
    inlines <- sub("G:\\\\Data","G:",inlines)
    
    #make sure data is available
    if(!file.exists(paste(inlines[5],"\\CS\\rtable_CSC_CCS.csv",sep=""))){
      next
    }
    
    for(s in 1:numSessions){
      
      #session 1 data
      #position data
      sdata<-read.table(paste(inlines[s+4],"\\CS\\rtable_tracking.csv",sep=""),sep=",",header=TRUE)
      sdata$Session <- factor(rep(toString(s),nrow(sdata)))
      sdata$Trial <- factor(rep(toString(trial),nrow(sdata)))
      # #spiking data
      # spdata<-read.table(paste(inlines[s+4],"\\CS\\rtable_TS.csv",sep=""),sep=",",header=TRUE)
      # spdata$CellID<-as.factor(spdata$CellID)
      # spdata$Session <- factor(rep(toString(s),nrow(spdata)))
      # spdata$Trial <- factor(rep(toString(trial),nrow(spdata)))
      #spectral data
      cscdata<-fread(paste(inlines[s+4],"\\CS\\rtable_CSC_CCS.csv",sep=""),sep=",",header=TRUE)
      cscdata<-as.data.frame(cscdata)
      cscdata <- cscdata[,colnames(cscdata)!="Wave_Phase"]
      # #assymetry indexes
      # asymmetry<-read.table(paste(inlines[s+4],"\\CS\\sym.csv",sep=""),sep=",",header=TRUE)
      # asymmetry$Session <- factor(rep(toString(s),nrow(asymmetry)))
      # asymmetry$Trial <- factor(rep(toString(trial),nrow(asymmetry)))
      
      gc()
      
      trial<-trial+1
      
      #smooth position samples and compute running speed/acceleration by finite differences
      # sdata$Position <- unwrap(sdata$Position)
      sp <- 1.5/max(sdata$Time)
      smPos <- loess(Position ~ Time,span = sp,data=sdata,family = "symmetric",degree=1)$fitted
      rs <- abs(diff(smPos))/diff(sdata$Time)
      rs <- c(rs[1],rs)
      smrs <- abs(loess(rs ~ Time,span = sp,data=sdata,family = "symmetric",degree=1)$fitted)
      acc <- diff(smrs)/diff(sdata$Time)
      acc <- c(acc[1],acc)
      
      # sdata$Position <- (smPos-pi) %% (2*pi) - pi
      sdata$Position <- smPos
      sdata$Running_Speed <- smrs
      sdata$Acceleration <- acc
      sdata$Day<-factor(rep(toString(d),length(sdata[,1])))
      
      sdata <- sdata[,colnames(sdata)!="Time"]
      sdata<-cbind(cscdata,sdata)
      
      #remove data close to reward location
      sdata<-sdata[sdata$Position>15 & sdata$Position<185,]
      
      # #transform/normalize response variables
      # for(f in 1:nf)
      # {
      #   # sdata[,f]<-log(sdata[,f])
      #   # mycdf<-kcde(sdata[,f],gridsize=1000,eval.points=sdata[,f],xmin=min(sdata[,f])-sd(sdata[,f]),xmax=max(sdata[,f])+sd(sdata[,f]))
      #   # sdata[,f]<-qnorm(mycdf$estimate)
      #   sdata[,f] <- scale(sdata[,f]) 
      #   # sdata[,f] <- sdata[,f]/mean(sdata[sdata$Running_Speed<2.5,f])
      # }
      
      if(s==1){
        # spkphase <- spdata
        mydata <- sdata
        # Asymm <- asymmetry
      } else {
        # spkphase <- rbind(spkphase,spdata)
        mydata <- rbind(mydata,sdata)
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
    
    # # transform/normalize response variables
    # for(f in 1:nf)
    # {
    #   # sdata[,f]<-log(sdata[,f])
    #   # mycdf<-kcde(sdata[,f],gridsize=1000,eval.points=sdata[,f],xmin=min(sdata[,f])-sd(sdata[,f]),xmax=max(sdata[,f])+sd(sdata[,f]))
    #   # sdata[,f]<-qnorm(mycdf$estimate)
    #   mydata[,f] <- scale(mydata[,f]) + mean(mydata[,f])
    #   # mydata[,f] <- mydata[,f]/mean(mydata[,f])
    # }
    
    #add factors
    # mydata$Day<-factor(rep(toString(day),length(mydata[,1])))
    mydata$Mouse<-factor(rep(toString(m),length(mydata[,1])))
    # Asymm$Day<-factor(rep(toString(day),length(Asymm[,1])))
    # Asymm$Mouse<-factor(rep(toString(m),length(Asymm[,1])))
    # if (m<=3)
    # {
    #   mydata$Strain<-factor(rep("hyb",length(mydata[,1])))
    #   # Asymm$Strain<-factor(rep("hyb",length(Asymm[,1])))
    # } else{
    #   mydata$Strain<-factor(rep("B6",length(mydata[,1])))
    #   # Asymm$Strain<-factor(rep("B6",length(Asymm[,1])))
    # }
    day<-day+1
    mydata <- mydata[mydata$Time <= tth,]  #set time threshold
    # Asymm <- Asymm[Asymm$Time <= tth,]
    
    if(first_day){
      mouseDATA<-mydata
      first_day <- FALSE
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

#Rescale running speed
MYDATA$Scale_Speed <- MYDATA$Running_Speed
MYDATA$Norm_Speed <- MYDATA$Running_Speed
for(d in 1:max(as.numeric(MYDATA$Day)))
{
  MYDATA$Scale_Speed[MYDATA$Day==toString(d)]<-scale(log(MYDATA$Running_Speed[MYDATA$Day==toString(d)]))
  MYDATA$Norm_Speed[MYDATA$Day==toString(d)]<-MYDATA$Running_Speed[MYDATA$Day==toString(d)]/max(MYDATA$Running_Speed[MYDATA$Day==toString(d)])
  MYDATA$Acceleration[MYDATA$Day==toString(d)]<-MYDATA$Acceleration[MYDATA$Day==toString(d)]/max(abs(MYDATA$Acceleration[MYDATA$Day==toString(d)]))
  
  #remove days that have bad tracking
  if(max(MYDATA$Running_Speed[MYDATA$Day==toString(d)])>100){
    MYDATA <- MYDATA[MYDATA$Day!=toString(d),]
  }
  
}

tdiff<-diff(MYDATA$Time)
tdiff<-c(0,tdiff)
MYDATA$ARStart<-rep(FALSE,length(MYDATA[,1]))
MYDATA$ARStart[abs(tdiff)>1]<-TRUE
MYDATA$ARStart[1]<-TRUE

MYDATA$MouseSession<-with(MYDATA,interaction(Mouse,Session))
MYDATA$DaySession<-with(MYDATA,interaction(Day,Session))

#save data frame
dir.create(file.path(cwd, "GAMS_full_CCS"), showWarnings = FALSE)
setwd(file.path(cwd, "GAMS_full_CCS"))
#save MYDATA
saveRDS(MYDATA,file="MYDATA_full_CCS.rds")