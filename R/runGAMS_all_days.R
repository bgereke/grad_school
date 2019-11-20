#runGAMS_all_days

#load my functions
setwd("C:/Users/Brian/Documents/R")

#load libraries  
library(mgcv)
library(RevoUtilsMath)
library(ks)

#go to data directory
setwd("C:/Data")
miceDirs<-readLines("mice.txt",warn=FALSE,skipNul=TRUE)
numMice<-length(miceDirs)

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
    
    #load and merge data from sessions
    inlines<-readLines(infile,warn=FALSE,skipNul=TRUE)
    
    s1data<-read.table(paste(inlines[5],"\\CS\\rtable.csv",sep=""),sep=",",header=TRUE)
    s1data$Session <- factor(rep("one",length(s1data$Time)))
    
    s2data<-read.table(paste(inlines[6],"\\CS\\rtable.csv",sep=""),sep=",",header=TRUE)
    s2data$Session <- factor(rep("two",length(s2data$Time)))
    
    s3data<-read.table(paste(inlines[7],"\\CS\\rtable.csv",sep=""),sep=",",header=TRUE)
    s3data$Session <- factor(rep("three",length(s3data$Time)))
    
    mydata <- rbind(s1data,s2data,s3data)
    rm(s1data,s2data,s3data)
    mydata$Running_Speed <- abs(mydata$Running_Speed)
    mydata$Theta<-rowMeans(mydata[,4:6])
    mydata$Day<-factor(rep(toString(d),length(mydata[,1])))
    
    #transform response variables
    for(f in 1:50)
    {
      mycdf<-kcde(mydata[,f],gridsize=1000,eval.points=mydata[,f],xmin=min(mydata[,f])-sd(mydata[,f]),xmax=max(mydata[,f])+sd(mydata[,f]))
      mydata[,f]<-qnorm(mycdf$estimate)
    }
    #   mycdf<-kcde(mydata$Theta,gridsize=1000,eval.points=mydata$Theta,xmin=min(mydata$Theta)-sd(mydata$Theta),xmax=max(mydata$Theta)+sd(mydata$Theta))
    #   mydata$Theta<-qnorm(mycdf$estimate)
 
    if(d==1){
      MYDATA = mydata
    } else{
      MYDATA<-rbind(MYDATA,mydata)
      rm(mydata)
    }
    
  }
  
  MYDATA <- MYDATA[MYDATA$Running_Speed <= quantile(MYDATA$Running_Speed,probs=0.995),]
  
  #run and save gam for each frequency
  setwd(miceDirs[m])
  cwd<-getwd()
  dir.create(file.path(cwd, "GAMS"), showWarnings = FALSE)
  setwd(file.path(cwd, "GAMS"))
  for(f in 1:50)
  {
    #     if(f<3||f>6){
    #       fmod <- gam(mydata[,f] ~ s(Running_Speed)+s(Theta_Phase,bs="cc")+s(Position,bs="cc")+s(Time,by=Session)+s(Theta)
    #                   +ti(Running_Speed,Theta_Phase,bs=c("tp","cc"))+ti(Running_Speed,Position,bs=c("tp","cc"))+ti(Theta_Phase,Position,bs=c("cc","cc")),
    #                   method="REML",family=gaussian(),data = mydata) 
    #     } else {
    fmod <- bam(MYDATA[,f] ~ s(Running_Speed,k=9,fx=TRUE,by=Day)+s(Theta_Phase,bs="cc",k=9,fx=TRUE,by=Day)+s(Time,by=Session,k=9,fx=TRUE)
                # +ti(Running_Speed,Theta_Phase,bs=c("tp","cc"),k=10,fx=TRUE)+ti(Running_Speed,Position,bs=c("tp","cc"),k=10,fx=TRUE)+ti(Theta_Phase,Position,bs=c("cc","cc"),k=10,fx=TRUE)
                ,method="fREML",family=gaussian(),data = MYDATA) 
    # }
    
    filename<-paste(names(MYDATA)[f],".rds",sep="")
    
    saveRDS(fmod,filename)
    
    #numSmooths<-length(fmod$smooth)
    
    print(paste("Completed: ",2*f," Hz"),quote=FALSE)
  }
  write.table(MYDATA,file="mydata.csv",sep=",",col.names=NA,qmethod="double")
}
