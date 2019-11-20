#runGAMS_all_mice

#go to library directory
setwd("/home1/02592/bgereke")

#load libraries  
#install.packages("ks",lib = "/work/02592/bgereke/R", repos = "http://cran.us.r-project.org")
library(mgcv)
library(parallel)
# library(ks,lib.loc = "/home1/02592/bgereke/R/x86_64-unknown-linux-gnu-library/3.0/")

#go to data directory
setwd("/work/02592/bgereke/Data")
basedir<-getwd()
# miceDirs<-readLines("mice.txt",warn=FALSE)
# numMice<-length(miceDirs)
# day<-1
# trial<-1
# rsth<-5 #running speed threshold
# tth<-600 #time threshold
# 
# for(m in 1:numMice)
#   {
#   setwd(paste(basedir,"/",substr(miceDirs[m],start=9,stop=nchar(miceDirs[m])),sep="")) #14 for rats
#   mousedir<-getwd()
#   print(c("Reading data for: ",getwd()),quote=FALSE)
#   
#   # read the dates.txt file
#   dates<-readLines("dates.txt",warn=FALSE)
#   numDates<-length(dates)
#   
#   for(d in 1:numDates)
#     {
#     print(dates[d])
#     setwd(paste(mousedir,"/",dates[d],sep=""))
#     infile<-"infile - Copy.txt"
#     #infile<-"infile.txt"
#     # infile<-"InFile.txt"
#     
#     #load and merge data from sessions
#     inlines<-readLines(infile,warn=FALSE)
#     
#     s1data<-read.table(paste(getwd(),"/",substr(inlines[5],nchar(inlines[5])-5,nchar(inlines[5])),"/CS/rtable.csv",sep=""),sep=",",header=TRUE)
#     s1data$Session <- factor(rep(toString(1),length(s1data$Time)))
#     s1data$Trial <- factor(rep(toString(trial),length(s1data$Time)))
#     trial<-trial+1
#     
#     s2data<-read.table(paste(getwd(),"/",substr(inlines[6],nchar(inlines[6])-5,nchar(inlines[5])),"/CS/rtable.csv",sep=""),sep=",",header=TRUE)
#     s2data$Session <- factor(rep(toString(2),length(s2data$Time)))
#     s2data$Trial <- factor(rep(toString(trial),length(s2data$Time)))
#     trial<-trial+1
#     
#     s3data<-read.table(paste(getwd(),"/",substr(inlines[7],nchar(inlines[7])-5,nchar(inlines[5])),"/CS/rtable.csv",sep=""),sep=",",header=TRUE)
#     s3data$Session <- factor(rep(toString(3),length(s3data$Time)))
#     s3data$Trial <- factor(rep(toString(trial),length(s3data$Time)))
#     trial<-trial+1
#     
#     mydata <- rbind(s1data,s2data,s3data)
#     rm(s1data,s2data,s3data)
#     mydata$Running_Speed <- abs(mydata$Running_Speed)
#     mydata$Theta<-rowMeans(mydata[,4:6])
#     mydata$Day<-factor(rep(toString(day),length(mydata[,1])))
#     mydata$Mouse<-rep(toString(m),length(mydata[,1]))
#     day<-day+1
#     mydata <- mydata[mydata$Running_Speed <= quantile(mydata$Running_Speed,probs=0.995),]
#     # mydata <- mydata[mydata$Running_Speed >= rsth,]  #set running speed threshold
#     mydata <- mydata[mydata$Time <= tth,]  #set time threshold 
#     # mydata$Running_Speed<-scale(mydata$Running_Speed,center=TRUE,scale=TRUE)
#     
#     #transform response variables
#     for(f in 1:50)
#     {
#       mycdf<-kcde(mydata[,f],gridsize=1000,eval.points=mydata[,f],xmin=min(mydata[,f])-sd(mydata[,f]),xmax=max(mydata[,f])+sd(mydata[,f]))
#       mydata[,f]<-qnorm(mycdf$estimate)
#     }
#     #   mycdf<-kcde(mydata$Theta,gridsize=1000,eval.points=mydata$Theta,xmin=min(mydata$Theta)-sd(mydata$Theta),xmax=max(mydata$Theta)+sd(mydata$Theta))
#     #   mydata$Theta<-qnorm(mycdf$estimate)
#       # mycdf<-kcde(mydata$Running_Speed,gridsize=1000,eval.points=mydata$Running_Speed,xmin=min(mydata$Theta)-sd(mydata$Running_Speed),xmax=max(mydata$Running_Speed)+sd(mydata$Running_Speed))
#       # mydata$Running_Speed<-qnorm(mycdf$estimate)
#  
#     if(d==1){
#       mouseDATA<-mydata
#     } else{
#       mouseDATA<-rbind(mouseDATA,mydata)
#       rm(mydata)
#     }
#     
#   }
#   
#   if(m==1){
#     MYDATA<-mouseDATA
#   } else{
#     MYDATA<-rbind(MYDATA,mouseDATA)
#     rm(mouseDATA)
#   }
#   
# }

#run and save gam for each frequency
#dir.create(file.path(basedir, "GAMS_full_wint_rsnorm"), showWarnings = FALSE)
setwd(file.path(basedir, "GAMS_full_wint_rsnorm"))
#save MYDATA
#saveRDS(MYDATA,file="MYDATA.rds")
# read MYDATA if not running previous code
MYDATA<-readRDS("MYDATA.rds ")

# Do if reading MYDATA and need to rescale running speed
for(d in 1:max(as.numeric(MYDATA$Day)))
{
  MYDATA$Running_Speed[MYDATA$Day==toString(d)]<-MYDATA$Running_Speed[MYDATA$Day==toString(d)]/max(MYDATA$Running_Speed[MYDATA$Day==toString(d)])
}

#run gam on all frequencies
mb<-"cr"
# for(f in 1:50)
# {
mclapply(1:50,function(f)
{
  fmod <- bam(MYDATA[,f] ~ s(Running_Speed,bs=mb,k=8,fx=TRUE,by=Day)+s(Theta_Phase,bs="cc",k=6,fx=TRUE,by=Day)+s(Time,bs=mb,k=9,fx=TRUE,by=Session)
              +ti(Running_Speed,Theta_Phase,bs=c(mb,"cc"),fx=TRUE,k=c(4,4),by=Day)+ti(Time,Running_Speed,bs=c(mb,mb),k=c(6,3),fx=TRUE,by=Session)
              ,method="fREML",family=gaussian(),data = MYDATA) 
  
  filename<-paste(names(MYDATA)[f],".rds",sep="")
  
  saveRDS(fmod,filename)
  
  print(paste("Completed: ",2*f," Hz"),quote=FALSE)
},mc.cores=10)

