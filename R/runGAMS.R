runGAMS <- function(infile){
  
  #load libraries  
  library(mgcv)
  library(RevoUtilsMath)
  library(ks)
  
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
  mydata <- mydata[mydata$Running_Speed <= quantile(mydata$Running_Speed,probs=0.995),]
  
#   mydata[,1:50]<-scale(mydata[,1:50],center=TRUE,scale=TRUE)
#   mydata$Theta<-scale(mydata$Theta,center=TRUE,scale=TRUE)
  
  #transform response variables
  for(f in 1:50)
  {
    mycdf<-kcde(mydata[,f],gridsize=1000,eval.points=mydata[,f],xmin=min(mydata[,f])-sd(mydata[,f]),xmax=max(mydata[,f])+sd(mydata[,f]))
    mydata[,f]<-qnorm(mycdf$estimate)
  }
  mycdf<-kcde(mydata$Theta,gridsize=1000,eval.points=mydata$Theta,xmin=min(mydata$Theta)-sd(mydata$Theta),xmax=max(mydata$Theta)+sd(mydata$Theta))
  mydata$Theta<-qnorm(mycdf$estimate)
  
  #run and save gam for each frequency
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
      fmod <- bam(mydata[,f] ~ s(Running_Speed,k=12,fx=TRUE)+s(Theta_Phase,bs="cc",k=10,fx=TRUE)+s(Time,by=Session,k=10,fx=TRUE)
                  # +ti(Running_Speed,Theta_Phase,bs=c("tp","cc"),k=10,fx=TRUE)+ti(Running_Speed,Position,bs=c("tp","cc"),k=10,fx=TRUE)+ti(Theta_Phase,Position,bs=c("cc","cc"),k=10,fx=TRUE)
                  ,method="fREML",family=gaussian(),data = mydata) 
    # }
    
    filename<-paste(names(mydata)[f],".rds",sep="")
    
    saveRDS(fmod,filename)
    
    #numSmooths<-length(fmod$smooth)
    
    print(paste("Completed: ",2*f," Hz"),quote=FALSE)
  }
  write.table(mydata,file="mydata.csv",sep=",",col.names=NA,qmethod="double")
  #mydata<-read.table("mydata.csv", header = TRUE, sep = ",", row.names = 1)
}
