#plot spectra

#load data#runGAMS_all_mice

#load my functions
setwd("C:/Users/Brian/Documents/R")

#load libraries  
library(RevoUtilsMath)
library(mgcv)
library(ks)
library(signal)
# library(lme4)

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
    print(dates[d])
    setwd(paste(miceDirs[m],"\\",dates[d],sep=""))
    infile<-"infile - Copy.txt"
    #infile<-"infile.txt"
    # infile<-"InFile.txt"
    
    #load and merge data from sessions
    inlines<-readLines(infile,warn=FALSE,skipNul=TRUE)
    
    s1data<-read.table(paste(inlines[5],"\\CS\\rtable_CSC.csv",sep=""),sep=",",header=TRUE)
    s1data$Session <- factor(rep(toString(1),length(s1data$Time)))
    s1data$Trial <- factor(rep(toString(trial),length(s1data$Time)))
    trial<-trial+1
    
    s2data<-read.table(paste(inlines[6],"\\CS\\rtable_CSC.csv",sep=""),sep=",",header=TRUE)
    s2data$Session <- factor(rep(toString(2),length(s2data$Time)))
    s2data$Trial <- factor(rep(toString(trial),length(s2data$Time)))
    trial<-trial+1
    
    s3data<-read.table(paste(inlines[7],"\\CS\\rtable_CSC.csv",sep=""),sep=",",header=TRUE)
    s3data$Session <- factor(rep(toString(3),length(s3data$Time)))
    s3data$Trial <- factor(rep(toString(trial),length(s3data$Time)))
    trial<-trial+1
    
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
    
    mydata <- rbind(s1data,s2data,s3data)
    rm(s1data,s2data,s3data)
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
    
    #transform response variables
    for(f in 1:50)
    {
      mydata[,f]<-log(mydata[,f])
      mydata[,f+4*nf]<-log(Mod(complex(real=mydata[,f+4*nf],imaginary=mydata[,f+5*f]))^2)
      colnames(mydata)[f+4*nf] <- paste0("PP_f",toString(f))
    }
    mydata <- mydata[,-((5*nf+1):(6*nf))]
    
    
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

#get spectra
# idx<-c(1:50) #PICS
idx<-c(101:150) #CICS
# idx<-c(201:250) #PP
freq<-10^seq(log10(2),log10(100),length.out=50)
nd = length(unique(MYDATA$Day))
nm = length(unique(MYDATA$Mouse))
CICSd<-matrix(0,nrow=nf*nd,ncol=4) #means by day
CICSm<-matrix(0,nrow=nf*nm,ncol=3) #means by mouse
CICSg<-matrix(0,nrow=nf*nm,ncol=2) #means by mouse

for(d in 1:nd)
  {
    CICSd[((d-1)*nf+1):(d*nf),1]<-colMeans(MYDATA[MYDATA$Day==toString(d),idx]) #normal mean
    CICSd[((d-1)*nf+1):(d*nf),2]<-rep(unique(MYDATA$Mouse[MYDATA$Day==toString(d)]),times=nf)
    CICSd[((d-1)*nf+1):(d*nf),3]<-rep(unique(MYDATA$Day[MYDATA$Day==toString(d)]),times=nf)
    CICSd[((d-1)*nf+1):(d*nf),4]<-freq
}

for(m in 1:nm)
{
  CICSm[((m-1)*nf+1):(m*nf),1]<-colMeans(MYDATA[MYDATA$Mouse==toString(m),idx]) #normal mean
  CICSm[((m-1)*nf+1):(m*nf),2]<-rep(m,times=nf)
  CICSm[((m-1)*nf+1):(m*nf),3]<-freq
}

CICSg[,1]<-colMeans(MYDATA[,idx]) #normal mean
CICSg[,2]<-freq

#days data frame
cicsddf<-data.frame(CICSd)
colnames(cicsddf)<-c("CICS","mouse","day","frequency")
cicsddf$mouse<-as.factor(cicsddf$mouse)
cicsddf$day<-as.factor(cicsddf$day)

#mouse data frame
cicsmdf<-data.frame(CICSm)
colnames(cicsmdf)<-c("CICS","mouse","frequency")
cicsmdf$mouse<-as.factor(cicsmdf$mouse)

#group data frame
cicsgdf<-data.frame(CICSg)
colnames(cicsgdf)<-c("CICS","frequency")

#plot CICS spectra
ggplot(cicsddf,aes(x=frequency,y=CICS,group=day))+
  geom_line(aes(colour=mouse))+
  scale_color_brewer(palette="Set2")+
  geom_line(data=cicsmdf,aes(x=frequency,y=CICS,group=mouse,color=mouse),size=1.5)+
  geom_line(data=cicsgdf,aes(x=frequency,y=CICS),size=2,color="black")+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())
