#load libraries  
library(RevoUtilsMath)
library(mgcv)
library(ks)
library(signal)
library(refund)
library(data.table)

#load data to get place fields boundaries, centers, and std's
setwd("C:/Data")
data <- read.table("rtable_PlaceCells.csv",sep=",",header=TRUE)

#set appropriate variables to factors
data$CELL <- as.factor(data$CELL)
data$MOUSE <- as.factor(data$MOUSE)
data$SESSION <- as.factor(data$SESSION) 
md<-with(data,interaction(MOUSE,DAY))
dd <- c(which(diff(data$DAY) != 0) +1, length(data$DAY)+1)
for(d in 1:(length(dd)-1)){data$DAY[dd[d]:(dd[d+1]-1)] <- d+1}
data$DAY <- as.factor(data$DAY)

#get cells that meet infield firing rate criterion
cells <- unique(data$CELL)
nc <- length(cells)
data$MURATES <- data$RATES
for(c in 1:nc){
  idx1 <- data$CELL==cells[c] & data$SESSION=="1"
  idx2 <- data$CELL==cells[c] & data$SESSION=="2"
  idx3 <- data$CELL==cells[c] & data$SESSION=="3"
  rate1 <- sum(data$NUMSPKS[idx1])/sum(data$ENTERTIME[idx1]+data$EXITTIME[idx1])
  rate2 <- sum(data$NUMSPKS[idx2])/sum(data$ENTERTIME[idx2]+data$EXITTIME[idx2])
  rate3 <- sum(data$NUMSPKS[idx3])/sum(data$ENTERTIME[idx3]+data$EXITTIME[idx3])
  data$MURATES[data$CELL==cells[c]] <- min(c(rate1,rate2,rate3))
}
data <- data[!is.nan(data$MURATES),]
data <- data[data$MURATES>1,] #remove low firing rate cells
cells <- unique(data$CELL) #cell list
nc <- length(cells)  #num cells
mice <- matrix(0,nrow=nc,ncol=1)
days <- mice
for(c in 1:nc){
  mice[c] <- data$MOUSE[data$CELL==cells[c]][1]
  days[c] <- data$DAY[data$CELL==cells[c]][1]
  }
mice <-as.factor(mice)
days <-as.factor(days)

#preallocate space for dataframe
tth<-600   #time threshold in seconds
nf<-50     #number of frequencies
dt<-0.001  #bin width in seconds
nt<-tth/dt #total number time points
nh<-round(0.2/dt) #number of spike history bins
zth <- 3 #place field boundaries (z-score)

#load mouse list
miceDirs<-readLines("mice.txt",warn=FALSE,skipNul=TRUE)
numMice<-length(miceDirs)

#build dataframe 
day<-1
trial<-1
cell<-1
rsth<-5    #running speed threshold (cm/sec)

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
    
    #smooth position samples and compute running speed by finite differences
    s1data$Position <- unwrap(s1data$Position)
    s2data$Position <- unwrap(s2data$Position)
    s3data$Position <- unwrap(s3data$Position)
    smPos1 <- bam(Position ~ s(Time,bs="cr",k=200,fx=TRUE),data = s1data)
    smPos2 <- bam(Position ~ s(Time,bs="cr",k=200,fx=TRUE),data = s2data)
    smPos3 <- bam(Position ~ s(Time,bs="cr",k=200,fx=TRUE),data = s3data)
    rm(s1data,s2data,s3data)
    eps <- 0.0001
    newd <- data.frame(Time=seq(dt,600,by=dt))
    newdt <- newd+eps
    Pos1 <- predict.bam(smPos1,newd)
    Vel1 <- abs(predict.bam(smPos1,newdt) - Pos1)/eps
    Pos2 <- predict.bam(smPos2,newd)
    Vel2 <- abs(predict.bam(smPos2,newdt) - Pos2)/eps
    Pos3 <- predict.bam(smPos3,newd)
    Vel3 <- abs(predict.bam(smPos3,newdt) - Pos3)/eps
    
    Pos1 <- (Pos1-pi) %% (2*pi) - pi
    Pos2 <- (Pos2-pi) %% (2*pi) - pi
    Pos3 <- (Pos3-pi) %% (2*pi) - pi
    
    Pos <- c(Pos1, Pos2, Pos3)
    Vel <- c(Vel1, Vel2, Vel3)
    
    rm(smPos1,smPos2,smPos3,Pos1,Pos2,Pos3,Vel1,Vel2,Vel3)
    gc()
    
    #bin spikes for each cell
    nc_d<-length(unique(sp1data$CellID))
    for(c in 1:nc_d)
    {
      if(sum(cells==toString(cell))>0){
        
        if(m==1 && d==1 && c==1){ 
          dfTrain <- data.frame("Train"=rep(0,nt*3))
          dfTrain$Session <- as.factor(rep(c("1","2","3"),each=nt))
          dfTrain$Time <- rep((0.5:nt)*dt,times=3)
          dfTrain$Position <- 0
          dfTrain$zPosition <- 0
          dfTrain$Running_Speed <- 0
          dfTrain$Cell <- as.factor(rep(toString(cell),times=3*nt))
          dfTrain$Mouse <- as.factor(rep(mice[m],times=3*nt))
          dfTrain$Day <- as.factor(rep(days[d],times=3*nt))
        } else{
          tmpTrain <- data.frame("Train"=rep(0,nt*3))
          tmpTrain$Session <- as.factor(rep(c("1","2","3"),each=nt))
          tmpTrain$Time <- rep((0.5:nt)*dt,times=3)
          tmpTrain$Position <- 0
          tmpTrain$zPosition <- 0
          tmpTrain$Running_Speed <- 0
          tmpTrain$Cell <- as.factor(rep(toString(cell),times=3*nt))
          tmpTrain$Mouse <- as.factor(rep(mice[m],times=3*nt))
          tmpTrain$Day <- as.factor(rep(days[d],times=3*nt))
          
          dfTrain <- rbind(dfTrain,tmpTrain)
          rm(tmpTrain)
          gc()
        }
        
        #spike times for each session
        s1times<-sp1data$SpTime[sp1data$CellID==toString(c)]
        s2times<-sp2data$SpTime[sp2data$CellID==toString(c)]
        s3times<-sp3data$SpTime[sp3data$CellID==toString(c)]

        #insert binned spikes, position, and velocity into dataframe
        cbool <- dfTrain$Cell==toString(cell)
        dfTrain$Train[cbool & dfTrain$Session=="1"] <- hist(s1times,breaks=seq(0,tth,by=dt),plot=FALSE)$counts
        dfTrain$Train[cbool & dfTrain$Session=="2"] <- hist(s2times,breaks=seq(0,tth,by=dt),plot=FALSE)$counts
        dfTrain$Train[cbool & dfTrain$Session=="3"] <- hist(s3times,breaks=seq(0,tth,by=dt),plot=FALSE)$counts
        dfTrain$Position[cbool] <- Pos
        dfTrain$Running_Speed[cbool] <- Vel/max(Vel)
        
        #rotate and scale position around place field center
        fcen <- data$CENTER[data$CELL==toString(cell)][1]
        fstd <- data$STD[data$CELL==toString(cell)][1]
        rotmat <- matrix(c(cos(-fcen),sin(-fcen),-sin(-fcen),cos(-fcen)),nrow = 2)
        rotp <- rotmat %*% t(cbind(cos(dfTrain$Position[cbool]),sin(dfTrain$Position[cbool])))
        dfTrain$Position[cbool] <- atan2(t(rotp[2,]),t(rotp[1,]))
        dfTrain$zPosition[cbool] <- dfTrain$Position[cbool]/fstd
        
        #rotate field boundaries
        left <- data$LEFT[data$CELL==toString(cell)][1]
        right <- data$RIGHT[data$CELL==toString(cell)][1]
        rotl <- rotmat %*% t(cbind(cos(left),sin(left)))
        rotr <- rotmat %*% t(cbind(cos(right),sin(right)))
        left <- as.numeric(atan2(t(rotl[2,]),t(rotl[1,])))
        right <- as.numeric(atan2(t(rotr[2,]),t(rotr[1,])))
        
        #get spike history
        bspk1 <- c(rep(0,times=nh), dfTrain$Train[cbool & dfTrain$Session=="1"])
        bspk2 <- c(rep(0,times=nh), dfTrain$Train[cbool & dfTrain$Session=="2"])
        bspk3 <- c(rep(0,times=nh), dfTrain$Train[cbool & dfTrain$Session=="3"])
        h1 <- embed(bspk1,nh+1)[,2:(nh+1)]
        h2 <- embed(bspk2,nh+1)[,2:(nh+1)]
        h3 <- embed(bspk3,nh+1)[,2:(nh+1)]
        if(m==1 && d==1 && c==1){
          History<-rbind(h1,h2,h3)
        }else{
          History<-rbind(History,h1,h2,h3)
        }
        
        #remove data outside place field
        # zbool <- dfTrain$zPosition>-zth & dfTrain$zPosition<zth
        # History <- History[zbool[1:nrow(History)],]
        # dfTrain <- dfTrain[zbool,]
        pbool <- dfTrain$Position<left | dfTrain$Position>right | dfTrain$Train==0 & cbool
        History <- History[!pbool,]
        dfTrain <- dfTrain[!pbool,]
      
        rm(s1times,s2times,s3times,rotp,bspk1,bspk2,bspk3,h1,h2,h3)
        gc()
      }
      cell <- cell + 1
    }
    rm(sp1data,sp2data,sp3data)
  }
}

zbool <- dfTrain$zPosition>-zth & dfTrain$zPosition<zth
History <- History[zbool[1:nrow(History)],]
dfTrain <- dfTrain[zbool,]

gc()

# for(i in 1:nrow(History)){History[i,]<-History[i,]*dfTrain$Running_Speed[i]}
# Hbins <- matrix(rep(seq(5,200,by=5),times=nrow(History)),nrow=nrow(History),ncol=ncol(History),byrow = TRUE)
# 
# mod3 <- pfr(Train ~ af(History,Hbins,k=c(5,7),bs=c('cr','cr')),
#                     family = poisson(link="log"),
#                     data = dfTrain)
# 
# mod <- bam(Train ~ Session*Cell +
#                    History +
#                    s(zPosition,bs='cr',by=Session) +
#                    s(Running_Speed,bs='cr') + 
#                    s(Time,bs='cr',by=Session) +
#                    ti(zPosition,Running_Speed,bs=c('cr','cr')) + 
#                    ti(Time,zPosition,bs=c('cr','cr'),by=Session),
#                    family = poisson(link="log"),
#                    data = dfTrain, samfrac = 0.1)
# 
# 
# mod2 <- bam(Train ~ Session + History*Running_Speed +
#              te(Time,zPosition,bs=c('cr','cr'),by=Session),
#              family = poisson(link="log"),
#              data = dfTrain, samfrac = 0.1)

df <- data.frame(Hist=as.vector(History),Vel=rep(dfTrain$Running_Speed,times=ncol(History)),Lags=rep(seq(1,200,by=1),each=nrow(History)))

mod <- bam(Hist~s(Lags,k=18,bs='cr',fx=TRUE) +
                ti(Lags,Vel,bs=c('cr','cr'),k=c(18,5),fx=TRUE),family=poisson(link="log"),data = df, samfrac = 0.1)

STA <- matrix(0,ncol=1,nrow=200)
Run <- matrix(0,ncol=200,nrow=100)
Place <- matrix(0,ncol=200,nrow=100)
newd <- data.frame(Running_Speed=seq(0,1,length.out=100),zPosition=seq(-3,3,length.out=100))
for(l in 1:200){
  lmod <- bam(History[,l]~s(Running_Speed,bs='cr') +
                           s(zPosition,bs='cr'),
                           data=dfTrain,family=poisson(link = "log"))
  STA[l] <- lmod$coefficients[1]
  pred <- predict(lmod,newd,type="terms")
  Run[,l] <- pred[,1]
  Place[,l] <- pred[,2]
}