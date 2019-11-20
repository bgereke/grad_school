# Calcium imaging place cell analysis
library(tidyverse)
library(mgcv)
library(wmtsa)
library(MASS)
library(data.table)
library(pracma)
library(signal)
library(circular)

#set directory
setwd("C:/Users/Brian/Downloads")

#load data
pdata <- read.table("position.csv",sep=",",header=FALSE)
rawdata <- fread("dFoF.csv")
rawdata <- as.data.frame(t(rawdata)) #want it in long format
numcells <- ncol(rawdata)

#interpolate to common time grid
Fp <- 0.032979438 #frame period
rawdata$time <- cumsum(rep(Fp,times = nrow(rawdata)))
rawdata <- rawdata[rawdata$time<=max(pdata$V1),] #remove extra imaging time points
rawdata$position <- approx(pdata$V1,pdata$V2,rawdata$time)$y

#rescale position and get running speed
unwrpos <- unwrap(2*pi*rawdata$position-pi)
sp = 3/max(rawdata$time) #smoothing span
unwrpos <- loess(unwrpos~time,span=sp,degree=1,data=rawdata)$fitted
vel <- abs(diff(unwrpos))/diff(rawdata$time)
rawdata$position <- 2*pi*rawdata$position-pi #radians
rawdata$vel <- c(vel[1],vel)/(2*pi)*200 #cm/sec

#baseline correct, denoise, detect transients, compute ratemaps, mvl, spatial info
procdata <- rawdata #processed data
bindata <- rawdata #binary data
bindata[,1:(ncol(rawdata)-3)] <- 0
numbins <- 200
grid <- seq(-pi,pi,length.out = numbins)
ratemaps <- matrix(0,nrow=numbins,ncol=numcells)
mvl <- rep(0,numcells)
spinfo <- rep(0,numcells)
placecells <- rep(0,numcells)
rmth <- 0.5
sp <- 30/max(rawdata$time) #smoothing span
bw <- 10 #smoothing bandwidth in cm for ratemaps
kappa <- log(2)/(1-cos(bw/200*pi)) #von Mises concentration parameter
posden<-density.circular(rawdata$position,
                         bw=kappa,
                         z=grid,
                         kernel="vonmises")  #circular density of position samples
px <- posden$y/sum(posden$y) #occupancy probability distribution
weights <- 1/px #weights for mvl
for(c in 1:numcells){
  #1st degree robust loess
  baseline <- loess(rawdata[,c]~time,span=sp,degree=1,family="symmetric",data=rawdata)$fitted
  procdata[,c] <- (rawdata[,c] - baseline)/baseline #dF/F
  #wavelet shrinkage denoising
  procdata[,c] <- wavShrink(procdata[,c],wavelet="s2",shrink.fun="mid",xform="modwt")
  #transient detection
  th <- abs(min(procdata[,c]))
  pks <- findpeaks(as.numeric(procdata[,c]),zero="-",minpeakheight=th)
  bindata[pks[,2],c] <- 1
  if(sum(bindata[,c])>3){ #want at least 3 spikes
    #ratemap
    spkden<-density.circular(rawdata$position[pks[,2]],
                             bw=kappa,
                             z=grid,
                             kernel="vonmises") #circular density of spikes
    ratemaps[,c] <- (spkden$n*spkden$y/sum(spkden$y))/(Fp*posden$n*posden$y/sum(posden$y))
    
    #mean vector length
    cweights <- approx(grid,weights,rawdata$position[pks[,2]])$y
    comp <- complex(arg = rawdata$position[pks[,2]])
    mvl[c] <- Mod(sum(cweights*comp)/sum(cweights))
    
    #spatial information
    trate <- sum(px*ratemaps[,c])
    spinfo[c] <- sum(ratemaps[,c]*log2(ratemaps[,c]/trate)*px)
    # spinfo[c] <- spinfo[c]/trate
    
    #is it a place cell?
    rmpks <- findpeaks(ratemaps[,c])
    if(sum(rmpks[,1]>rmth)>0){
      placecells[c] <- 1
    }
  }
}


# #make data frame for detected transients
# c <- 1
# pos <- circular(2*pi*bindata$position[bindata[,c]==1]/200-pi,type="angles",units="radians")
# med <- as.numeric(mean.circular(pos))
# pos <- as.numeric(pos)
# pos <- (unwrap(pos)-med-pi) %% (2*pi) - pi
# detdata <- data.frame(time=bindata$time[bindata[,c]==1],
#                       position=pos,
#                       cell=c)
# for(c in 2:numcells){
#   if(sum(bindata[,c])>10){
#     pos <- circular(2*pi*bindata$position[bindata[,c]==1]/200-pi,type="angles",units="radians")
#     med <- as.numeric(mean.circular(pos))
#     pos <- as.numeric(pos)
#     pos <- (unwrap(pos)-med-pi) %% (2*pi) - pi 
#     detdata <- rbind(detdata,
#                      data.frame(time=bindata$time[bindata[,c]==1],
#                                 position=pos,
#                                 cell=c))
#   }
# }
# ggplot(data = detdata, aes(x=time,y=position)) +
#   geom_point() +
#   geom_density2d() +
#   geom_smooth(method="gam",formula=y~s(x,bs="cr")) 