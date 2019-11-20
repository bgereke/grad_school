library(mgcv)
library(RevoUtilsMath)
library(ggplot2)
library(MASS)
library(itsadug)
library(RColorBrewer)
library(directlabels)
library(scales)
library(akima)

# if(require("RevoUtilsMath")){
#   setMKLthreads(3)
#   print("threads set")
# }

setwd("C:/Data/GAMS_CV")

#load my data
cwd<-getwd()
# setwd(file.path(cwd, "GAMS"))
# MYDATA<-read.table("MYDATA.csv", header = TRUE, sep = ",", row.names = 1)
MYDATA<-readRDS("MYDATAPP_full.rds ")
#MYDATA<-MYDATA[c(101:150,(ncol(MYDATA)-10):ncol(MYDATA))] #CICS only
#MYDATA<-MYDATA[-c(51:(ncol(MYDATA)-11))] #PICS only
# MYDATA<-MYDATA[c(201:250,(ncol(MYDATA)-10):ncol(MYDATA))] #PP only

#Set model prediction vars
days<-unique(MYDATA$Day)
nd<-length(days)
trials<-unique(MYDATA$Trial)
nt<-length(trials)
nb<-60
nf<-50
mint<-0
maxt<-600
freq<-10^seq(log10(2),log10(100),length.out=nf)
tbins<-seq(mint,maxt,length.out=nb)
minrs<-min(MYDATA$Running_Speed)
maxrs<-max(MYDATA$Running_Speed)
rsbins<-seq(minrs,maxrs,length.out=nb)
minp<- -pi
maxp<- pi
pbins<-seq(minp,maxp,length.out=nb)

#create new data for running speed plot
rsdat<-with(MYDATA,data.frame(Running_Speed=rep(rsbins,times=nd),Day=rep(days,each=nb),
                              Time=rep(300,nd*nb),FTphase=rep(0,nd*nb),Session=rep("1",nd*nb),
                              Trial=rep("1",nd*nb)))

#create new data for theta phase plot
tpdat<-with(MYDATA,data.frame(FTphase=rep(pbins,times=nd),Day=rep(days,each=nb),
                              Time=rep(300,nd*nb),Running_Speed=rep(0.3,nd*nb),Session=rep("1",nd*nb),
                              Trial=rep("1",nd*nb)))

#create new data for theta phase plot
tdat<-with(MYDATA,data.frame(Time=rep(tbins,times=nt),Trial=rep(trials,each=nb),
                              FTphase=rep(0,nt*nb),Running_Speed=rep(0.3,nt*nb),Session=rep("1",nt*nb),
                              Day=rep("1",nt*nb)))

#preallocation space for plotting data
rsdf <- matrix(0,nrow = nf*nd*nb, ncol = 4)
tpdf <- matrix(0,nrow = nf*nd*nb, ncol = 4)
tdf <- matrix(0,nrow = nf*nt*nb, ncol = 4)

#set known values
rsdf[,2] <- rep(freq, each = nd*nb)
rsdf[,3] <- rep(rep(days, each = nb), times = nf)
rsdf[,4] <- rep(rsbins, times = nf*nd)

tpdf[,2] <- rep(freq, each = nd*nb)
tpdf[,3] <- rep(rep(days, each = nb), times = nf)
tpdf[,4] <- rep(pbins, times = nf*nd)

tdf[,2] <- rep(freq, each = nt*nb)
tdf[,3] <- rep(rep(trials, each = nb), times = nf)
tdf[,4] <- rep(tbins, times = nf*nt)

setwd("C:/Data/GAMS_full_rho")

#prediction using "iterms"
for(f in 1:nf)
{
  
  filename<-paste0("Base_RT_TP_f",f,"_full_rho_PP.rds")
  fmod<-readRDS(filename)

  rspred<-predict.bam(fmod, rsdat, type = "terms")
  tppred<-predict.bam(fmod, tpdat, type = "terms")
  tpred<-predict.bam(fmod, tdat, type = "terms")

  rswant <- grep("s\\(Running\\_Speed\\,Day\\)",colnames(rspred)) 
  tpwant <- grep("s\\(FTphase\\,Day\\)",colnames(tppred))
  twant <- grep("s\\(Time\\,Trial\\)",colnames(tpred))
  
  rsdf[((f-1)*nd*nb+1):(f*nd*nb),1] <- rspred[,rswant]
  tpdf[((f-1)*nd*nb+1):(f*nd*nb),1] <- tppred[,tpwant]
  tdf[((f-1)*nt*nb+1):(f*nt*nb),1] <- tpred[,twant]
    
  print(paste("Completed: ",f," of ",nf),quote=FALSE)
  
}


#make running speed plot
rsdf<-as.data.frame(rsdf)
colnames(rsdf)<-c("fit","frequency","day","runningspeed")

#single frequency
ggplot(rsdf[rsdf$frequency==freq[1],],aes(x=runningspeed,y=fit,group=day))+
geom_line(aes(color=day))+geom_hline(yintercept = 0,color="black")

#interpolate
irsbins<-seq(min(rsbins),max(rsbins),length=250)
fbins<-seq(min(freq),max(freq),length=250)

didx <- rsdf$day == days[1]
fit <- interp(x=rsdf$runningspeed[didx],y=rsdf$frequency[didx],rsdf$fit[didx],xo=irsbins,yo=fbins)
fitdf <- data.frame(expand.grid(runningspeed = fit$x, frequency = fit$y),fit = c(fit$z))
fitdf$day <- rep(factor(days[1]),nrow(fitdf))

for(d in 2:nd){
  
  didx <- rsdf$day == days[d]
  fit <- interp(x=rsdf$runningspeed[didx],y=rsdf$frequency[didx],rsdf$fit[didx],xo=irsbins,yo=fbins)
  tempdf <- data.frame(expand.grid(runningspeed = fit$x, frequency = fit$y),fit = c(fit$z))
  tempdf$day <- rep(factor(days[d]),nrow(tempdf))
  fitdf <- rbind(fitdf, tempdf)
  
}

#display plot
cl<-0.25*max(abs(fitdf$fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
ggplot(fitdf,aes(x=runningspeed,y=frequency,fill=fit,z=fit))+ 
  geom_raster()+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,1,by=0.5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  facet_wrap(~day,nrow=5)+theme(aspect.ratio = 1)+
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))


#make theta phase plot
tpdf<-as.data.frame(tpdf)
colnames(tpdf)<-c("fit","frequency","day","phase")

#single frequency
ggplot(tpdf[tpdf$frequency==freq[2],],aes(x=phase,y=fit,group=day))+
  geom_line(aes(color=day))+geom_hline(yintercept = 0,color="black")

#interpolate
ipbins<-seq(min(pbins),max(pbins),length=250)
fbins<-seq(min(freq),max(freq),length=250)

didx <- tpdf$day == days[1]
fit <- interp(x=tpdf$phase[didx],y=tpdf$frequency[didx],tpdf$fit[didx],xo=ipbins,yo=fbins)
fitdf <- data.frame(expand.grid(phase = fit$x, frequency = fit$y),fit = c(fit$z))
fitdf$day <- rep(factor(days[1]),nrow(fitdf))

for(d in 2:nd){
  
  didx <- tpdf$day == days[d]
  fit <- interp(x=tpdf$phase[didx],y=tpdf$frequency[didx],tpdf$fit[didx],xo=ipbins,yo=fbins)
  tempdf <- data.frame(expand.grid(phase = fit$x, frequency = fit$y),fit = c(fit$z))
  tempdf$day <- rep(factor(days[d]),nrow(tempdf))
  fitdf <- rbind(fitdf, tempdf)
  
}

#display plot
cl<-0.7*max(abs(fitdf$fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
ggplot(fitdf,aes(x=phase,y=frequency,fill=fit,z=fit))+ 
  geom_raster()+
  scale_x_continuous(breaks=seq(-pi,pi,by=pi),labels=c("-pi","0","pi"),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  facet_wrap(~day,nrow=5)+theme(aspect.ratio = 1)+
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))


#make time plot
tdf<-as.data.frame(tdf)
colnames(tdf)<-c("fit","frequency","trial","time")

#interpolate
itbins<-seq(min(tbins),max(tbins),length=250)
fbins<-seq(min(freq),max(freq),length=250)

tidx <- tdf$trial == trials[1]
fit <- interp(x=tdf$time[tidx],y=tdf$frequency[tidx],tdf$fit[tidx],xo=itbins,yo=fbins)
fitdf <- data.frame(expand.grid(time = fit$x, frequency = fit$y),fit = c(fit$z))
fitdf$trial <- rep(factor(trials[1]),nrow(fitdf))

for(t in 2:nt){
  
  tidx <- tdf$trial == trials[t]
  fit <- interp(x=tdf$time[tidx],y=tdf$frequency[tidx],tdf$fit[tidx],xo=itbins,yo=fbins)
  tempdf <- data.frame(expand.grid(time = fit$x, frequency = fit$y),fit = c(fit$z))
  tempdf$trial <- rep(factor(trials[t]),nrow(tempdf))
  fitdf <- rbind(fitdf, tempdf)
  
}

#display plot
cl<-0.7*max(abs(fitdf$fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
ggplot(fitdf,aes(x=time,y=frequency,fill=fit,z=fit))+ 
  geom_raster()+
  scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  facet_wrap(~trial,ncol=11)+theme(aspect.ratio = 1)+
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"))
