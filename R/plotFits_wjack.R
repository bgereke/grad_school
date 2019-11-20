library(mgcv)
library(RevoUtilsMath)
library(ggplot2)
library(MASS)
library(itsadug)
library(RColorBrewer)
library(directlabels)
library(scales)
library(akima)
# library(doParallel)

# if(require("RevoUtilsMath")){
#   setMKLthreads(3)
#   print("threads set")
# }

setwd("C:/Data/GAMS_full")

#load my data
cwd<-getwd()
# setwd(file.path(cwd, "GAMS"))
# MYDATA<-read.table("MYDATA.csv", header = TRUE, sep = ",", row.names = 1)
MYDATA<-readRDS("MYDATA_full.rds ")
MYDATA<-MYDATA[c(101:150,(ncol(MYDATA)-10):ncol(MYDATA))] #CICS only
#MYDATA<-MYDATA[-c(51:(ncol(MYDATA)-11))] #PICS only
#MYDATA<-MYDATA[c(201:250,(ncol(MYDATA)-10):ncol(MYDATA))] #PP only

#Set model prediction vars
numdays<-length(unique(MYDATA$Day))
days<-seq(1,numdays)
nb<-100
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


#create new data
dgrid<-expand.grid(rsbins,tbins,unique(MYDATA$Session))
colnames(dgrid)<-c("runningspeed","time","session")
newdat<-with(MYDATA,data.frame(Running_Speed=dgrid$runningspeed,Time=dgrid$time,Session=dgrid$session,Trial=rep("1",length(dgrid[,1]))))

#preallocate predictions for full data
frsfit<-matrix(0,nrow=nb*nf,ncol=3)
ftimefit<-matrix(0,nrow=nb*nf*3,ncol=4)
fintfit<-matrix(0,nrow=nb*nb*nf*3,ncol=5)
#fRsq<-matrix(0,nrow=1,ncol=nf)

#set known values
frsfit[,2]<-rep(freq,each=nb) 
frsfit[,3]<-rep(rsbins,times=nf)
ftimefit[,2]<-rep(unique(MYDATA$Session),each=nb*nf)
ftimefit[,3]<-rep(rep(freq,each=nb),times=3) 
ftimefit[,4]<-rep(tbins,times=3*nf)
fintfit[,2]<-rep(unique(MYDATA$Session),each=nb*nb*nf)
fintfit[,3]<-rep(rep(freq,each=nb*nb),times=3)
fintfit[,4]<-rep(rep(tbins,each=nb),times=nf*3)
fintfit[,5]<-rep(rsbins,times=nb*nf*3)

# rsse<-rsfit
# timese<-timefit
# intse<-intfit

#get full data predictions
for(f in 1:nf)
{
  filename<-paste0("VTbSIbS_f",f,"_full_CICS.rds")
  fmod<-readRDS(filename)
  #fRsq[f]<-summary(fmod)$r.sq
  pred<-predict.bam(fmod, newdat, type = "terms")

  #running speed plot
  frsfit[((f-1)*nb+1):(f*nb),1]<-pred[1:nb,1]
  # rsse[((f-1)*nb+1):(f*nb),1]<-pred$se.fit[1:nb,1]

  #time by session plot
  ftimefit[((f-1)*nb+1):(f*nb),1]<-pred[seq(1,nb*nb,by=nb),2] #session 1
  ftimefit[((f-1)*nb+1+nb*nf):(f*nb+nb*nf),1]<-pred[seq(nb*nb+1,2*nb*nb,by=nb),3] #session 2
  ftimefit[((f-1)*nb+1+2*nb*nf):(f*nb+2*nb*nf),1]<-pred[seq(2*nb*nb+1,3*nb*nb,by=nb),4] #session 3
  # timese[((f-1)*nb+1):(f*nb),1]<-pred$se.fit[seq(1,nb*nb,by=nb),2] #session 1
  # timese[((f-1)*nb+1+nb*nf):(f*nb+nb*nf),1]<-pred$se.fit[seq(nb*nb+1,2*nb*nb,by=nb),3] #session 2
  # timese[((f-1)*nb+1+2*nb*nf):(f*nb+2*nb*nf),1]<-pred$se.fit[seq(2*nb*nb+1,3*nb*nb,by=nb),4] #session 3

  #time/running speed interactions by session plot
  fintfit[((f-1)*nb*nb+1):(f*nb*nb),1]<-pred[1:(nb*nb),5] #session 1
  fintfit[((f-1)*nb*nb+1+nb*nb*nf):(f*nb*nb+nb*nb*nf),1]<-pred[(nb*nb+1):(2*nb*nb),6] #session 2
  fintfit[((f-1)*nb*nb+1+2*nb*nb*nf):(f*nb*nb+2*nb*nb*nf),1]<-pred[(2*nb*nb+1):(3*nb*nb),7] #session 3
  # intse[((f-1)*nb*nb+1):(f*nb*nb),1]<-pred$se.fit[1:(nb*nb),5] #session 1
  # intse[((f-1)*nb*nb+1+nb*nb*nf):(f*nb*nb+nb*nb*nf),1]<-pred$se.fit[(nb*nb+1):(2*nb*nb),6] #session 2
  # intse[((f-1)*nb*nb+1+2*nb*nb*nf):(f*nb*nb+2*nb*nb*nf),1]<-pred$se.fit[(2*nb*nb+1):(3*nb*nb),7] #session 3

  print(paste("Completed: ",f," of ",nf),quote=FALSE)
}

#get day predictions in parallel
# cl <- makeCluster(detectCores(), type='PSOCK')
# registerDoParallel(cl)

setwd("C:/Data/GAMS_CV")
count<-0
results<-list(0)

# results<-foreach(d = 1:numdays, .packages='mgcv') %dopar%
# {
for(d in 1:numdays)
{
  #change trial
  newdat$Trial<-MYDATA$Trial[MYDATA$Day!=toString(d)][1]
  
  #preallocate plot data
  drsfit<-matrix(0,nrow=nb*nf,ncol=1)
  dtimefit<-matrix(0,nrow=nb*nf*3,ncol=1)
  dintfit<-matrix(0,nrow=nb*nb*nf*3,ncol=1)
  # dRsq<-matrix(0,nrow=1,ncol=nf)
  
  # rsse<-rsfit
  # timese<-timefit
  # intse<-intfit
  
  
  for(f in 1:nf)
  {
    filename<-paste0("VTbSIbS_f",f,"_d",d,"_CICS.rds")
    fmod<-readRDS(filename)
    # dRsq[f]<-summary(fmod)$r.sq
    pred<-predict.bam(fmod, newdat, type = "terms")
    
    #running speed plot
    drsfit[((f-1)*nb+1):(f*nb)]<-pred[1:nb,1]
    # rsse[((f-1)*nb+1):(f*nb)]<-pred$se.fit[1:nb,1]
    
    #time by session plot
    dtimefit[((f-1)*nb+1):(f*nb)]<-pred[seq(1,nb*nb,by=nb),2] #session 1
    dtimefit[((f-1)*nb+1+nb*nf):(f*nb+nb*nf)]<-pred[seq(nb*nb+1,2*nb*nb,by=nb),3] #session 2
    dtimefit[((f-1)*nb+1+2*nb*nf):(f*nb+2*nb*nf)]<-pred[seq(2*nb*nb+1,3*nb*nb,by=nb),4] #session 3
    # timese[((f-1)*nb+1):(f*nb)]<-pred$se.fit[seq(1,nb*nb,by=nb),2] #session 1
    # timese[((f-1)*nb+1+nb*nf):(f*nb+nb*nf)]<-pred$se.fit[seq(nb*nb+1,2*nb*nb,by=nb),3] #session 2
    # timese[((f-1)*nb+1+2*nb*nf):(f*nb+2*nb*nf)]<-pred$se.fit[seq(2*nb*nb+1,3*nb*nb,by=nb),4] #session 3
    
    #time/running speed interactions by session plot
    dintfit[((f-1)*nb*nb+1):(f*nb*nb)]<-pred[1:(nb*nb),5] #session 1
    dintfit[((f-1)*nb*nb+1+nb*nb*nf):(f*nb*nb+nb*nb*nf)]<-pred[(nb*nb+1):(2*nb*nb),6] #session 2
    dintfit[((f-1)*nb*nb+1+2*nb*nb*nf):(f*nb*nb+2*nb*nb*nf)]<-pred[(2*nb*nb+1):(3*nb*nb),7] #session 3
    # intse[((f-1)*nb*nb+1):(f*nb*nb)]<-pred$se.fit[1:(nb*nb),5] #session 1
    # intse[((f-1)*nb*nb+1+nb*nb*nf):(f*nb*nb+nb*nb*nf)]<-pred$se.fit[(nb*nb+1):(2*nb*nb),6] #session 2
    # intse[((f-1)*nb*nb+1+2*nb*nb*nf):(f*nb*nb+2*nb*nb*nf)]<-pred$se.fit[(2*nb*nb+1):(3*nb*nb),7] #session 3
    
    count<-count+1
    print(paste("Completed: ",count," of ",nf*numdays),quote=FALSE)
  }
  results<-append(results,list(drsfit, dtimefit, dintfit))
}

# registerDoSEQ()

#get bias-corrected jackknife estimates
results<-results[2:length(results)] #only run once!
jrsfit<-matrix(0,nrow=nb*nf,ncol=1)
jtimefit<-matrix(0,nrow=nb*nf*3,ncol=1)
jintfit<-matrix(0,nrow=nb*nb*nf*3,ncol=1)
n<-length(MYDATA[,1])
for(d in 1:numdays){
  h<-sum(MYDATA$Day==toString(d))/n
  jrsfit<-jrsfit + (1-h)*results[[3*(d-1)+1]]
  jtimefit<-jtimefit + (1-h)*results[[3*(d-1)+2]]
  jintfit<-jintfit + (1-h)*results[[3*(d-1)+3]]
}
frsfit[,1]<-numdays*fresults[[1]][,1] - jrsfit
ftimefit[,1]<-numdays*fresults[[2]][,1] - jtimefit
fintfit[,1]<-numdays*fresults[[3]][,1] - jintfit

frsfit[,1]<-fresults[[1]][,1]
ftimefit[,1]<-fresults[[2]][,1]
fintfit[,1]<-fresults[[3]][,1]

#get bias-corrected jackknife ci's
jrsse<-matrix(0,nrow=nb*nf,ncol=1)
jtimese<-matrix(0,nrow=nb*nf*3,ncol=1)
jintse<-matrix(0,nrow=nb*nb*nf*3,ncol=1)
for(d in 1:numdays){
  h<-n/sum(MYDATA$Day==toString(d))
  jrsse<-jrsse+(1/(h-1))*(results[[3*(d-1)+1]]-frsfit[,1])^2
  jtimese<-jtimese+(1/(h-1))*(results[[3*(d-1)+2]]-ftimefit[,1])^2
  jintse<-jintse+(1/(h-1))*(results[[3*(d-1)+3]]-fintfit[,1])^2
}
jrsse<-sqrt(jrsse/numdays)
jtimese<-sqrt(jtimese/numdays)
jintse<-sqrt(jintse/numdays)

#r-squared plot
rsqdf<-data.frame(Frequency=freq,Rsquared=as.numeric(fRsq))
p<-ggplot(rsqdf,aes(x=Frequency,y=Rsquared))
p+geom_line()+geom_path(size=1)+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,0.2,by=0.025),name="R^2")+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())

#running speed plot
rbins<-seq(rsbins[1],rsbins[length(rsbins)],length=250)
fbins<-seq(2,100,length=250)
fitinterp<-interp(x=frsfit[,3],y=frsfit[,2],frsfit[,1],xo=rbins,yo=fbins)
highinterp<-interp(x=frsfit[,3],y=frsfit[,2],frsfit[,1]+2*jrsse,xo=rbins,yo=fbins)
lowinterp<-interp(x=frsfit[,3],y=frsfit[,2],frsfit[,1]-2*jrsse,xo=rbins,yo=fbins)
rfit <- data.frame(expand.grid(runningspeed = fitinterp$x, frequency = fitinterp$y), 
                   fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
cl<-max(abs(rfit$fit))
p<-ggplot(rfit,aes(x=runningspeed,y=frequency,fill=fit,z=fit))+
  geom_tile()+scale_x_continuous(breaks=seq(round(minrs),round(maxrs),by=0.2),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  # stat_contour(mapping=aes(x=runningspeed,y=frequency,z=high),colour="black",breaks=c(0),alpha=0.3)+
  # stat_contour(mapping=aes(x=runningspeed,y=frequency,z=low),colour="black",breaks=c(0),alpha=0.3)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

ggsave("RSfit.eps",p,width=8)


#time by session
tbins<-seq(tbins[1],tbins[length(tbins)],length=200)
fbins<-seq(2,100,length=200)
idxone<-ftimefit[,2] == "1"
idxtwo<-ftimefit[,2] == "2"
idxthree<-ftimefit[,2] == "3"
fitsone<-interp(x=ftimefit[idxone,4],y=ftimefit[idxone,3],ftimefit[idxone,1],xo=tbins,yo=fbins)
fitstwo<-interp(x=ftimefit[idxtwo,4],y=ftimefit[idxtwo,3],ftimefit[idxtwo,1],xo=tbins,yo=fbins)
fitsthree<-interp(x=ftimefit[idxthree,4],y=ftimefit[idxthree,3],ftimefit[idxthree,1],xo=tbins,yo=fbins)
tfitone <- data.frame(expand.grid(time = fitsone$x, frequency = fitsone$y), fit = c(fitsone$z))
tfittwo <- data.frame(expand.grid(time = fitstwo$x, frequency = fitstwo$y), fit = c(fitstwo$z))
tfitthree <- data.frame(expand.grid(time = fitsthree$x, frequency = fitsthree$y), fit = c(fitsthree$z))
tfitone$session<-rep(factor("1"),times=length(tfitone[,1]))
tfittwo$session<-rep(factor("2"),times=length(tfittwo[,1]))
tfitthree$session<-rep(factor("3"),times=length(tfitthree[,1]))
tfit<-rbind(tfitone,tfittwo,tfitthree)

cl<-max(abs(tfit$fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
ggplot(tfit,aes(x=time,y=frequency,fill=fit,z=fit))+ 
  geom_tile()+scale_x_continuous(breaks=seq(0,600,by=600),labels=seq(0,10,by=10),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",alpha=0.3,binwidth=0.05,size=0.05,show.legend = TRUE)+
  # stat_contour(mapping=aes(x=time,y=frequency,z=high),colour="black",breaks=c(0),alpha=0.3)+
  # stat_contour(mapping=aes(x=time,y=frequency,z=low),colour="black",breaks=c(0),alpha=0.3)+
  facet_wrap(~session,nrow=1)+theme(aspect.ratio = 1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank())

#differences time by session
tfit$diffs<-tfit$fit
tfit$diffs[tfit$session=="1"]<-tfit$fit[tfit$session=="1"] - tfit$fit[tfit$session=="2"]
tfit$diffs[tfit$session=="2"]<-tfit$fit[tfit$session=="1"] - tfit$fit[tfit$session=="3"]
tfit$diffs[tfit$session=="3"]<-tfit$fit[tfit$session=="2"] - tfit$fit[tfit$session=="3"]
tfit$dsessions<-as.factor(c(rep("1-2",times=200*200),rep("1-3",times=200*200),rep("2-3",times=200*200)))
cl<-max(abs(tfit$diffs))
p<-ggplot(tfit,aes(x=time,y=frequency,fill=diffs,z=diffs))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_tile()+scale_x_continuous(breaks=seq(0,600,by=600),labels=seq(0,10,by=10),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",alpha=0.3,binwidth=0.05,size=0.05,show.legend = TRUE)+
  # stat_contour(mapping=aes(x=time,y=frequency,z=high),colour="black",breaks=c(0),alpha=0.3)+
  # stat_contour(mapping=aes(x=time,y=frequency,z=low),colour="black",breaks=c(0),alpha=0.3)+
  facet_wrap(~dsessions,nrow=1)+theme(aspect.ratio = 1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank())  

#time/running speed interactions by frequency session 1
trsfit<-as.data.frame(fintfit)
colnames(trsfit)<-c("fit","session","frequency","time","runningspeed")
trsfit$frequency<-round(trsfit$frequency,digits=1)
trsfit$frequency<-as.ordered(trsfit$frequency)
for(f in 1:nf){
  levels(trsfit$frequency)[f]<-paste0(levels(trsfit$frequency)[f]," Hz")
  }
trsfit<-trsfit[trsfit$session==1,]
cl<-0.7*max(abs(trsfit$fit))
p<-ggplot(trsfit,aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=600),labels=seq(0,10,by=10),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  facet_wrap(~frequency,nrow=10) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"),
        strip.text.x = element_text(size = 7))

#time/running speed interactions by frequency session 2
trsfit<-as.data.frame(fintfit)
colnames(trsfit)<-c("fit","session","frequency","time","runningspeed")
trsfit$frequency<-round(trsfit$frequency,digits=1)
trsfit$frequency<-as.ordered(trsfit$frequency)
for(f in 1:nf){
  levels(trsfit$frequency)[f]<-paste0(levels(trsfit$frequency)[f]," Hz")
}
trsfit<-trsfit[trsfit$session==2,]
#cl<-0.7*max(abs(trsfit$fit))
p<-ggplot(trsfit,aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=600),labels=seq(0,10,by=10),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  facet_wrap(~frequency,nrow=10) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"),
        strip.text.x = element_text(size = 7))

#time/running speed interactions by frequency session 3
trsfit<-as.data.frame(fintfit)
colnames(trsfit)<-c("fit","session","frequency","time","runningspeed")
trsfit$frequency<-round(trsfit$frequency,digits=1)
trsfit$frequency<-as.ordered(trsfit$frequency)
for(f in 1:nf){
  levels(trsfit$frequency)[f]<-paste0(levels(trsfit$frequency)[f]," Hz")
}
trsfit<-trsfit[trsfit$session==3,]
#cl<-0.7*max(abs(trsfit$fit))
p<-ggplot(trsfit,aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  facet_wrap(~frequency,nrow=10) + 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"),
        strip.text.x = element_text(size = 7))

#time/running speed interactions by frequency session 1-2
trsfit<-as.data.frame(fintfit)
colnames(trsfit)<-c("fit","session","frequency","time","runningspeed")
trsfit$frequency<-round(trsfit$frequency,digits=1)
trsfit$frequency<-as.ordered(trsfit$frequency)
for(f in 1:nf){
  levels(trsfit$frequency)[f]<-paste0(levels(trsfit$frequency)[f]," Hz")
}
trsfit$fit[trsfit$session==1]<-trsfit$fit[trsfit$session==1]-trsfit$fit[trsfit$session==2]
trsfit<-trsfit[trsfit$session==1,]
cl<-0.7*max(abs(trsfit$fit))
p<-ggplot(trsfit,aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  facet_wrap(~frequency,nrow=10)+ 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"),
        strip.text.x = element_text(size = 7))

#time/running speed interactions by frequency session 1-3
trsfit<-as.data.frame(fintfit)
colnames(trsfit)<-c("fit","session","frequency","time","runningspeed")
trsfit$frequency<-round(trsfit$frequency,digits=1)
trsfit$frequency<-as.ordered(trsfit$frequency)
for(f in 1:nf){
  levels(trsfit$frequency)[f]<-paste0(levels(trsfit$frequency)[f]," Hz")
}
trsfit$fit[trsfit$session==1]<-trsfit$fit[trsfit$session==1]-trsfit$fit[trsfit$session==3]
trsfit<-trsfit[trsfit$session==1,]
#cl<-0.7*max(abs(trsfit$fit))
p<-ggplot(trsfit,aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  facet_wrap(~frequency,nrow=10)+ 
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"),
        strip.text.x = element_text(size = 7))

#time/running speed interactions by frequency session 2-3
trsfit<-as.data.frame(fintfit)
colnames(trsfit)<-c("fit","session","frequency","time","runningspeed")
trsfit$frequency<-round(trsfit$frequency,digits=1)
trsfit$frequency<-as.ordered(trsfit$frequency)
for(f in 1:nf){
  levels(trsfit$frequency)[f]<-paste0(levels(trsfit$frequency)[f]," Hz")
}
trsfit$fit[trsfit$session==1]<-trsfit$fit[trsfit$session==2]-trsfit$fit[trsfit$session==3]
trsfit<-trsfit[trsfit$session==1,]
#cl<-0.7*max(abs(trsfit$fit))
p<-ggplot(trsfit,aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  facet_wrap(~frequency,nrow=10)+
  theme(aspect.ratio = 1,panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(color="black",fill=NA,size=0.5),
        strip.background=element_rect(color="black",fill="white"),
        strip.text.x = element_text(size = 7))