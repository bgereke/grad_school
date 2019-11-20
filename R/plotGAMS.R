
  #load libraries  
  library(mgcv)
  library(RevoUtilsMath)
  library(ggplot2)
  library(MASS)
  library(itsadug)
  library(RColorBrewer)
  library(directlabels)
  library(scales)
  
  if(require("RevoUtilsMath")){
    setMKLthreads(3)
    print("threads set")
  }
  
  setwd("C:/Data/GAMS_full_wint_rsnorm_discrete")
  
  #load my data
  cwd<-getwd()
  # setwd(file.path(cwd, "GAMS"))
  # MYDATA<-read.table("MYDATA.csv", header = TRUE, sep = ",", row.names = 1)
  MYDATA<-readRDS("MYDATA.rds ")
  
  #only run to normalize running speed between 0 and 1
  for(d in 1:max(as.numeric(MYDATA$Day)))
  {
    MYDATA$Running_Speed[MYDATA$Day==toString(d)]<-MYDATA$Running_Speed[MYDATA$Day==toString(d)]/max(MYDATA$Running_Speed[MYDATA$Day==toString(d)])
  }
  
  #Set model prediciton vars
  numdays<-length(unique(MYDATA$Day))
  nb<-200
  nf<-50
  mint<-0
  maxt<-600
  day<-toString(10)
  minrs<-min(MYDATA$Running_Speed[MYDATA$Day==day])
  maxrs<-max(MYDATA$Running_Speed[MYDATA$Day==day])
  rsbins<-seq(minrs,maxrs,length.out=nb)
  tbins<-seq(0,600,length.out=nb)
  thbins<-seq(-pi,pi,length.out=nb)
  
  #make new data to predict responses from
  newdat<-with(MYDATA,data.frame(Running_Speed=c(rsbins,rep(mean(Running_Speed),4*nb)),
                                 Theta_Phase=c(rep(mean(Theta_Phase),nb),thbins,rep(mean(Theta_Phase),3*nb)),
                                 Time=c(rep(mean(Time),2*nb),rep(tbins,3)),
                                 Session=c(rep("2",2*nb),rep("1",nb),rep("2",nb),rep("3",nb)),
                                 Day=rep(day,5*nb)))
  #add grids for interaction terms
  inb = 200
  irsbins<-seq(minrs,maxrs,length.out=inb)
  ithbins<-seq(-pi,pi,length.out=inb)
  itbins<-seq(0,600,length.out=inb)
  trsgrid<-expand.grid(itbins,irsbins,unique(MYDATA$Session))
  colnames(trsgrid)<-c("t","rs","ses")
  rsthgrid<-expand.grid(irsbins,ithbins)
  colnames(rsthgrid)<-c("rs","thp")
  len_trs<-length(trsgrid$rs)
  len_rsth<-length(rsthgrid$rs)
  intdat<-with(MYDATA,data.frame(Running_Speed=c(trsgrid$rs,rsthgrid$rs),
                                 Theta_Phase=c(rep(mean(Theta_Phase),len_trs),rsthgrid$thp),
                                 Time=c(trsgrid$t,rep(mean(Time),len_rsth)),
                                 Session=c(trsgrid$ses,rep("2",len_rsth)),
                                 Day=rep(day,len_trs+len_rsth)))
  newdat<-rbind(newdat,intdat)
  len_trs<-len_trs/length(unique(MYDATA$Session))
  
  #preallocate predictions
  runningspeedFR<-matrix(nrow=nf,ncol=nb)
  thetaphaseFR<-matrix(nrow=nf,ncol=nb)
  sessiononeFR<-matrix(nrow=nf,ncol=nb)
  sessiontwoFR<-matrix(nrow=nf,ncol=nb)
  sessionthreeFR<-matrix(nrow=nf,ncol=nb)
  timerunningspeedoneFR<-matrix(nrow=nf,ncol=len_trs)
  timerunningspeedtwoFR<-matrix(nrow=nf,ncol=len_trs)
  timerunningspeedthreeFR<-matrix(nrow=nf,ncol=len_trs)
  runningspeedthetaphaseFR<-matrix(nrow=nf,ncol=len_rsth)
  fRsq<-matrix(nrow=1,ncol=nf) 
  
  #examine ACF plots
  # filename<-paste("f",2,".rds",sep="")
  # fmod<-readRDS(filename)
  # a<-acf_resid(fmod)
  # acf<-matrix(nrow=nf,ncol=length(a))
  # acf[1,]<-a
  # rm(a)
  # #load model for each frequency
  # for(f in 2:nf)
  # {
  #   filename<-paste("f",2*f,".rds",sep="")
  #   fmod<-readRDS(filename)
  #   acf[f,]<-acf_resid(fmod)
  # fRsq[f]<-summary(fmod)$r.sq
  # }
  
  #load model for each frequency
  for(f in 1:nf)
  {
    filename<-paste("f",2*f,".rds",sep="")
    fmod<-readRDS(filename)
    
    #assign(paste0("fmod", 2*f), fmod)
    # gam.check(fmod)
    fRsq[f]<-summary(fmod)$r.sq
    # fDev[f]<-summary(fmod)$dev.exp
    
    #lpmatrix for simultaneous confidence intervals
    newpred<-predict(fmod, newdat, type = "lpmatrix",discrete=TRUE)
    coefs<-coef(fmod)
    coefs<-as.matrix(coefs)
    
    #get Running Speed fits
    want<-grep("Running_Speed",colnames(newpred))
    dwant<-diff(want)
    didx<-which(dwant>1)
    want<-want[1:didx[1]]
    fits<-newpred[,want] %*% coefs[want,]
    runningspeedFR[f,]<-fits[(1):(nb),]
    
    #get Theta Phase fits
    want<-grep("Theta_Phase",colnames(newpred))
    dwant<-diff(want)
    didx<-which(dwant>1)
    want<-want[1:didx[1]]
    fits<-newpred[,want] %*% coefs[want,]
    thetaphaseFR[f,]<-fits[(nb+1):(2*nb),]
    
    #get Session One fits
    want<-grep("Time",colnames(newpred))
    dwant<-diff(want)
    didx<-which(dwant>1)
    want<-want[1:didx[1]]
    lwant<-length(want)
    numcols<-lwant/3
    want1<-want[1:numcols]
    fits<-newpred[,want1] %*% coefs[want1,]
    sessiononeFR[f,]<-fits[(2*nb+1):(3*nb),]
    #get Session Two fits
    want2<-want[(numcols+1):(2*numcols)]
    fits<-newpred[,want2] %*% coefs[want2,]
    sessiontwoFR[f,]<-fits[(3*nb+1):(4*nb),]
    #get Session Three fits
    want3<-want[(2*numcols+1):(3*numcols)]
    fits<-newpred[,want3] %*% coefs[want3,]
    sessionthreeFR[f,]<-fits[(4*nb+1):(5*nb),]
    
    #get interactions
    
    #get Time x Running_Speed Session 1
    want<-grep("Time,Running_Speed",colnames(newpred))
    lwant<-length(want)
    numcols<-lwant/3
    want1<-want[1:numcols]
    fits<-newpred[,want1] %*% coefs[want1,]
    timerunningspeedoneFR[f,]<-fits[(5*nb+1):(5*nb+len_trs),]
    #get Time x Running_Speed Session 2
    want2<-want[(numcols+1):(2*numcols)]
    fits<-newpred[,want2] %*% coefs[want2,]
    timerunningspeedtwoFR[f,]<-fits[(5*nb+1+len_trs):(5*nb+2*len_trs),]
    #get Time x Running_Speed Session 3
    want3<-want[(2*numcols+1):(3*numcols)]
    fits<-newpred[,want3] %*% coefs[want3,]
    timerunningspeedthreeFR[f,]<-fits[(5*nb+1+2*len_trs):(5*nb+3*len_trs),]
    
    #get Theta Phase x Running_Speed
    want<-grep("Running_Speed,Theta_Phase",colnames(newpred))
    lwant<-length(want)
    numcols<-lwant/3
    want1<-want[1:numcols]
    fits<-newpred[,want1] %*% coefs[want1,]
    runningspeedthetaphaseFR[f,]<-fits[(5*nb+3*len_trs+1):(5*nb+3*len_trs+len_rsth),]
    
    print(paste("Completed: ",2*f," Hz"),quote=FALSE)
  }
  
  #set tick values
  freq<-seq(2,100,by=2)
  time<-seq(0,600,length.out=nb)
  phase<-seq(-pi,pi,length.out=nb)
  # rsbins<-rsbins[rsbins<=max(MYDATA$Running_Speed[MYDATA$Day == day]) &
  #                  rsbins>=min(MYDATA$Running_Speed[MYDATA$Day == day])]
  itime<-seq(0,600,length.out=inb)
  iphase<-seq(-pi,pi,length.out=inb)
  # rsthgrid<-rsthgrid[rsthgrid$rs<=max(MYDATA$Running_Speed[MYDATA$Day == day]) &
  #                    rsthgrid$rs>=min(MYDATA$Running_Speed[MYDATA$Day == day])]
  
  #create data frames for ggplot2 figures
  
  #r-squared
  rsqdf<-data.frame(Frequency=freq,Rsquared=as.numeric(fRsq))
  p<-ggplot(rsqdf,aes(x=Frequency,y=Rsquared))
  p+geom_line()+geom_path(size=1)+
    scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,0.2,by=0.025),name="R^2")+
    theme(axis.line.x = element_line(color="black", size = 1),
          axis.line.y = element_line(color="black", size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank())
  
  #running speed
  runningspeeddf<-as.data.frame(matrix(t(runningspeedFR),nrow=length(rsbins)*nf,ncol=1))
  colnames(runningspeeddf)<-c("response")
  runningspeeddf$frequency<-rep(freq,rep(length(rsbins),nf))
  runningspeeddf$runningspeed<-rep(rsbins,nf)
  p<-ggplot(runningspeeddf,aes(x=runningspeed,y=frequency,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_tile()+scale_x_continuous(breaks=seq(round(minrs),round(maxrs),by=0.2),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-max(abs(runningspeeddf$response)),max(abs(runningspeeddf$response))),oob=squish)+
    stat_contour(colour="black",binwidth=0.075,size=0.05,alpha=0.3,show.legend = TRUE)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank())
  
  #theta phase
  thetaphasedf<-as.data.frame(matrix(t(thetaphaseFR),nrow=length(phase)*nf,ncol=1))
  colnames(thetaphasedf)<-c("response")
  thetaphasedf$frequency<-rep(freq,rep(length(phase),nf))
  thetaphasedf$thetaphase<-rep(phase,nf)
  p<-ggplot(thetaphasedf,aes(x=thetaphase,y=frequency,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_tile()+scale_x_continuous(breaks=seq(-pi,pi,pi/2),labels=c("-pi","-pi/2","0","pi/2","pi"),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-max(abs(thetaphasedf$response)),max(abs(thetaphasedf$response))),oob=squish)+
    stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank())
  
  #session one
  cm <- max(abs(timeonedf$response))
  
  timeonedf<-as.data.frame(matrix(t(sessiononeFR),nrow=length(time)*nf,ncol=1))
  colnames(timeonedf)<-c("response")
  timeonedf$frequency<-rep(freq,rep(length(time),nf))
  timeonedf$Time<-rep(time,nf)
  p<-ggplot(timeonedf,aes(x=Time,y=frequency,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_tile()+scale_x_continuous(breaks=seq(0,600,60),labels=seq(0,10,by=1),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cm,cm),oob=squish)+
    stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank())
  
  #session two
  timetwodf<-as.data.frame(matrix(t(sessiontwoFR),nrow=length(time)*nf,ncol=1))
  colnames(timetwodf)<-c("response")
  timetwodf$frequency<-rep(freq,rep(length(time),nf))
  timetwodf$Time<-rep(time,nf)
  p<-ggplot(timetwodf,aes(x=Time,y=frequency,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_tile()+scale_x_continuous(breaks=seq(0,600,60),labels=seq(0,10,by=1),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cm,cm),oob=squish)+
    stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank())
  
  #session three
  timethreedf<-as.data.frame(matrix(t(sessionthreeFR),nrow=length(time)*nf,ncol=1))
  colnames(timethreedf)<-c("response")
  timethreedf$frequency<-rep(freq,rep(length(time),nf))
  timethreedf$Time<-rep(time,nf)
  p<-ggplot(timethreedf,aes(x=Time,y=frequency,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_tile()+scale_x_continuous(breaks=seq(0,600,60),labels=seq(0,10,by=1),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cm,cm),oob=squish)+
    stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank())
  
  #session differences
  timeonetwodf<-as.data.frame(matrix(t(sessiononeFR-sessiontwoFR),nrow=length(time)*nf,ncol=1))
  timeonethreedf<-as.data.frame(matrix(t(sessiononeFR-sessionthreeFR),nrow=length(time)*nf,ncol=1))
  timetwothreedf<-as.data.frame(matrix(t(sessiontwoFR-sessionthreeFR),nrow=length(time)*nf,ncol=1))
  colnames(timeonetwodf)<-c("response")
  colnames(timeonethreedf)<-c("response")
  colnames(timetwothreedf)<-c("response")
  timeonetwodf$frequency<-rep(freq,rep(length(time),nf))
  timeonetwodf$Time<-rep(time,nf)
  timeonethreedf$frequency<-rep(freq,rep(length(time),nf))
  timeonethreedf$Time<-rep(time,nf)
  timetwothreedf$frequency<-rep(freq,rep(length(time),nf))
  timetwothreedf$Time<-rep(time,nf)
  cl<-max(abs(rbind(timeonetwodf$response,timeonethreedf$response)))
  p<-ggplot(timeonetwodf,aes(x=Time,y=frequency,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_tile()+scale_x_continuous(breaks=seq(0,600,60),labels=seq(0,10,by=1),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank())+
    ggtitle("Session One - Session Two")
  p<-ggplot(timeonethreedf,aes(x=Time,y=frequency,fill=response,z=response))
  p+geom_tile()+scale_x_continuous(breaks=seq(0,600,60),labels=seq(0,10,by=1),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank())+
    ggtitle("Session One - Session Three")
  p<-ggplot(timetwothreedf,aes(x=Time,y=frequency,fill=response,z=response))
  p+geom_tile()+scale_x_continuous(breaks=seq(0,600,60),labels=seq(0,10,by=1),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank())+
    ggtitle("Session Two - Session Three")
  
  #plot interactions
  
  #running speed x theta phase tiled by frequency
  runningspeedthetaphasedf<-as.data.frame(matrix(t(runningspeedthetaphaseFR),ncol=1))
  colnames(runningspeedthetaphasedf)<-c("response")
  runningspeedthetaphasedf$frequency<-rep(freq,rep(length(iphase)*length(irsbins),nf))
  runningspeedthetaphasedf$runningspeed<-rep(irsbins,length(iphase)*nf)
  runningspeedthetaphasedf$thetaphase<-rep(rep(iphase,rep(length(irsbins),length(iphase))),nf)
  cl<-max(abs(runningspeedthetaphasedf$response))
  p<-ggplot(runningspeedthetaphasedf,aes(x=runningspeed,y=thetaphase,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_raster()+scale_x_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
    scale_y_continuous(breaks=seq(-pi,pi,pi/2),labels=c("-pi","-pi/2","0","pi/2","pi"),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    facet_wrap(~frequency,nrow=5)+theme(aspect.ratio = 1)
  
  #running speed x theta phase tiled by running speed
  runningspeedthetaphasedf<-as.data.frame(matrix(t(runningspeedthetaphaseFR),ncol=1))
  colnames(runningspeedthetaphasedf)<-c("response")
  runningspeedthetaphasedf$frequency<-rep(freq,rep(length(iphase)*length(irsbins),nf))
  runningspeedthetaphasedf$runningspeed<-rep(irsbins,length(iphase)*nf)
  runningspeedthetaphasedf$thetaphase<-rep(rep(iphase,rep(length(irsbins),length(iphase))),nf)
  keep<-(round((1:nb/nb)/0.05,digits=3) %% 1) == 0
  keep[1]<-TRUE
  keep<-rep(keep,length(iphase)*nf)
  runningspeedthetaphasedf<-runningspeedthetaphasedf[keep,]
  cl<-max(abs(runningspeedthetaphasedf$response))
  p<-ggplot(runningspeedthetaphasedf,aes(x=thetaphase,y=frequency,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "RdBu")), space="Lab")
  p+geom_raster()+scale_x_continuous(breaks=seq(-pi,pi,pi/2),labels=c("-pi","-pi/2","0","pi/2","pi"),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    facet_wrap(~runningspeed)

  #session one x running speed tiled by frequency
  timerunningspeedonedf<-as.data.frame(matrix(t(timerunningspeedoneFR),ncol=1))
  colnames(timerunningspeedonedf)<-c("response")
  timerunningspeedonedf$frequency<-rep(freq,rep(length(itime)*length(irsbins),nf))
  timerunningspeedonedf$runningspeed<-rep(irsbins,length(itime)*nf)
  timerunningspeedonedf$time<-rep(rep(itime,rep(length(irsbins),length(itime))),nf)
  cl<-0.75*max(abs(timerunningspeedonedf$response))
  p<-ggplot(timerunningspeedonedf,aes(x=time,y=runningspeed,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_raster()+scale_x_continuous(breaks=seq(0,600,120),labels=seq(0,10,by=2),expand=c(0,0))+
    scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    facet_wrap(~frequency,nrow=5) + theme(aspect.ratio = 1)
  
  #session one x running speed tiled by time
  timerunningspeedonedf<-as.data.frame(matrix(t(timerunningspeedoneFR),ncol=1))
  colnames(timerunningspeedonedf)<-c("response")
  timerunningspeedonedf$frequency<-rep(freq,rep(length(itime)*length(irsbins),nf))
  timerunningspeedonedf$runningspeed<-rep(irsbins,length(itime)*nf)
  timerunningspeedonedf$time<-rep(rep(itime,rep(length(irsbins),length(itime))),nf)
  keep<-(round((1:nb/nb)/0.05,digits=3) %% 1) == 0
  keep[1]<-TRUE
  keep<-rep(rep(keep,rep(length(irsbins),length(keep))),nf)
  timerunningspeedonedf<-timerunningspeedonedf[keep,]
  timerunningspeedonedf$time<-round(timerunningspeedonedf$time)
  cl<-0.9*max(abs(timerunningspeedonedf$response))
  p<-ggplot(timerunningspeedonedf,aes(x=runningspeed,y=frequency,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_raster()+scale_x_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,100,by=20),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    facet_wrap(~time,ncol=3)+ theme(aspect.ratio = 1)
  
  #session two x running speed tiled by frequency
  timerunningspeedtwodf<-as.data.frame(matrix(t(timerunningspeedtwoFR),ncol=1))
  colnames(timerunningspeedtwodf)<-c("response")
  timerunningspeedtwodf$frequency<-rep(freq,rep(length(itime)*length(irsbins),nf))
  timerunningspeedtwodf$runningspeed<-rep(irsbins,length(itime)*nf)
  timerunningspeedtwodf$time<-rep(rep(itime,rep(length(irsbins),length(itime))),nf)
  cl<-0.75*max(abs(timerunningspeedonedf$response))
  p<-ggplot(timerunningspeedtwodf,aes(x=time,y=runningspeed,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_raster()+scale_x_continuous(breaks=seq(0,600,120),labels=seq(0,10,by=2),expand=c(0,0))+
    scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.2),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    facet_wrap(~frequency,nrow=5) + theme(aspect.ratio = 1)
  
  #session two x running speed tiled by time
  timerunningspeedtwodf<-as.data.frame(matrix(t(timerunningspeedtwoFR),ncol=1))
  colnames(timerunningspeedtwodf)<-c("response")
  timerunningspeedtwodf$frequency<-rep(freq,rep(length(itime)*length(irsbins),nf))
  timerunningspeedtwodf$runningspeed<-rep(irsbins,length(itime)*nf)
  timerunningspeedtwodf$time<-rep(rep(itime,rep(length(irsbins),length(itime))),nf)
  keep<-(round((1:nb/nb)/0.05,digits=3) %% 1) == 0
  keep[1]<-TRUE
  keep<-rep(rep(keep,rep(length(irsbins),length(keep))),nf)
  timerunningspeedtwodf<-timerunningspeedtwodf[keep,]
  timerunningspeedtwodf$time<-round(timerunningspeedtwodf$time)
  cl<-0.9*max(abs(timerunningspeedonedf$response))
  p<-ggplot(timerunningspeedtwodf,aes(x=runningspeed,y=frequency,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_raster()+scale_x_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,100,by=20),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    facet_wrap(~time,ncol=3)+ theme(aspect.ratio = 1)
  
  #session three x running speed tiled by frequency
  timerunningspeedthreedf<-as.data.frame(matrix(t(timerunningspeedthreeFR),ncol=1))
  colnames(timerunningspeedthreedf)<-c("response")
  timerunningspeedthreedf$frequency<-rep(freq,rep(length(itime)*length(irsbins),nf))
  timerunningspeedthreedf$runningspeed<-rep(irsbins,length(itime)*nf)
  timerunningspeedthreedf$time<-rep(rep(itime,rep(length(irsbins),length(itime))),nf)
  cl<-0.75*max(abs(timerunningspeedonedf$response))
  p<-ggplot(timerunningspeedthreedf,aes(x=time,y=runningspeed,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_raster()+scale_x_continuous(breaks=seq(0,600,120),labels=seq(0,10,by=2),expand=c(0,0))+
    scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.2),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    facet_wrap(~frequency,nrow=5) + theme(aspect.ratio = 1)
  
  #session three x running speed tiled by time
  timerunningspeedthreedf<-as.data.frame(matrix(t(timerunningspeedthreeFR),ncol=1))
  colnames(timerunningspeedthreedf)<-c("response")
  timerunningspeedthreedf$frequency<-rep(freq,rep(length(itime)*length(irsbins),nf))
  timerunningspeedthreedf$runningspeed<-rep(irsbins,length(itime)*nf)
  timerunningspeedthreedf$time<-rep(rep(itime,rep(length(irsbins),length(itime))),nf)
  keep<-(round((1:nb/nb)/0.05,digits=3) %% 1) == 0
  keep[1]<-TRUE
  keep<-rep(rep(keep,rep(length(irsbins),length(keep))),nf)
  timerunningspeedthreedf<-timerunningspeedthreedf[keep,]
  timerunningspeedthreedf$time<-round(timerunningspeedthreedf$time)
  cl<-0.9*max(abs(timerunningspeedonedf$response))
  p<-ggplot(timerunningspeedthreedf,aes(x=runningspeed,y=frequency,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_raster()+scale_x_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,100,by=20),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    facet_wrap(~time,ncol=3)+ theme(aspect.ratio = 1)
  
  
  #session one-two x running speed tiled by frequency
  timerunningspeedonetwodf<-as.data.frame(matrix(t(timerunningspeedoneFR-timerunningspeedtwoFR),ncol=1))
  colnames(timerunningspeedonetwodf)<-c("response")
  timerunningspeedonetwodf$frequency<-rep(freq,rep(length(itime)*length(irsbins),nf))
  timerunningspeedonetwodf$runningspeed<-rep(irsbins,length(itime)*nf)
  timerunningspeedonetwodf$time<-rep(rep(itime,rep(length(irsbins),length(itime))),nf)
  cl<-0.75*max(abs(timerunningspeedonetwodf$response))
  p<-ggplot(timerunningspeedonetwodf,aes(x=time,y=runningspeed,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_raster()+scale_x_continuous(breaks=seq(0,600,120),labels=seq(0,10,by=2),expand=c(0,0))+
    scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    facet_wrap(~frequency,nrow=5) + theme(aspect.ratio = 1)
  
  #session one-two x running speed tiled by time
  timerunningspeedonetwodf<-as.data.frame(matrix(t(timerunningspeedoneFR-timerunningspeedtwoFR),ncol=1))
  colnames(timerunningspeedonetwodf)<-c("response")
  timerunningspeedonetwodf$frequency<-rep(freq,rep(length(itime)*length(irsbins),nf))
  timerunningspeedonetwodf$runningspeed<-rep(irsbins,length(itime)*nf)
  timerunningspeedonetwodf$time<-rep(rep(itime,rep(length(irsbins),length(itime))),nf)
  keep<-(round((1:nb/nb)/0.05,digits=3) %% 1) == 0
  keep[1]<-TRUE
  keep<-rep(rep(keep,rep(length(irsbins),length(keep))),nf)
  timerunningspeedonetwodf<-timerunningspeedonetwodf[keep,]
  timerunningspeedonetwodf$time<-round(timerunningspeedonetwodf$time)
  cl<-0.9*max(abs(timerunningspeedonetwodf$response))
  p<-ggplot(timerunningspeedonetwodf,aes(x=runningspeed,y=frequency,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_raster()+scale_x_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,100,by=20),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    facet_wrap(~time,ncol=3)+ theme(aspect.ratio = 1)
  
  #session one-three x running speed tiled by frequency
  timerunningspeedonethreedf<-as.data.frame(matrix(t(timerunningspeedoneFR-timerunningspeedthreeFR),ncol=1))
  colnames(timerunningspeedonethreedf)<-c("response")
  timerunningspeedonethreedf$frequency<-rep(freq,rep(length(itime)*length(irsbins),nf))
  timerunningspeedonethreedf$runningspeed<-rep(irsbins,length(itime)*nf)
  timerunningspeedonethreedf$time<-rep(rep(itime,rep(length(irsbins),length(itime))),nf)
  cl<-0.75*max(abs(timerunningspeedonetwodf$response))
  p<-ggplot(timerunningspeedonethreedf,aes(x=time,y=runningspeed,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_raster()+scale_x_continuous(breaks=seq(0,600,120),labels=seq(0,10,by=2),expand=c(0,0))+
    scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    facet_wrap(~frequency,nrow=5) + theme(aspect.ratio = 1)
  
  #session one-three x running speed tiled by time
  timerunningspeedonethreedf<-as.data.frame(matrix(t(timerunningspeedoneFR-timerunningspeedthreeFR),ncol=1))
  colnames(timerunningspeedonethreedf)<-c("response")
  timerunningspeedonethreedf$frequency<-rep(freq,rep(length(itime)*length(irsbins),nf))
  timerunningspeedonethreedf$runningspeed<-rep(irsbins,length(itime)*nf)
  timerunningspeedonethreedf$time<-rep(rep(itime,rep(length(irsbins),length(itime))),nf)
  keep<-(round((1:nb/nb)/0.05,digits=3) %% 1) == 0
  keep[1]<-TRUE
  keep<-rep(rep(keep,rep(length(irsbins),length(keep))),nf)
  timerunningspeedonethreedf<-timerunningspeedonethreedf[keep,]
  timerunningspeedonethreedf$time<-round(timerunningspeedonethreedf$time)
  cl<-0.9*max(abs(timerunningspeedonetwodf$response))
  p<-ggplot(timerunningspeedonethreedf,aes(x=runningspeed,y=frequency,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_raster()+scale_x_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,100,by=20),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    facet_wrap(~time,ncol=3)+ theme(aspect.ratio = 1)
  
  #session two-three x running speed tiled by frequency
  timerunningspeedtwothreedf<-as.data.frame(matrix(t(timerunningspeedtwoFR-timerunningspeedthreeFR),ncol=1))
  colnames(timerunningspeedtwothreedf)<-c("response")
  timerunningspeedtwothreedf$frequency<-rep(freq,rep(length(itime)*length(irsbins),nf))
  timerunningspeedtwothreedf$runningspeed<-rep(irsbins,length(itime)*nf)
  timerunningspeedtwothreedf$time<-rep(rep(itime,rep(length(irsbins),length(itime))),nf)
  cl<-0.75*max(abs(timerunningspeedonetwodf$response))
  p<-ggplot(timerunningspeedtwothreedf,aes(x=time,y=runningspeed,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_raster()+scale_x_continuous(breaks=seq(0,600,120),labels=seq(0,10,by=2),expand=c(0,0))+
    scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    facet_wrap(~frequency,nrow=5) + theme(aspect.ratio = 1)
  
  #session two-three x running speed tiled by time
  timerunningspeedtwothreedf<-as.data.frame(matrix(t(timerunningspeedtwoFR-timerunningspeedthreeFR),ncol=1))
  colnames(timerunningspeedtwothreedf)<-c("response")
  timerunningspeedtwothreedf$frequency<-rep(freq,rep(length(itime)*length(irsbins),nf))
  timerunningspeedtwothreedf$runningspeed<-rep(irsbins,length(itime)*nf)
  timerunningspeedtwothreedf$time<-rep(rep(itime,rep(length(irsbins),length(itime))),nf)
  keep<-(round((1:nb/nb)/0.05,digits=3) %% 1) == 0
  keep[1]<-TRUE
  keep<-rep(rep(keep,rep(length(irsbins),length(keep))),nf)
  timerunningspeedtwothreedf<-timerunningspeedtwothreedf[keep,]
  timerunningspeedtwothreedf$time<-round(timerunningspeedtwothreedf$time)
  cl<-0.9*max(abs(timerunningspeedonetwodf$response))
  p<-ggplot(timerunningspeedtwothreedf,aes(x=runningspeed,y=frequency,fill=response,z=response))
  myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
  p+geom_raster()+scale_x_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,100,by=20),expand=c(0,0))+
    scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
    #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
    facet_wrap(~time,ncol=3)+ theme(aspect.ratio = 1)
  
  
# }
