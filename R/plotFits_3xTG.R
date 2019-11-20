library(mgcv)
library(ggplot2)
library(MASS)
library(itsadug)
library(RColorBrewer)
library(wesanderson)
library(directlabels)
library(scales)
library(akima)

setwd("E:/For Brian/GAMS_full/")

#load my data
cwd<-getwd()
MYDATA<-readRDS("MYDATAPOW_full.rds")
# MYDATA$Mouse<-ordered(MYDATA$Mouse)
# MYDATA<-MYDATA[MYDATA$Mouse>3,]

#Set model prediction vars
numdays<-length(unique(MYDATA$Day))
days<-seq(1,numdays)
nb<-45
nf<-50
freq<-10^seq(log10(2),log10(100),length.out=nf)
mint<-0
maxt<-600
tbins<-seq(mint,maxt,length.out=nb)
minrs<-0
maxrs<-0.9
rsbins<-seq(minrs,maxrs,length.out=nb)

#create new data
dgrid<-expand.grid(rsbins,tbins,unique(MYDATA$Session))
colnames(dgrid)<-c("runningspeed","time","session")
newdat<-with(MYDATA,data.frame(Running_Speed=dgrid$runningspeed,Time=dgrid$time,
                               Session=dgrid$session))

#define indices and remove unused predictions
pridx<-which(!duplicated(newdat[,c('Running_Speed')]))
ptidx<-which(!duplicated(newdat[,c('Time','Session')]))
used<-c(pridx,ptidx)
newdat<-newdat[sort(unique(used)),]
newdat$Day<-MYDATA$Day[1]
newdat$Trial<-MYDATA$Trial[1]
newdat$Mouse<-MYDATA$Mouse[1]

#define new indices
pi1idx<-which(!duplicated(newdat[,c('Session')]))[1]
pi2idx<-which(!duplicated(newdat[,c('Session')]))[2]
pi3idx<-which(!duplicated(newdat[,c('Session')]))[3]
pridx<-which(!duplicated(newdat[,c('Running_Speed')]))
ptidx<-which(!duplicated(newdat[,c('Time','Session')]))
pt1idx<-ptidx[newdat$Session[ptidx]=="1"]
pt2idx<-ptidx[newdat$Session[ptidx]=="2"]
pt3idx<-ptidx[newdat$Session[ptidx]=="3"]

#preallocate plot data for iterms
lags<-30
acf<-matrix(0,nrow=nf,ncol=lags)
acfrho<-acf
rhofit<-matrix(0,nrow=nf,ncol=1)
rsfit<-matrix(0,nrow=nb*nf,ncol=3)
timefit<-matrix(0,nrow=nb*nf*3,ncol=4)
AlzRsq<-matrix(0,nrow=1,ncol=nf)
HybRsq<-matrix(0,nrow=1,ncol=nf)

#set known values
rsfit[,2]<-rep(freq,each=nb) #frequency
rsfit[,3]<-rep(rsbins,times=nf) #runnings speed
timefit[,2]<-rep(unique(MYDATA$Session),each=nb*nf) #session
timefit[,3]<-rep(rep(freq,each=nb),times=3) #frequency
timefit[,4]<-rep(tbins,times=3*nf) #time

#prediction using "lpmatrix"
ns<-1500
Alzisim<-matrix(0,nrow = 3*nf,ncol = ns)
Alzrsim<-matrix(0,nrow = nf*nb,ncol = ns)
Alztimesim<-matrix(0,nrow = 3*nf*nb,ncol = ns)
Hybisim<-matrix(0,nrow = 3*nf,ncol = ns)
Hybrsim<-matrix(0,nrow = nf*nb,ncol = ns)
Hybtimesim<-matrix(0,nrow = 3*nf*nb,ncol = ns)
fnum<-1:nf

for(f in 1:nf)
{

  Alzname<-paste0("Alz_f",fnum[f],"_POW.rds")
  Alzmod<-readRDS(Alzname)
  Hybname<-paste0("Hyb_f",fnum[f],"_POW.rds")
  Hybmod<-readRDS(Hybname)

  AlzRsq[f] <- 1 - Alzmod$deviance/Alzmod$null.deviance
  HybRsq[f] <- 1 - Hybmod$deviance/Hybmod$null.deviance

  #get mean coefficient values and variance-covariance matrix
  Alzcoefs<-coef(Alzmod)
  Alzvc<-vcov(Alzmod)
  Hybcoefs<-coef(Hybmod)
  Hybvc<-vcov(Hybmod)
  #draw random coefficient values from the model
  Alzsim<-mvrnorm(ns,mu=Alzcoefs,Sigma=Alzvc)
  Hybsim<-mvrnorm(ns,mu=Hybcoefs,Sigma=Hybvc)

  #get predictions
  newdat$Mouse <- MYDATA$Mouse[MYDATA$Strain=="3xTG"][1] 
  newdat$Day <- MYDATA$Day[MYDATA$Strain=="3xTG"][1] 
  newdat$Trial <- MYDATA$Trial[MYDATA$Strain=="3xTG"][1] 
  Alzpred<-predict.bam(Alzmod, newdat, type = "lpmatrix", newdata.guaranteed = TRUE)
  
  newdat$Mouse <- MYDATA$Mouse[MYDATA$Strain=="hyb"][1] 
  newdat$Day <- MYDATA$Day[MYDATA$Strain=="hyb"][1] 
  newdat$Trial <- MYDATA$Trial[MYDATA$Strain=="hyb"][1] 
  Hybpred<-predict.bam(Hybmod, newdat, type = "lpmatrix", newdata.guaranteed = TRUE)

  rm(Alzcoefs,Alzvc,Hybcoefs,Hybvc)
  gc()

  #define indices
  ridx<-((f-1)*nb+1):(f*nb)
  t1idx<-((f-1)*nb+1):(f*nb)
  t2idx<-((f-1)*nb+1+nb*nf):(f*nb+nb*nf)
  t3idx<-((f-1)*nb+1+2*nb*nf):(f*nb+2*nb*nf)
  
  #intercepts
  Alzisim[f,] <- Alzpred[pi1idx, 1] %*% t(Alzsim[, 1]) #session 1
  Alzisim[nf+f,] <- Alzpred[pi2idx, 2] %*% t(Alzsim[, 2]) + Alzisim[f,] #session 2
  Alzisim[2*nf+f,] <- Alzpred[pi3idx, 3] %*% t(Alzsim[, 3]) + Alzisim[f,] #session 3
  Hybisim[f,] <- Hybpred[pi1idx, 1] %*% t(Hybsim[, 1]) #session 1
  Hybisim[nf+f,] <- Hybpred[pi2idx, 2] %*% t(Hybsim[, 2]) + Hybisim[f,] #session 2
  Hybisim[2*nf+f,] <- Hybpred[pi3idx, 3] %*% t(Hybsim[, 3]) + Hybisim[f,] #session 3

  #runningspeed
  want<-grep("s\\(Running\\_Speed\\)",names(Alzmod$coefficients))
  Alzrsim[ridx,]<-Alzpred[pridx,want] %*% t(Alzsim[,want])
  want<-grep("s\\(Running\\_Speed\\)",names(Hybmod$coefficients))
  Hybrsim[ridx,]<-Hybpred[pridx,want] %*% t(Hybsim[,want])

  #time by session
  want1<-grep("s\\(Time\\)\\:Session1",names(Alzmod$coefficients))
  want2<-grep("s\\(Time\\)\\:Session2",names(Alzmod$coefficients))
  want3<-grep("s\\(Time\\)\\:Session3",names(Alzmod$coefficients))
  Alztimesim[t1idx,] <- Alzpred[pt1idx, want1] %*% t(Alzsim[, want1]) #session 1
  Alztimesim[t2idx,] <- Alzpred[pt2idx, want2] %*% t(Alzsim[, want2]) #session 2
  Alztimesim[t3idx,] <- Alzpred[pt3idx, want3] %*% t(Alzsim[, want3]) #session 3
  want1<-grep("s\\(Time\\)\\:Session1",names(Hybmod$coefficients))
  want2<-grep("s\\(Time\\)\\:Session2",names(Hybmod$coefficients))
  want3<-grep("s\\(Time\\)\\:Session3",names(Hybmod$coefficients))
  Hybtimesim[t1idx,] <- Hybpred[pt1idx, want1] %*% t(Hybsim[, want1]) #session 1
  Hybtimesim[t2idx,] <- Hybpred[pt2idx, want2] %*% t(Hybsim[, want2]) #session 2
  Hybtimesim[t3idx,] <- Hybpred[pt3idx, want3] %*% t(Hybsim[, want3]) #session 3

  rm(Alzpred,Alzsim,Hybpred,Hybsim)
  gc()
  
  print(paste("Completed: ",f," of ",nf),quote=FALSE)
  
}
rm(newdat)
gc()

#define indices for plot terms
isimidx<-1:nrow(Alzisim) #intercepts
rsimidx<-(max(isimidx)+1):(max(isimidx)+nrow(Alzrsim)) #running speed
timesimidx<-(max(rsimidx)+1):(max(rsimidx)+nrow(Alztimesim)) #time by session

#get simultaneous CI's for main terms
Alzsims<-rbind(Alzisim,Alzrsim)
rm(Alzisim,Alzrsim)
gc()
Alzsims<-rbind(Alzsims,Alztimesim)
rm(Alztimesim)
gc()
Hybsims<-rbind(Hybisim,Hybrsim)
rm(Hybisim,Hybrsim)
gc()
Hybsims<-rbind(Hybsims,Hybtimesim)
rm(Hybtimesim)
gc()
dAHsims <- Alzsims - Hybsims

Alzci<-matrix(0,nrow = nrow(Alzsims),ncol = 3) #matrix to hold results
Alzci[,2]<-rowMeans(Alzsims)
Alzsims_sd<-rep(0,times=nrow(Alzsims))
for(i in 1:nrow(Alzsims)){
  Alzsims_sd[i]<-sd(Alzsims[i,])
  Alzsims[i,]<-abs(Alzsims[i,]-Alzci[i,2])/Alzsims_sd[i]
}
gc()
maxima<-rep(0,times=ns)
for(i in 1:ns){maxima[i]<-max(Alzsims[,i])}
m <- quantile(maxima,probs=0.95,type=8)
Alzci[,1]<--m*Alzsims_sd
Alzci[,3]<-m*Alzsims_sd
rm(Alzsims)
gc()

Hybci<-matrix(0,nrow = nrow(Hybsims),ncol = 3) #matrix to hold results
Hybci[,2]<-rowMeans(Hybsims)
Hybsims_sd<-rep(0,times=nrow(Hybsims))
for(i in 1:nrow(Hybsims)){
  Hybsims_sd[i]<-sd(Hybsims[i,])
  Hybsims[i,]<-abs(Hybsims[i,]-Hybci[i,2])/Hybsims_sd[i]
}
gc()
maxima<-rep(0,times=ns)
for(i in 1:ns){maxima[i]<-max(Hybsims[,i])}
m <- quantile(maxima,probs=0.95,type=8)
Hybci[,1]<--m*Hybsims_sd
Hybci[,3]<-m*Hybsims_sd
rm(Hybsims)
gc()

#get simultaneous CI's for between group difference terms
dAHci<-matrix(0,nrow = nrow(dAHsims),ncol = 3) #matrix to hold results
dAHci[,2]<-rowMeans(dAHsims)
dAHsims_sd<-rep(0,times=nrow(dAHsims))
for(i in 1:nrow(dAHsims)){
  dAHsims_sd[i]<-sd(dAHsims[i,])
  dAHsims[i,]<-abs(dAHsims[i,]-dAHci[i,2])/dAHsims_sd[i]
}
gc()
maxima<-rep(0,times=ns)
for(i in 1:ns){maxima[i]<-max(dAHsims[,i])}
m <- quantile(maxima,probs=0.95,type=8)
dAHci[,1]<--m*dAHsims_sd
dAHci[,3]<-m*dAHsims_sd

rm(dAHsims,maxima)
gc()

# #get between session difference terms
# Alzisim <- exp(Alzisim)
# Alztimesim <- exp(Alztimesim)
# dAlzisim<-matrix(0,nrow = nrow(Alzisim),ncol = ncol(Alzisim))
# dAlzisim[1:nf,]<-Alzisim[1:nf,]-Alzisim[(1+nf):(2*nf),]
# dAlzisim[(1+nf):(2*nf),]<-Alzisim[1:nf,]-Alzisim[(1+2*nf):(3*nf),]
# dAlzisim[(1+2*nf):(3*nf),]<-Alzisim[(1+nf):(2*nf),]-Alzisim[(1+2*nf):(3*nf),]
# rm(Alzisim)
# gc()
# dAlztimesim<-matrix(0,nrow = nrow(Alztimesim),ncol = ncol(Alztimesim))
# dAlztimesim[1:(nf*nb),]<-Alztimesim[1:(nf*nb),]-Alztimesim[(1+nb*nf):(2*nb*nf),]
# dAlztimesim[(1+nb*nf):(2*nb*nf),]<-Alztimesim[1:(nf*nb),]-Alztimesim[(1+2*nb*nf):(3*nb*nf),]
# dAlztimesim[(1+2*nb*nf):(3*nb*nf),]<-Alztimesim[(1+nb*nf):(2*nb*nf),]-Alztimesim[(1+2*nb*nf):(3*nb*nf),]
# rm(Alztimesim)
# gc()
# Hybisim <- exp(Hybisim)
# Hybtimesim <- exp(Hybtimesim)
# dHybisim<-matrix(0,nrow = nrow(Hybisim),ncol = ncol(Hybisim))
# dHybisim[1:nf,]<-Hybisim[1:nf,]-Hybisim[(1+nf):(2*nf),]
# dHybisim[(1+nf):(2*nf),]<-Hybisim[1:nf,]-Hybisim[(1+2*nf):(3*nf),]
# dHybisim[(1+2*nf):(3*nf),]<-Hybisim[(1+nf):(2*nf),]-Hybisim[(1+2*nf):(3*nf),]
# rm(Hybisim)
# gc()
# dHybtimesim<-matrix(0,nrow = nrow(Hybtimesim),ncol = ncol(Hybtimesim))
# dHybtimesim[1:(nf*nb),]<-Hybtimesim[1:(nf*nb),]-Hybtimesim[(1+nb*nf):(2*nb*nf),]
# dHybtimesim[(1+nb*nf):(2*nb*nf),]<-Hybtimesim[1:(nf*nb),]-Hybtimesim[(1+2*nb*nf):(3*nb*nf),]
# dHybtimesim[(1+2*nb*nf):(3*nb*nf),]<-Hybtimesim[(1+nb*nf):(2*nb*nf),]-Hybtimesim[(1+2*nb*nf):(3*nb*nf),]
# rm(Hybtimesim)
# gc()
# 
# #get simultaneous CI's for difference terms
# Alzsims<-rbind(dAlzisim,Alzrsim)
# rm(dAlzisim,Alzrsim)
# gc()
# Alzsims<-rbind(Alzsims,dAlztimesim)
# rm(dAlztimesim)
# gc()
# Hybsims<-rbind(dHybisim,Hybrsim)
# rm(dHybisim,Hybrsim)
# gc()
# Hybsims<-rbind(Hybsims,dHybtimesim)
# rm(dHybtimesim)
# gc()
# 
# dAlzci<-matrix(0,nrow = nrow(Alzsims),ncol = 3) #matrix to hold results
# dAlzci[,2]<-rowMeans(Alzsims)
# Alzsims_sd<-rep(0,times=nrow(Alzsims))
# for(i in 1:nrow(Alzsims)){
#   Alzsims_sd[i]<-sd(Alzsims[i,])
#   Alzsims[i,]<-abs(Alzsims[i,]-dAlzci[i,2])/Alzsims_sd[i]
# }
# gc()
# maxima<-rep(0,times=ns)
# for(i in 1:ns){maxima[i]<-max(Alzsims[,i])}
# m <- quantile(maxima,probs=0.95,type=8)
# dAlzci[,1]<--m*Alzsims_sd
# dAlzci[,3]<-m*Alzsims_sd
# rm(Alzsims,maxima)
# gc()
# 
# dHybci<-matrix(0,nrow = nrow(Hybsims),ncol = 3) #matrix to hold results
# dHybci[,2]<-rowMeans(Hybsims)
# Hybsims_sd<-rep(0,times=nrow(Hybsims))
# for(i in 1:nrow(Hybsims)){
#   Hybsims_sd[i]<-sd(Hybsims[i,])
#   Hybsims[i,]<-abs(Hybsims[i,]-dHybci[i,2])/Hybsims_sd[i]
# }
# gc()
# maxima<-rep(0,times=ns)
# for(i in 1:ns){maxima[i]<-max(Hybsims[,i])}
# m <- quantile(maxima,probs=0.95,type=8)
# dHybci[,1]<--m*Hybsims_sd
# dHybci[,3]<-m*Hybsims_sd
# rm(Hybsims,maxima)
# gc()

#r-squared plot
rsqdf<-data.frame(Frequency=rep(freq,times=2),
                  Rsquared=c(as.numeric(HybRsq),as.numeric(AlzRsq)),
                  Strain = c(rep("Hyb",times=nf),rep("Alz",times=nf)))
p<-ggplot(rsqdf,aes(x=Frequency,y=Rsquared,group=Strain,color=Strain)) 
p + geom_line()+geom_path(size=1)+
    scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,0.2,by=0.025),name="R^2")+
    theme(axis.line.x = element_line(color="black", size = 1),
          axis.line.y = element_line(color="black", size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank())

#rho plot
rhodf<-data.frame(Frequency=rep(freq,times=2),
                  rho=c(acf[,1],acfrho[,1]),
                  model=rep(c("Normal","AR1"),each=nf))
p<-ggplot(rhodf,aes(x=Frequency,y=rho,group=model))
p+geom_line(size=1,aes(colour=model))+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,1,by=0.25),labels=seq(0,1,by=0.25),limits=c(-0.05,1),expand=c(0,0))+
  # geom_hline(aes(yintercept=0.05))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.ticks = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA))

#acf plot
acfdf<-as.data.frame(c(as.vector(acf),as.vector(acfrho)))
colnames(acfdf)<-c("acf")
acfdf$lags<-rep(rep(1:lags,each=nf),2)
acfdf$frequency<-rep(freq,times=2*lags)
acfdf$rho<-rep(c("normal","AR1"),each=nf*lags)

#interpolate
lbins<-1:30
fbins<-seq(2,100,length=250)
idxnorm<-acfdf$rho=="normal"
idxrho<-acfdf$rho=="AR1"
acfinterp<-interp(x=acfdf$lags[idxnorm],y=acfdf$frequency[idxnorm],acfdf$acf[idxnorm],xo=lbins,yo=fbins)
acfrhointerp<-interp(x=acfdf$lags[idxrho],y=acfdf$frequency[idxrho],acfdf$acf[idxrho],xo=lbins,yo=fbins)
acfnormdf <- data.frame(expand.grid(lags = acfinterp$x, frequency = acfinterp$y), acf = c(acfinterp$z))
acfrhodf <- data.frame(expand.grid(lags = acfrhointerp$x, frequency = acfrhointerp$y), acf = c(acfrhointerp$z))
acfnormdf$rho<-rep(factor("normal"),times=length(acfnormdf[,1]))
acfrhodf$rho<-rep(factor("AR1"),times=length(acfrhodf[,1]))
acfabsrhodf<-acfrhodf
acfabsrhodf$acf<-abs(acfnormdf$acf)-abs(acfrhodf$acf)
acfabsrhodf$rho<-"|Normal|-|AR1|"
acfdf<-rbind(acfnormdf,acfrhodf,acfabsrhodf)


#display plot
myPalette <- colorRampPalette(rev(brewer.pal(10, "PRGn")), space="Lab")
cl<-1
ggplot(acfdf,aes(x=lags,y=frequency,fill=acf,z=acf))+ 
  geom_raster()+scale_x_continuous(breaks=c(1,10,20,30),labels=c(1,10,20,30),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",alpha=0.3,binwidth=0.05,size=0.05,show.legend = TRUE)+
  facet_wrap(~rho,nrow=1)+theme(aspect.ratio = 1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA,size=1))

#intercepts by session
#run for "lpmatrix"
ci <- Alzci
ifit<-as.data.frame(ci[isimidx,])
colnames(ifit)<-c("low","fit","high")
ifit$low<-ifit$low+ifit$fit
ifit$high<-ifit$high+ifit$fit
ifit$session<-rep(rep(c("1","2","3")),each=nf)
ifit$frequency<-rep(freq,times=3) 
ifit$fit <- exp(ifit$fit)
ifit$low <- exp(ifit$low)
ifit$high <- exp(ifit$high)

p<-ggplot(ifit,aes(x=frequency,y=fit))
p+geom_ribbon(aes(ymin=low,ymax=high))+
  geom_line(size=1)+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  facet_wrap(~session,nrow=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.ticks = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA))

ggsave("ampsession_Hyb.eps",device="eps",width=14,height=14,units="in")

#running speed plot
#run for "lpmatrix"
ci <- Hybci
rbins<-seq(min(rsbins),max(rsbins),length=250)
fbins<-seq(min(freq),max(freq),length=250)
fitinterp<-akima::interp(x=rsfit[,3],y=rsfit[,2],z=ci[rsimidx,2],xo=rbins,yo=fbins)
highinterp<-akima::interp(x=rsfit[,3],y=rsfit[,2],z=ci[rsimidx,2]+ci[rsimidx,3],xo=rbins,yo=fbins)
lowinterp<-akima::interp(x=rsfit[,3],y=rsfit[,2],z=ci[rsimidx,2]+ci[rsimidx,1],xo=rbins,yo=fbins)
rfit <- data.frame(expand.grid(runningspeed = fitinterp$x, frequency = fitinterp$y), 
                   fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z))
rfit$fit <- exp(rfit$fit)
rfit$low <- exp(rfit$low)
rfit$high <- exp(rfit$high)

#display plot
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
cl<-0.5
ggplot(rfit,aes(x=runningspeed,y=frequency,fill=fit,z=fit))+
  geom_raster()+scale_x_continuous(breaks=seq(round(minrs),maxrs,by=0.2),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(1-cl,1+cl),oob=squish)+
  stat_contour(mapping=aes(x=runningspeed,y=frequency,z=high),colour="blue",breaks=c(1))+
  stat_contour(mapping=aes(x=runningspeed,y=frequency,z=low),colour="red",breaks=c(1))+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("speed_Alz.eps",device="eps",width=5,height=5,units="in")


#time by session
#run for "lpmatrix"
ci <- Hybci
tfit<-as.data.frame(ci[timesimidx,])
colnames(tfit)<-c("low","fit","high")
tfit$low<-tfit$low+tfit$fit
tfit$high<-tfit$high+tfit$fit
tfit$session<-rep(rep(c("1","2","3")),each=nb*nf)
tfit$frequency<-rep(rep(freq,each=nb),times=3) 
tfit$time<-rep(seq(mint,maxt,length.out=nb),times=3*nf)

#interpolate
tibins<-seq(mint,maxt,length=250)
fbins<-seq(min(freq),max(freq),length=250)
idxone<-tfit$session == "1"
idxtwo<-tfit$session == "2"
idxthree<-tfit$session == "3"
fitsone<-akima::interp(x=tfit$time[idxone],y=tfit$frequency[idxone],tfit$fit[idxone],xo=tibins,yo=fbins)
fitstwo<-akima::interp(x=tfit$time[idxtwo],y=tfit$frequency[idxtwo],tfit$fit[idxtwo],xo=tibins,yo=fbins)
fitsthree<-akima::interp(x=tfit$time[idxthree],y=tfit$frequency[idxthree],tfit$fit[idxthree],xo=tibins,yo=fbins)
highone<-akima::interp(x=tfit$time[idxone],y=tfit$frequency[idxone],tfit$high[idxone],xo=tibins,yo=fbins)
hightwo<-akima::interp(x=tfit$time[idxtwo],y=tfit$frequency[idxtwo],tfit$high[idxtwo],xo=tibins,yo=fbins)
highthree<-akima::interp(x=tfit$time[idxthree],y=tfit$frequency[idxthree],tfit$high[idxthree],xo=tibins,yo=fbins)
lowone<-akima::interp(x=tfit$time[idxone],y=tfit$frequency[idxone],tfit$low[idxone],xo=tibins,yo=fbins)
lowtwo<-akima::interp(x=tfit$time[idxtwo],y=tfit$frequency[idxtwo],tfit$low[idxtwo],xo=tibins,yo=fbins)
lowthree<-akima::interp(x=tfit$time[idxthree],y=tfit$frequency[idxthree],tfit$low[idxthree],xo=tibins,yo=fbins)
tfitone <- data.frame(expand.grid(time = fitsone$x, frequency = fitsone$y), 
                      fit = c(fitsone$z), high = c(highone$z), low = c(lowone$z))
tfittwo <- data.frame(expand.grid(time = fitstwo$x, frequency = fitstwo$y), 
                      fit = c(fitstwo$z), high = c(hightwo$z), low = c(lowtwo$z))
tfitthree <- data.frame(expand.grid(time = fitsthree$x, frequency = fitsthree$y), 
                        fit = c(fitsthree$z), high = c(highthree$z), low = c(lowthree$z))
tfitone$session<-rep(factor("1"),times=length(tfitone[,1]))
tfittwo$session<-rep(factor("2"),times=length(tfittwo[,1]))
tfitthree$session<-rep(factor("3"),times=length(tfitthree[,1]))
tfit<-rbind(tfitone,tfittwo,tfitthree)
tfit$fit <- exp(tfit$fit)
tfit$low <- exp(tfit$low)
tfit$high <- exp(tfit$high)

#display plot
cl<-0.5
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
ggplot(tfit,aes(x=time,y=frequency,fill=fit,z=fit))+ 
  geom_raster()+scale_x_continuous(breaks=seq(0,600,by=60),labels=seq(0,10,by=1),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(1-cl,1+cl),oob=squish)+
  #stat_contour(colour="black",alpha=0.3,binwidth=0.05,size=0.05,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=frequency,z=high),colour="blue",breaks=c(1),alpha=1)+
  stat_contour(mapping=aes(x=time,y=frequency,z=low),colour="red",breaks=c(1),alpha=1)+
  facet_wrap(~session,nrow=1)+theme(aspect.ratio = 1)+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("timesession_Hyb.eps",device="eps",width=14,height=14,units="in")

#differences time by session

#run for "lpmatrix"
dci <- dHybci
tfit<-as.data.frame(dci[timesimidx,])
colnames(tfit)<-c("low","fit","high")
tfit$low<-tfit$low+tfit$fit
tfit$high<-tfit$high+tfit$fit
tfit$session<-rep(rep(c("1","2","3")),each=nb*nf)
tfit$frequency<-rep(rep(freq,each=nb),times=3) 
tfit$time<-rep(seq(mint,maxt,length.out=nb),times=3*nf)

#interpolate for "lpmatrix"
tibins<-seq(tbins[1],tbins[length(tbins)],length=250)
fbins<-seq(min(freq),max(freq),length=250)
idxone<-tfit$session == "1"
idxtwo<-tfit$session == "2"
idxthree<-tfit$session == "3"
fitsone<-akima::interp(x=tfit$time[idxone],y=tfit$frequency[idxone],tfit$fit[idxone],xo=tibins,yo=fbins)
fitstwo<-akima::interp(x=tfit$time[idxtwo],y=tfit$frequency[idxtwo],tfit$fit[idxtwo],xo=tibins,yo=fbins)
fitsthree<-akima::interp(x=tfit$time[idxthree],y=tfit$frequency[idxthree],tfit$fit[idxthree],xo=tibins,yo=fbins)
highone<-akima::interp(x=tfit$time[idxone],y=tfit$frequency[idxone],tfit$high[idxone],xo=tibins,yo=fbins)
hightwo<-akima::interp(x=tfit$time[idxtwo],y=tfit$frequency[idxtwo],tfit$high[idxtwo],xo=tibins,yo=fbins)
highthree<-akima::interp(x=tfit$time[idxthree],y=tfit$frequency[idxthree],tfit$high[idxthree],xo=tibins,yo=fbins)
lowone<-akima::interp(x=tfit$time[idxone],y=tfit$frequency[idxone],tfit$low[idxone],xo=tibins,yo=fbins)
lowtwo<-akima::interp(x=tfit$time[idxtwo],y=tfit$frequency[idxtwo],tfit$low[idxtwo],xo=tibins,yo=fbins)
lowthree<-akima::interp(x=tfit$time[idxthree],y=tfit$frequency[idxthree],tfit$low[idxthree],xo=tibins,yo=fbins)
tfitone <- data.frame(expand.grid(time = fitsone$x, frequency = fitsone$y), 
                      fit = c(fitsone$z), high = c(highone$z), low = c(lowone$z))
tfittwo <- data.frame(expand.grid(time = fitstwo$x, frequency = fitstwo$y), 
                      fit = c(fitstwo$z), high = c(hightwo$z), low = c(lowtwo$z))
tfitthree <- data.frame(expand.grid(time = fitsthree$x, frequency = fitsthree$y), 
                        fit = c(fitsthree$z), high = c(highthree$z), low = c(lowthree$z))
tfitone$session<-rep(factor("1"),times=length(tfitone[,1]))
tfittwo$session<-rep(factor("2"),times=length(tfittwo[,1]))
tfitthree$session<-rep(factor("3"),times=length(tfitthree[,1]))
tfit<-rbind(tfitone,tfittwo,tfitthree)

#display plot
tfit$dsessions<-as.factor(c(rep("1-2",times=250*250),rep("1-3",times=250*250),rep("2-3",times=250*250)))
cl<-0.3
p<-ggplot(tfit,aes(x=time,y=frequency,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster()+scale_x_continuous(breaks=seq(0,600,by=60),labels=seq(0,10,by=1),expand=c(0,0))+
  scale_y_continuous(breaks=seq(10,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",alpha=0.3,binwidth=0.05,size=0.05,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=frequency,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=frequency,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~dsessions,nrow=1)+theme(aspect.ratio = 1)+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dtimesession_Hyb.eps",device="eps",width=14,height=14,units="in")  

#time/running speed interactions by frequency session 1
#run for lpmatrix
trsfit<-as.data.frame(ci[rtintsimidx,])
colnames(trsfit)<-c("low","fit","high")
trsfit$low<-trsfit$low+trsfit$fit
trsfit$high<-trsfit$high+trsfit$fit
trsfit$session<-rtintfit[,2]
trsfit$frequency<-rtintfit[,3]
trsfit$time<-rtintfit[,4]
trsfit$runningspeed<-rtintfit[,5]

#set NA's
th = 1;
b = 1;
MYDATA$Time<-MYDATA$Time-min(MYDATA$Time)
MYDATA$Time<-MYDATA$Time/max(MYDATA$Time)
trsfit$time<-trsfit$time-min(trsfit$time)
trsfit$time<-trsfit$time/max(trsfit$time)
require(gplots)
h1<-hist2d(MYDATA$Time[MYDATA$Session=="1"],MYDATA$Running_Speed[MYDATA$Session=="1"],nbins=nb/b,show=FALSE)
h2<-hist2d(MYDATA$Time[MYDATA$Session=="2"],MYDATA$Running_Speed[MYDATA$Session=="2"],nbins=nb/b,show=FALSE)
h3<-hist2d(MYDATA$Time[MYDATA$Session=="3"],MYDATA$Running_Speed[MYDATA$Session=="3"],nbins=nb/b,show=FALSE)
counts1<-matrix(0,nb,nb)
counts2<-matrix(0,nb,nb)
counts3<-matrix(0,nb,nb)
for(i in 1:30){
  counts1[2*i,]<-rep(h1$counts[i,],each=b)
  counts1[2*i-1,]<-rep(h1$counts[i,],each=b)
  counts2[2*i,]<-rep(h2$counts[i,],each=b)
  counts2[2*i-1,]<-rep(h2$counts[i,],each=b)
  counts3[2*i,]<-rep(h3$counts[i,],each=b)
  counts3[2*i-1,]<-rep(h3$counts[i,],each=b)
}
counts1<-rep(as.vector(t(counts1)),times=nf)
counts2<-rep(as.vector(t(counts2)),times=nf)
counts3<-rep(as.vector(t(counts3)),times=nf)
counts<-c(counts1,counts2,counts3)
trsfit$fit[counts<th]<-NA
trsfit$low[counts<th]<-NA
trsfit$high[counts<th]<-NA

#display plot
trsfit$frequency<-round(trsfit$frequency,digits=1)
trsfit$frequency<-as.ordered(trsfit$frequency)
for(f in 1:nf){
  levels(trsfit$frequency)[f]<-paste0(levels(trsfit$frequency)[f]," Hz")
}
cl<-0.7*max(abs(trsfit$fit),na.rm=TRUE)
p<-ggplot(trsfit[trsfit$session=="1",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish,na.value="white")+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("speedtimesession1_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time/running speed interactions by frequency session 2
p<-ggplot(trsfit[trsfit$session=="2",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish,na.value="white")+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("speedtimesession2_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time/running speed interactions by frequency session 3
p<-ggplot(trsfit[trsfit$session=="3",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish,na.value="white")+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("speedtimesession3_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time/running speed interactions by frequency session 1-2
trsfit<-as.data.frame(dci[rtintsimidx,])
colnames(trsfit)<-c("low","fit","high")
trsfit$low<-trsfit$low+trsfit$fit
trsfit$high<-trsfit$high+trsfit$fit
trsfit$session<-rep(c("1","2","3"),each=nb*nb*nf)
trsfit$frequency<-rep(rep(freq,each=nb*nb),times=3)
trsfit$time<-rep(rep(seq(mint,maxt,length.out=nb),each=nb),times=nf*3)
trsfit$runningspeed<-rep(rsbins,times=nb*nf*3)

#set NA's
th = 1;
b = 1;
MYDATA$Time<-MYDATA$Time-min(MYDATA$Time)
MYDATA$Time<-MYDATA$Time/max(MYDATA$Time)
trsfit$time<-trsfit$time-min(trsfit$time)
trsfit$time<-trsfit$time/max(trsfit$time)
require(gplots)
h1<-hist2d(MYDATA$Time[MYDATA$Session=="1"],MYDATA$Running_Speed[MYDATA$Session=="1"],nbins=nb/b,show=FALSE)
h2<-hist2d(MYDATA$Time[MYDATA$Session=="2"],MYDATA$Running_Speed[MYDATA$Session=="2"],nbins=nb/b,show=FALSE)
h3<-hist2d(MYDATA$Time[MYDATA$Session=="3"],MYDATA$Running_Speed[MYDATA$Session=="3"],nbins=nb/b,show=FALSE)
counts1<-matrix(0,nb,nb)
counts2<-matrix(0,nb,nb)
counts3<-matrix(0,nb,nb)
for(i in 1:30){
  counts1[2*i,]<-rep(h1$counts[i,],each=b)
  counts1[2*i-1,]<-rep(h1$counts[i,],each=b)
  counts2[2*i,]<-rep(h2$counts[i,],each=b)
  counts2[2*i-1,]<-rep(h2$counts[i,],each=b)
  counts3[2*i,]<-rep(h3$counts[i,],each=b)
  counts3[2*i-1,]<-rep(h3$counts[i,],each=b)
}
counts1<-rep(as.vector(t(counts1)),times=nf)
counts2<-rep(as.vector(t(counts2)),times=nf)
counts3<-rep(as.vector(t(counts3)),times=nf)
nas1<-rep(0,times=length(counts1))
nas2<-nas1
nas3<-nas1
nas1[counts1<th | counts2<th]<-NA
nas2[counts1<th | counts3<th]<-NA
nas3[counts2<th | counts3<th]<-NA
nas<-c(nas1,nas2,nas3)
trsfit$fit[is.na(nas)]<-NA
trsfit$low[is.na(nas)]<-NA
trsfit$high[is.na(nas)]<-NA

trsfit$frequency<-round(trsfit$frequency,digits=1)
trsfit$frequency<-as.ordered(trsfit$frequency)
for(f in 1:nf){
  levels(trsfit$frequency)[f]<-paste0(levels(trsfit$frequency)[f]," Hz")
}
cl<-0.8*max(abs(trsfit$fit),na.rm=TRUE)
p<-ggplot(trsfit[trsfit$session=="1",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish,na.value="white")+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="blue",breaks=c(0))+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="red",breaks=c(0))+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dspeedtimesession1_Hyb.eps",device="eps",width=5.5,height=11,units="in")

#time/running speed interactions by frequency session 1-3
p<-ggplot(trsfit[trsfit$session=="2",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish,na.value="white")+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="blue",breaks=c(0))+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="red",breaks=c(0))+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dspeedtimesession2_Hyb.eps",device="eps",width=5.5,height=11,units="in")

#time/running speed interactions by frequency session 2-3
p<-ggplot(trsfit[trsfit$session=="3",],aes(x=time,y=runningspeed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(round(minrs),round(maxrs),by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish,na.value="white")+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=high),colour="blue",breaks=c(0))+
  stat_contour(mapping=aes(x=time,y=runningspeed,z=low),colour="red",breaks=c(0))+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dspeedtimesession3_Hyb.eps",device="eps",width=5.5,height=11,units="in")

#running speed theta phase interaction tiled by frequency
#run for lpmatrix
rspfit<-as.data.frame(ci[rpintsimidx,])
colnames(rspfit)<-c("low","fit","high")
rspfit$low<-rspfit$low+rspfit$fit
rspfit$high<-rspfit$high+rspfit$fit
rspfit$frequency<-rpintfit[,2]
rspfit$phase<-rpintfit[,3]
rspfit$runningspeed<-rpintfit[,4]

# rspfit$phase<-rspfit$phase+pi
# rspfit$phase[rspfit$phase>pi]<-rspfit$phase[rspfit$phase>pi]-2*pi-median(diff(pbins))

#display plot
rspfit$frequency<-round(rspfit$frequency,digits=1)
rspfit$frequency<-as.ordered(rspfit$frequency)
for(f in 1:nf){
  levels(rspfit$frequency)[f]<-paste0(levels(rspfit$frequency)[f]," Hz")
}
cl<-max(abs(rspfit$fit))
p<-ggplot(rspfit,aes(x=runningspeed,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,1,by=0.5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi),labels=c(expression(paste("-",pi)),"0",expression(paste(pi))),
                     expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=runningspeed,y=phase,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=runningspeed,y=phase,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=5) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("speedphase_Alz.eps",device="eps",width=14,height=14,units="in")

#running speed theta phase interaction tiled by running speed
#run for lpmatrix
rspfit<-as.data.frame(ci[rpintsimidx,])
colnames(rspfit)<-c("low","fit","high")
rspfit$fit<-rspfit$fit + rep(ci[psimidx,2],each=nb)
rspfit$low<-rspfit$low+rspfit$fit
rspfit$high<-rspfit$high+rspfit$fit
rspfit$frequency<-rpintfit[,2]
rspfit$phase<-rpintfit[,3]
rspfit$runningspeed<-rpintfit[,4]
keep<-(round((1:nb/nb)/0.2,digits=3) %% 1) == 0
keep[1]<-TRUE
keep<-rep(keep,nb*nf)
rspfit<-rspfit[keep,]

#interpolate for each running speed
#run for "lpmatrix"
pbins<-seq(pbins[1],pbins[length(pbins)],length=250)
fbins<-seq(min(freq),max(freq),length=250)
rs<-unique(rspfit$runningspeed)
fitinterp<-interp(x=rspfit$phase[rspfit$runningspeed==rs[1]],
                  y=rspfit$frequency[rspfit$runningspeed==rs[1]],
                  rspfit$fit[rspfit$runningspeed==rs[1]],
                  xo=pbins,yo=fbins)
highinterp<-interp(x=rspfit$phase[rspfit$runningspeed==rs[1]],
                   y=rspfit$frequency[rspfit$runningspeed==rs[1]],
                   rspfit$high[rspfit$runningspeed==rs[1]],
                   xo=pbins,yo=fbins)
lowinterp<-interp(x=rspfit$phase[rspfit$runningspeed==rs[1]],
                  y=rspfit$frequency[rspfit$runningspeed==rs[1]],
                  rspfit$low[rspfit$runningspeed==rs[1]],
                  xo=pbins,yo=fbins)
rfit <- data.frame(expand.grid(phase = fitinterp$x, frequency = fitinterp$y), 
                   fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z))
for(r in 2:length(rs)){
  fitinterp<-interp(x=rspfit$phase[rspfit$runningspeed==rs[r]],
                    y=rspfit$frequency[rspfit$runningspeed==rs[r]],
                    rspfit$fit[rspfit$runningspeed==rs[r]],
                    xo=pbins,yo=fbins)
  highinterp<-interp(x=rspfit$phase[rspfit$runningspeed==rs[r]],
                     y=rspfit$frequency[rspfit$runningspeed==rs[r]],
                     rspfit$high[rspfit$runningspeed==rs[r]],
                     xo=pbins,yo=fbins)
  lowinterp<-interp(x=rspfit$phase[rspfit$runningspeed==rs[r]],
                    y=rspfit$frequency[rspfit$runningspeed==rs[r]],
                    rspfit$low[rspfit$runningspeed==rs[r]],
                    xo=pbins,yo=fbins)
  rfit <- rbind(rfit,data.frame(expand.grid(phase = fitinterp$x, frequency = fitinterp$y), 
                     fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z)))
  }
rfit$runningspeed<-round(rep(rs,each=250^2),digits = 2)

#display plot
cl<-max(abs(rfit$fit))
p<-ggplot(rfit,aes(x=phase,y=frequency,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster()+scale_x_continuous(breaks=seq(-pi,pi,pi/2),labels=c("-pi","-pi/2","0","pi/2","pi"),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,250,by=25),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  # stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  # stat_contour(mapping=aes(x=phase,y=frequency,z=high),colour="purple",breaks=c(0),alpha=1)+
  # stat_contour(mapping=aes(x=phase,y=frequency,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~runningspeed)


#time theta phase interaction session1
#run for lpmatrix
tpfit<-as.data.frame(ci[tpintsimidx,])
colnames(tpfit)<-c("low","fit","high")
tpfit$low<-tpfit$low+tpfit$fit
tpfit$high<-tpfit$high+tpfit$fit
tpfit$session<-as.factor(tpintfit[,2])
tpfit$frequency<-tpintfit[,3]
tpfit$phase<-tpintfit[,4]
tpfit$time<-tpintfit[,5]

# tpfit$diff<-tpfit$fit
# tpfit$diff[tpfit$session=="1"]<-tpfit$fit[tpfit$session=="1"]-tpfit$fit[tpfit$session=="2"]
# tpfit$diff[tpfit$session=="2"]<-tpfit$fit[tpfit$session=="1"]-tpfit$fit[tpfit$session=="3"]
# tpfit$diff[tpfit$session=="3"]<-tpfit$fit[tpfit$session=="2"]-tpfit$fit[tpfit$session=="3"]
# tpfit$fit<-tpfit$diff

# tpfit$phase<-tpfit$phase+pi
# tpfit$phase[tpfit$phase>pi]<-tpfit$phase[tpfit$phase>pi]-2*pi-median(diff(pbins))

#display plot
tpfit$frequency<-round(tpfit$frequency,digits=1)
tpfit$frequency<-as.ordered(tpfit$frequency)
for(f in 1:nf){
  levels(tpfit$frequency)[f]<-paste0(levels(tpfit$frequency)[f]," Hz")
}
cl<-max(abs(tpfit$fit))
p<-ggplot(tpfit[tpfit$session=="1",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi),
                     labels=c(expression(paste("-",pi)),"0",expression(paste(pi))),
                     expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("phasetimesession1_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time theta phase interaction session2
p<-ggplot(tpfit[tpfit$session=="2",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi),
                     labels=c(expression(paste("-",pi)),"0",expression(paste(pi))),
                     expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("phasetimesession2_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time theta phase interaction session3
p<-ggplot(tpfit[tpfit$session=="3",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi),
                     labels=c(expression(paste("-",pi)),"0",expression(paste(pi))),
                     expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("phasetimesession3_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time/theta phase interactions by frequency session 1-2
tpfit<-as.data.frame(dci[tpintsimidx,])
colnames(tpfit)<-c("low","fit","high")
tpfit$low<-tpfit$low+tpfit$fit
tpfit$high<-tpfit$high+tpfit$fit
tpfit$session<-as.factor(tpintfit[,2])
tpfit$frequency<-tpintfit[,3]
tpfit$phase<-tpintfit[,4]
tpfit$time<-tpintfit[,5]

# tpfit$phase<-tpfit$phase+pi
# tpfit$phase[tpfit$phase>pi]<-tpfit$phase[tpfit$phase>pi]-2*pi-median(diff(pbins))

tpfit$frequency<-round(tpfit$frequency,digits=1)
tpfit$frequency<-as.ordered(tpfit$frequency)
for(f in 1:nf){
  levels(tpfit$frequency)[f]<-paste0(levels(tpfit$frequency)[f]," Hz")
}
cl<-max(abs(tpfit$fit))
p<-ggplot(tpfit[tpfit$session=="1",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi),
                     labels=c(expression(paste("-",pi)),"0",expression(paste(pi))),
                     expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dphasetimesession1_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time/theta phase interactions by frequency session 1-3
p<-ggplot(tpfit[tpfit$session=="2",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi),
                     labels=c(expression(paste("-",pi)),"0",expression(paste(pi))),
                     expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dphasetimesession2_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time/theta phase interactions by frequency session 2-3
p<-ggplot(tpfit[tpfit$session=="3",],aes(x=time,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi),
                     labels=c(expression(paste("-",pi)),"0",expression(paste(pi))),
                     expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  stat_contour(mapping=aes(x=time,y=phase,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=phase,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dphasetimesession3_Alz.eps",device="eps",width=5.5,height=11,units="in")

#time theta phase interaction tiled by time
#run for lpmatrix
tpfit<-as.data.frame(ci[tpintsimidx,])
colnames(tpfit)<-c("low","fit","high")
tpfit$low<-tpfit$low+tpfit$fit
tpfit$high<-tpfit$high+tpfit$fit
tpfit$session<-as.factor(tpintfit[,2])
tpfit$frequency<-tpintfit[,3]
tpfit$phase<-tpintfit[,4]
tpfit$time<-tpintfit[,5]
tpfit<-tpfit[tpfit$session=="1",]
tpfit$fit<-tpfit$fit + rep(ci[psimidx,2],each=nb)
keep<-(round((1:nb/nb)/0.2,digits=3) %% 1) == 0
keep[1]<-TRUE
keep<-rep(keep,nb*nf)
tpfit<-tpfit[keep,]

#interpolate for each time 
#run for "lpmatrix"
pbins<-seq(pbins[1],pbins[length(pbins)],length=250)
fbins<-seq(min(freq),max(freq),length=250)
t<-unique(tpfit$time)
fitinterp<-interp(x=tpfit$phase[tpfit$time==t[1]],
                  y=tpfit$frequency[tpfit$time==t[1]],
                  tpfit$fit[tpfit$time==t[1]],
                  xo=pbins,yo=fbins)
highinterp<-interp(x=tpfit$phase[tpfit$time==t[1]],
                   y=tpfit$frequency[tpfit$time==t[1]],
                   tpfit$high[tpfit$time==t[1]],
                   xo=pbins,yo=fbins)
lowinterp<-interp(x=tpfit$phase[tpfit$time==t[1]],
                  y=tpfit$frequency[tpfit$time==t[1]],
                  tpfit$low[tpfit$time==t[1]],
                  xo=pbins,yo=fbins)
tfit <- data.frame(expand.grid(phase = fitinterp$x, frequency = fitinterp$y), 
                   fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z))
for(r in 2:length(t)){
  fitinterp<-interp(x=tpfit$phase[tpfit$time==t[r]],
                    y=tpfit$frequency[tpfit$time==t[r]],
                    tpfit$fit[tpfit$time==t[r]],
                    xo=pbins,yo=fbins)
  highinterp<-interp(x=tpfit$phase[tpfit$time==t[r]],
                     y=tpfit$frequency[tpfit$time==t[r]],
                     tpfit$high[tpfit$time==t[r]],
                     xo=pbins,yo=fbins)
  lowinterp<-interp(x=tpfit$phase[tpfit$time==t[r]],
                    y=tpfit$frequency[tpfit$time==t[r]],
                    tpfit$low[tpfit$time==t[r]],
                    xo=pbins,yo=fbins)
  tfit <- rbind(tfit,data.frame(expand.grid(phase = fitinterp$x, frequency = fitinterp$y), 
                                fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z)))
}
tfit$time<-round(rep(t,each=250^2),digits = 2)

#display plot
cl<-max(abs(tfit$fit))
p<-ggplot(tfit,aes(x=phase,y=frequency,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster()+scale_x_continuous(breaks=seq(-pi,pi,pi/2),labels=c("-pi","-pi/2","0","pi/2","pi"),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,250,by=25),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  # stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  # stat_contour(mapping=aes(x=phase,y=frequency,z=high),colour="purple",breaks=c(0),alpha=1)+
  # stat_contour(mapping=aes(x=phase,y=frequency,z=low),colour="orange",breaks=c(0),alpha=1)+
  facet_wrap(~time)