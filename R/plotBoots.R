### plot boots ###

library(mgcv)
library(ggplot2)
library(MASS)
library(RColorBrewer)
library(wesanderson)
library(directlabels)
library(scales)
library(akima)

setwd("C:/Data/GAMS_full")

#load my data
cwd<-getwd()
MYDATA<-readRDS("MYDATAPOW_full.rds")
# bootlist<-readRDS("boots_full_POW_90172.rds ")
fits<-readRDS("boots_full_POW.rds")
fits<-as.vector(t(fits))

#concatenate bootstrap fits into matrix
nboots<-100
bootmat<-NULL
tmp<-as.vector(t(readRDS("boot_1.rds")))
bootmat<-matrix(0,nrow=length(tmp),ncol=nboots)
for(b in 1:nboots){
  bootmat[,b]<-as.vector(t(readRDS(paste0("boot_",b,".rds"))))
}
rm(tmp)
gc()

#set prediction grids
numdays<-length(unique(MYDATA$Day))
days<-seq(1,numdays)
nb<-60
nf<-50
freq<-10^seq(log10(2),log10(100),length.out=nf)
mint<-0
maxt<-600
tbins<-seq(mint,maxt,length.out=nb)
minrs<-0
maxrs<-1
rbins<-seq(minrs,maxrs,length.out=nb)
minp<- -pi
maxp<- pi
pbins<-seq(minp,maxp,length.out=nb)

#set plot indices
iidx <- 1:(3*nf)
ridx <- (max(iidx)+1):(max(iidx)+nf*nb)
pidx <- (max(ridx)+1):(max(ridx)+nf*nb)
rpidx <- (max(pidx)+1):(max(pidx)+nf*nb^2)
tsidx <- (max(rpidx)+1):(max(rpidx)+3*nf*nb)
rtsidx <- (max(tsidx)+1):(max(tsidx)+3*nf*nb^2)
ptsidx <- (max(rtsidx)+1):(max(rtsidx)+3*nf*nb^2)

########ONLY RUN WHEN PLOTTING DIFFERENCE TERMS!####################
#convert fits to session differences
#intercepts
dfits<-fits
dfits[1:nf]<-fits[1:nf]-fits[(nf+1):(2*nf)]
dfits[(nf+1):(2*nf)]<-fits[1:nf]-fits[(2*nf+1):(3*nf)]
dfits[(2*nf+1):(3*nf)]<-fits[(nf+1):(2*nf)]-fits[(2*nf+1):(3*nf)]
#time by session
dfits[(max(rpidx)+1):(max(rpidx)+nf*nb)]<-fits[(max(rpidx)+1):(max(rpidx)+nf*nb)]-
  fits[(max(rpidx)+nf*nb+1):(max(rpidx)+2*nf*nb)]
dfits[(max(rpidx)+nf*nb+1):(max(rpidx)+2*nf*nb)]<-fits[(max(rpidx)+1):(max(rpidx)+nf*nb)]-
  fits[(max(rpidx)+2*nf*nb+1):(max(rpidx)+3*nf*nb)]
dfits[(max(rpidx)+2*nf*nb+1):(max(rpidx)+3*nf*nb)]<-fits[(max(rpidx)+nf*nb+1):(max(rpidx)+2*nf*nb)]-
  fits[(max(rpidx)+2*nf*nb+1):(max(rpidx)+3*nf*nb)]
#runnings speed by time by session
dfits[(max(tsidx)+1):(max(tsidx)+nf*nb^2)]<-fits[(max(tsidx)+1):(max(tsidx)+nf*nb^2)]-
  fits[(max(tsidx)+nf*nb^2+1):(max(tsidx)+2*nf*nb^2)]
dfits[(max(tsidx)+nf*nb^2+1):(max(tsidx)+2*nf*nb^2)]<-fits[(max(tsidx)+1):(max(tsidx)+nf*nb^2)]-
  fits[(max(tsidx)+2*nf*nb^2+1):(max(tsidx)+3*nf*nb^2)]
dfits[(max(tsidx)+2*nf*nb^2+1):(max(tsidx)+3*nf*nb^2)]<-fits[(max(tsidx)+nf*nb^2+1):(max(tsidx)+2*nf*nb^2)]-
  fits[(max(tsidx)+2*nf*nb^2+1):(max(tsidx)+3*nf*nb^2)]
#theta phase by time by session
dfits[(max(rtsidx)+1):(max(rtsidx)+nf*nb^2)]<-fits[(max(rtsidx)+1):(max(rtsidx)+nf*nb^2)]-
  fits[(max(rtsidx)+nf*nb^2+1):(max(rtsidx)+2*nf*nb^2)]
dfits[(max(rtsidx)+nf*nb^2+1):(max(rtsidx)+2*nf*nb^2)]<-fits[(max(rtsidx)+1):(max(rtsidx)+nf*nb^2)]-
  fits[(max(rtsidx)+2*nf*nb^2+1):(max(rtsidx)+3*nf*nb^2)]
dfits[(max(rtsidx)+2*nf*nb^2+1):(max(rtsidx)+3*nf*nb^2)]<-fits[(max(rtsidx)+nf*nb^2+1):(max(rtsidx)+2*nf*nb^2)]-
  fits[(max(rtsidx)+2*nf*nb^2+1):(max(rtsidx)+3*nf*nb^2)]
fits<-dfits
rm(dfits)
gc()

#convert bootfits to session differences
#intercepts
dbootmat<-bootmat
dbootmat[1:nf,]<-bootmat[1:nf,]-bootmat[(nf+1):(2*nf),]
dbootmat[(nf+1):(2*nf),]<-bootmat[1:nf,]-bootmat[(2*nf+1):(3*nf),]
dbootmat[(2*nf+1):(3*nf),]<-bootmat[(nf+1):(2*nf),]-bootmat[(2*nf+1):(3*nf),]
#time by session
dbootmat[(max(rpidx)+1):(max(rpidx)+nf*nb),]<-bootmat[(max(rpidx)+1):(max(rpidx)+nf*nb),]-
  bootmat[(max(rpidx)+nf*nb+1):(max(rpidx)+2*nf*nb),]
dbootmat[(max(rpidx)+nf*nb+1):(max(rpidx)+2*nf*nb),]<-bootmat[(max(rpidx)+1):(max(rpidx)+nf*nb),]-
  bootmat[(max(rpidx)+2*nf*nb+1):(max(rpidx)+3*nf*nb),]
dbootmat[(max(rpidx)+2*nf*nb+1):(max(rpidx)+3*nf*nb),]<-bootmat[(max(rpidx)+nf*nb+1):(max(rpidx)+2*nf*nb),]-
  bootmat[(max(rpidx)+2*nf*nb+1):(max(rpidx)+3*nf*nb),]
#runnings speed by time by session
dbootmat[(max(tsidx)+1):(max(tsidx)+nf*nb^2),]<-bootmat[(max(tsidx)+1):(max(tsidx)+nf*nb^2),]-
  bootmat[(max(tsidx)+nf*nb^2+1):(max(tsidx)+2*nf*nb^2),]
dbootmat[(max(tsidx)+nf*nb^2+1):(max(tsidx)+2*nf*nb^2),]<-bootmat[(max(tsidx)+1):(max(tsidx)+nf*nb^2),]-
  bootmat[(max(tsidx)+2*nf*nb^2+1):(max(tsidx)+3*nf*nb^2),]
dbootmat[(max(tsidx)+2*nf*nb^2+1):(max(tsidx)+3*nf*nb^2),]<-bootmat[(max(tsidx)+nf*nb^2+1):(max(tsidx)+2*nf*nb^2),]-
  bootmat[(max(tsidx)+2*nf*nb^2+1):(max(tsidx)+3*nf*nb^2),]
#theta phase by time by session
dbootmat[(max(rtsidx)+1):(max(rtsidx)+nf*nb^2),]<-bootmat[(max(rtsidx)+1):(max(rtsidx)+nf*nb^2),]-
  bootmat[(max(rtsidx)+nf*nb^2+1):(max(rtsidx)+2*nf*nb^2),]
dbootmat[(max(rtsidx)+nf*nb^2+1):(max(rtsidx)+2*nf*nb^2),]<-bootmat[(max(rtsidx)+1):(max(rtsidx)+nf*nb^2),]-
  bootmat[(max(rtsidx)+2*nf*nb^2+1):(max(rtsidx)+3*nf*nb^2),]
dbootmat[(max(rtsidx)+2*nf*nb^2+1):(max(rtsidx)+3*nf*nb^2),]<-bootmat[(max(rtsidx)+nf*nb^2+1):(max(rtsidx)+2*nf*nb^2),]-
  bootmat[(max(rtsidx)+2*nf*nb^2+1):(max(rtsidx)+3*nf*nb^2),]
bootmat<-dbootmat
rm(dbootmat)
gc()
##############END OF STUFF FOR DIFFERENCE TERMS#########################################

#compute simultaneous confidence intervals
ci<-matrix(0,nrow = nrow(bootmat),ncol = 3) #matrix to hold results
boot_mu<-rowMeans(bootmat)
ci[,2]<-fits
# boot_sd<-rep(0,times=nrow(bootmat))
# for(i in 1:nrow(bootmat)){
#   boot_sd[i]<-sd(bootmat[i,])
#   bootmat[i,]<-abs(bootmat[i,]-boot_mu[i])/boot_sd[i]
# }
rm(fits)
gc()
# rankmat<-bootmat
# for(i in 1:nrow(bootmat)){
#   rankmat[i,]<-rank(bootmat[i,]-ci[i,2])/nboots
# }
# maxima<-rep(0,times=nboots)
# minima<-rep(0,times=nboots)
# for(i in 1:nboots){
#   maxima[i]<-max(rankmat[,i])
#   minima[i]<-min(rankmat[,i])
#   }
# mx <- quantile(maxima,probs=0.975,type=8)
# mn <- quantile(minima,probs=0.025,type=8)
for(i in 1:nrow(bootmat)){
  ci[i,1]<-quantile(bootmat[i,]-boot_mu[i],probs=0.025,type=8)
  ci[i,3]<-quantile(bootmat[i,]-boot_mu[i],probs=0.975,type=8)
}
# ci[,1]<--m*boot_sd
# ci[,3]<-m*boot_sd
rm(bootmat,rankmat,minima,maxima)
gc()

#create data frames
idf <- data.frame(frequency = rep(freq,times=3),
                  session = rep(c("1","2","3"),each=nf),
                  fit = ci[iidx,2],
                  low = ci[iidx,1] + ci[iidx,2],
                  high = ci[iidx,3] + ci[iidx,2])
rdf <- data.frame(frequency = rep(freq,times=nb),
                  speed = rep(rbins,each=nf),
                  fit = ci[ridx,2],
                  low = ci[ridx,1] + ci[ridx,2],
                  high = ci[ridx,3] + ci[ridx,2])
pdf <- data.frame(frequency = rep(freq,times=nb),
                  phase = rep(pbins,each=nf),
                  fit = ci[pidx,2],
                  low = ci[pidx,1] + ci[pidx,2],
                  high = ci[pidx,3] + ci[pidx,2])
rpdf <- data.frame(frequency = rep(freq,times=nb^2),
                   speed = rep(rep(rbins,each=nf),times=nb),
                   phase = rep(pbins,each=nf*nb),
                   fit = ci[rpidx,2],
                   low = ci[rpidx,1] + ci[rpidx,2],
                   high = ci[rpidx,3] + ci[rpidx,2])
tsdf <- data.frame(frequency = rep(freq,times=3*nb),
                   time = rep(rep(tbins,each=nf),times=3),
                   session = rep(c("1","2","3"),each=nf*nb),
                   fit = ci[tsidx,2],
                   low = ci[tsidx,1] + ci[tsidx,2],
                   high = ci[tsidx,3] + ci[tsidx,2])
rtsdf <- data.frame(frequency = rep(freq,times=3*nb^2),
                    speed = rep(rep(rbins,each=nf),times=3*nb),
                    time = rep(rep(tbins,each=nf*nb),times=3),
                    session = rep(c("1","2","3"),each=nf*nb^2),
                    fit = ci[rtsidx,2],
                    low = ci[rtsidx,1] + ci[rtsidx,2],
                    high = ci[rtsidx,3] + ci[rtsidx,2])
ptsdf <- data.frame(frequency = rep(freq,times=3*nb^2),
                    phase = rep(rep(pbins,each=nf),times=3*nb),
                    time = rep(rep(tbins,each=nf*nb),times=3),
                    session = rep(c("1","2","3"),each=nf*nb^2),
                    fit = ci[ptsidx,2],
                    low = ci[ptsidx,1] + ci[ptsidx,2],
                    high = ci[ptsidx,3] + ci[ptsidx,2])

#make plots
#intercepts by session
p<-ggplot(idf,aes(x=frequency,y=fit))
p+geom_ribbon(aes(ymin=low,ymax=high))+
  geom_line(size=1)+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  geom_hline(aes(yintercept=0),linetype="dashed")+
  facet_wrap(~session,nrow=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.ticks = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA))

#running speed plot
#interpolate
irbins<-seq(min(rbins),max(rbins),length=250)
ifbins<-seq(min(freq),max(freq),length=250)
fitinterp<-interp(x=rdf$speed,y=rdf$frequency,rdf$fit,xo=irbins,yo=ifbins)
highinterp<-interp(x=rdf$speed,y=rdf$frequency,rdf$high,xo=irbins,yo=ifbins)
lowinterp<-interp(x=rdf$speed,y=rdf$frequency,rdf$low,xo=irbins,yo=ifbins)
rfit <- data.frame(expand.grid(speed = fitinterp$x, frequency = fitinterp$y), 
                   fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z))

#display plot
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
cl<-max(abs(rfit$fit))
ggplot(rfit,aes(x=speed,y=frequency,fill=fit,z=fit))+
  geom_raster()+scale_x_continuous(breaks=seq(0,1,by=0.2),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  stat_contour(mapping=aes(x=speed,y=frequency,z=high),colour="blue",breaks=c(0))+
  stat_contour(mapping=aes(x=speed,y=frequency,z=low),colour="red",breaks=c(0))+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("speed.eps",device="eps",width=5,height=5,units="in")

#theta phase plot
#interpolate
ipbins<-seq(pbins[1],pbins[length(pbins)],length=250)
ifbins<-seq(min(freq),max(freq),length=250)
fitinterp<-interp(x=pdf$phase,y=pdf$frequency,pdf$fit,xo=ipbins,yo=ifbins)
highinterp<-interp(x=pdf$phase,y=pdf$frequency,pdf$high,xo=ipbins,yo=ifbins)
lowinterp<-interp(x=pdf$phase,y=pdf$frequency,pdf$low,xo=ipbins,yo=ifbins)
pfit <- data.frame(expand.grid(phase = fitinterp$x, frequency = fitinterp$y), 
                   fit = c(fitinterp$z), high = c(highinterp$z), low = c(lowinterp$z))

#display plot
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
cl<-max(abs(pfit$fit))
ggplot(pfit,aes(x=phase,y=frequency,fill=fit,z=fit))+
  geom_raster()+scale_x_continuous(breaks=seq(-pi,pi,by=pi/2),
                                   labels=c(expression(paste("-",pi)),expression(paste("-",pi,"/2")),
                                            "0",expression(paste(pi,"/2")),expression(paste(pi))),
                                   expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  stat_contour(mapping=aes(x=phase,y=frequency,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=phase,y=frequency,z=low),colour="red",breaks=c(0),alpha=1)+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("phase.eps",device="eps",width=5,height=5,units="in")

#running speed theta phase interaction tiled by frequency
#display plot
rpdf$frequency<-round(rpdf$frequency,digits=1)
rpdf$frequency<-as.ordered(rpdf$frequency)
for(f in 1:nf){
  levels(rpdf$frequency)[f]<-paste0(levels(rpdf$frequency)[f]," Hz")
}
cl<-max(abs(rpdf$fit))
p<-ggplot(rpdf,aes(x=speed,y=phase,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,1,by=0.5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-pi,pi,by=pi),labels=c(expression(paste("-",pi)),"0",expression(paste(pi))),
                     expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  stat_contour(mapping=aes(x=speed,y=phase,z=high),colour="blue",breaks=c(0),size=0.5)+
  stat_contour(mapping=aes(x=speed,y=phase,z=low),colour="red",breaks=c(0),size=0.5)+
  facet_wrap(~frequency,nrow=5) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("speedphase.eps",device="eps",width=14,height=14,units="in")

#time by session
#interpolate
itbins<-seq(min(tbins),max(tbins),length=250)
ifbins<-seq(min(freq),max(freq),length=250)
idxone<-tsdf$session == "1"
idxtwo<-tsdf$session == "2"
idxthree<-tsdf$session == "3"
fitsone<-interp(x=tsdf$time[idxone],y=tsdf$frequency[idxone],tsdf$fit[idxone],xo=itbins,yo=ifbins)
fitstwo<-interp(x=tsdf$time[idxtwo],y=tsdf$frequency[idxtwo],tsdf$fit[idxtwo],xo=itbins,yo=ifbins)
fitsthree<-interp(x=tsdf$time[idxthree],y=tsdf$frequency[idxthree],tsdf$fit[idxthree],xo=itbins,yo=ifbins)
highone<-interp(x=tsdf$time[idxone],y=tsdf$frequency[idxone],tsdf$high[idxone],xo=itbins,yo=ifbins)
hightwo<-interp(x=tsdf$time[idxtwo],y=tsdf$frequency[idxtwo],tsdf$high[idxtwo],xo=itbins,yo=ifbins)
highthree<-interp(x=tsdf$time[idxthree],y=tsdf$frequency[idxthree],tsdf$high[idxthree],xo=itbins,yo=ifbins)
lowone<-interp(x=tsdf$time[idxone],y=tsdf$frequency[idxone],tsdf$low[idxone],xo=itbins,yo=ifbins)
lowtwo<-interp(x=tsdf$time[idxtwo],y=tsdf$frequency[idxtwo],tsdf$low[idxtwo],xo=itbins,yo=ifbins)
lowthree<-interp(x=tsdf$time[idxthree],y=tsdf$frequency[idxthree],tsdf$low[idxthree],xo=itbins,yo=ifbins)
tsdfone <- data.frame(expand.grid(time = fitsone$x, frequency = fitsone$y), 
                      fit = c(fitsone$z), high = c(highone$z), low = c(lowone$z))
tsdftwo <- data.frame(expand.grid(time = fitstwo$x, frequency = fitstwo$y), 
                      fit = c(fitstwo$z), high = c(hightwo$z), low = c(lowtwo$z))
tsdfthree <- data.frame(expand.grid(time = fitsthree$x, frequency = fitsthree$y), 
                        fit = c(fitsthree$z), high = c(highthree$z), low = c(lowthree$z))
tsdfone$session<-rep(factor("1"),times=length(tsdfone[,1]))
tsdftwo$session<-rep(factor("2"),times=length(tsdftwo[,1]))
tsdfthree$session<-rep(factor("3"),times=length(tsdfthree[,1]))
tfit<-rbind(tsdfone,tsdftwo,tsdfthree)

#display plot
cl<-max(abs(tfit$fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
ggplot(tfit,aes(x=time,y=frequency,fill=fit,z=fit))+ 
  geom_raster()+scale_x_continuous(breaks=seq(0,600,by=60),labels=seq(0,10,by=1),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  stat_contour(mapping=aes(x=time,y=frequency,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=frequency,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~session,nrow=1)+theme(aspect.ratio = 1)+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dtimesession.eps",device="eps",width=14,height=14,units="in")

#time/running speed interactions by frequency session 1
#set NA's
th = 1;
b = 1;
MYDATA$Time<-MYDATA$Time-min(MYDATA$Time)
MYDATA$Time<-MYDATA$Time/max(MYDATA$Time)
rtsdf$time<-rtsdf$time-min(rtsdf$time)
rtsdf$time<-rtsdf$time/max(rtsdf$time)
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
counts1<-rep(as.vector(t(counts1)),each=nf)
counts2<-rep(as.vector(t(counts2)),each=nf)
counts3<-rep(as.vector(t(counts3)),each=nf)
counts<-c(counts1,counts2,counts3)
rtsdf$fit[counts<th]<-NA
rtsdf$low[counts<th]<-NA
rtsdf$high[counts<th]<-NA
#####ONLY RUN WHEN PLOTTING SESSION DIFFERENCES##########
nas1<-rep(0,times=length(counts1))
nas2<-nas1
nas3<-nas1
nas1[counts1<th | counts2<th]<-NA
nas2[counts1<th | counts3<th]<-NA
nas3[counts2<th | counts3<th]<-NA
nas<-c(nas1,nas2,nas3)
rtsdf$fit[is.na(nas)]<-NA
rtsdf$low[is.na(nas)]<-NA
rtsdf$high[is.na(nas)]<-NA
########################################################

#display plot
rtsdf$frequency<-round(rtsdf$frequency,digits=1)
rtsdf$frequency<-as.ordered(rtsdf$frequency)
for(f in 1:nf){
  levels(rtsdf$frequency)[f]<-paste0(levels(rtsdf$frequency)[f]," Hz")
}
cl<-0.6*max(abs(rtsdf$fit),na.rm=TRUE)
p<-ggplot(rtsdf[rtsdf$session=="1",],aes(x=time,y=speed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,1,by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish,na.value="white")+
  stat_contour(mapping=aes(x=time,y=speed,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=speed,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dspeedtimesession1.eps",device="eps",width=5.5,height=11,units="in")

#time/running speed interactions by frequency session 2
p<-ggplot(rtsdf[rtsdf$session=="2",],aes(x=time,y=speed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,1,by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish,na.value="white")+
  stat_contour(mapping=aes(x=time,y=speed,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=speed,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dspeedtimesession2.eps",device="eps",width=5.5,height=11,units="in")

#time/running speed interactions by frequency session 3
p<-ggplot(rtsdf[rtsdf$session=="3",],aes(x=time,y=speed,fill=fit,z=fit))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_raster(interpolate=TRUE)+scale_x_continuous(breaks=seq(0,1,by=0.5),labels=seq(0,10,by=5),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,1,by=0.5),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish,na.value="white")+
  stat_contour(mapping=aes(x=time,y=speed,z=high),colour="blue",breaks=c(0),alpha=1)+
  stat_contour(mapping=aes(x=time,y=speed,z=low),colour="red",breaks=c(0),alpha=1)+
  facet_wrap(~frequency,nrow=10) + 
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 12,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

ggsave("dspeedtimesession3.eps",device="eps",width=5.5,height=11,units="in")

#time theta phase interaction session1
#display plot
ptsdf$frequency<-round(ptsdf$frequency,digits=1)
ptsdf$frequency<-as.ordered(ptsdf$frequency)
for(f in 1:nf){
  levels(ptsdf$frequency)[f]<-paste0(levels(ptsdf$frequency)[f]," Hz")
}
cl<-max(abs(ptsdf$fit))
p<-ggplot(ptsdf[ptsdf$session=="1",],aes(x=time,y=phase,fill=fit,z=fit))
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

ggsave("dphasetimesession1.eps",device="eps",width=5.5,height=11,units="in")

#time theta phase interaction session2
p<-ggplot(ptsdf[ptsdf$session=="2",],aes(x=time,y=phase,fill=fit,z=fit))
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

ggsave("dphasetimesession2.eps",device="eps",width=5.5,height=11,units="in")

#time theta phase interaction session3
p<-ggplot(ptsdf[ptsdf$session=="3",],aes(x=time,y=phase,fill=fit,z=fit))
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

ggsave("dphasetimesession3.eps",device="eps",width=5.5,height=11,units="in")