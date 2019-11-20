#AA_coupling plots

library(ggplot2)
library(MASS)
library(itsadug)
library(RColorBrewer)
library(directlabels)
library(scales)

setwd("C:/Data/GAMS_AA_coupling")

nf<-50
jobid<-"7509645"
freq<-10^seq(log10(2),log10(100),length.out=nf)
mapRsq<-matrix(NA,nrow=nf*nf,ncol=3)
nb<-100
nv<-3
mina<- -4
maxa<- 4
abins<-seq(mina,maxa,length.out=nb)
dataAmp<-matrix(NA,nrow=nf*nb,ncol=nv*nf)

for(f in 1:nf)
{
  filename<-paste0(getwd(),"/AARsq_",jobid,"_",f,".rds")
  if(file.exists(filename))
  {
    Rsq<-readRDS(filename)
    mapRsq[((f-1)*nf+1):(f*nf),1]<-Rsq
    mapRsq[((f-1)*nf+1):(f*nf),2]<-rep(freq[f],nf)
    mapRsq[((f-1)*nf+1):(f*nf),3]<-freq
  }
  
  filename<-paste0(getwd(),"/modpred_",jobid,"_",f,".rds")
  if(file.exists(filename))
  {
    Amp<-readRDS(filename)
    dataAmp[((f-1)*nb+1):(f*nb),]<-Amp[1:nb,]
    
    h<-c(freq[1]-(freq[2]-freq[1])/2,0*freq)
    for(ff in 1:length(freq))
    {
      h[ff+1]<-freq[ff]+freq[ff]-h[ff]
    }
    w<-diff(abins)
    w<-rep(w[1],nb*nf)
    h<-rep(diff(h),each=nb)
    
    #plot amplitude response to each frequency
    mapAmp<-as.data.frame(matrix(Amp[1:nb,seq(1,nv*nf-nv+1,by=nv)],nrow=nb*nf,ncol=1))
    colnames(mapAmp)<-c("response_amplitude")
    mapAmp$frequency<-rep(freq,rep(nb,nf))
    mapAmp$predictor_amplitude<-rep(abins,nf)
    cl<-0.4
    p<-ggplot(mapAmp,aes(x=predictor_amplitude,y=frequency,fill=response_amplitude,z=response_amplitude))
    myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
    print(p+geom_tile(height=h,width=w)+scale_x_continuous(breaks=seq(round(mina),round(maxa),by=1),expand=c(0,0))+
            scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
            scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
            #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.ticks = element_blank())+
            ggtitle(paste0(toString(round(freq[f],1))," Hz")))
  }
}




#plot R-squared map
fbins<-seq(2,100,length=200)
mapRsq[is.na(mapRsq)]<-1
fitinterp<-interp(x=mapRsq[,3],y=mapRsq[,2],mapRsq[,1],xo=fbins,yo=fbins)
fit <- data.frame(expand.grid(pf = fitinterp$x, rf = fitinterp$y), rsq = c(fitinterp$z))
cl<-0.05
mapRsq<-as.data.frame(mapRsq)
p<-ggplot(fit,aes(x=pf,y=rf,fill=rsq,z=rsq))
myPalette <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")), space="Lab")
p+geom_tile()+scale_x_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(0,cl),oob=squish)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        aspect.ratio = 1)


#set pixel hieghts and widths for log scaling
h<-c(freq[1]-(freq[2]-freq[1])/2,0*freq)
for(ff in 1:length(freq))
{
  h[ff+1]<-freq[ff]+freq[ff]-h[ff]
}
w<-diff(abins)
w<-rep(w[1],nb*nf)
h<-rep(diff(h),each=nb)

#plot amplitude response to each frequency
dataAmp<-as.data.frame(matrix(dataAmp[,seq(1,nv*nf-nv+1,by=nv)],ncol=1))
colnames(dataAmp)<-c("response_amplitude")
dataAmp$response_frequency<-rep(rep(round(freq,1),each=nb),nf)
dataAmp$predictor_frequency<-rep(freq,each=nb*nf)
dataAmp$predictor_amplitude<-rep(rep(abins,nf),nf)
cl<-0.4
p<-ggplot(dataAmp,aes(x=predictor_amplitude,y=predictor_frequency,fill=response_amplitude,z=response_amplitude))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
p+geom_tile(height=h,width=w)+scale_x_continuous(breaks=seq(round(mina),round(maxa),by=2),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=25),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  #stat_contour(colour="black",binwidth=0.05,size=0.05,alpha=0.3,show.legend = TRUE)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())+
  facet_wrap(~response_frequency,nrow=8)+theme(aspect.ratio = 1)