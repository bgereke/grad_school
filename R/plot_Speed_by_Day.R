library(mgcv)
library(ggplot2)
library(MASS)
library(itsadug)
library(RColorBrewer)
library(wesanderson)
library(directlabels)
library(scales)
library(akima)

setwd("G:/Rats/GAMS_full_CA3")

#load my data
cwd<-getwd()
MYDATA<-readRDS("MYDATAPOW_full_CA3.rds")

#Set model prediction vars
numdays<-length(unique(MYDATA$Day))
days<-seq(1,numdays)
nb<-45
nf<-50
freq<-10^seq(log10(2),log10(100),length.out=nf)
minrs<-min(MYDATA$Running_Speed)
maxrs<-max(MYDATA$Running_Speed)
rsbins<-seq(minrs,maxrs,length.out=nb)

#create new data
dgrid<-expand.grid(rsbins,unique(MYDATA$Day))
colnames(dgrid)<-c("runningspeed","day")
newdat<-with(MYDATA,data.frame(Running_Speed=dgrid$runningspeed,Day=dgrid$day))

#preallocate space
pred <- matrix(0,nrow = numdays*nb*nf,ncol = 1)

for(f in 1:nf)
{
  filename<-paste0("CA3_f",f,"_POW.rds")
  fmod<-readRDS(filename)

  #get predictions
  pred[((f-1)*numdays*nb+1):(f*numdays*nb)]<-predict.bam(fmod, newdat[,], newdata.guaranteed = TRUE)
  
  print(paste("Completed: ",f," of ",nf),quote=FALSE)
}

#build data frame
dfspeed <- data.frame("Pred" = pred,
                      "Speed" = rep(rsbins,times = nf*numdays),
                      "Freq" = rep(1:nf,each = nb*numdays),
                      "Day" = rep(newdat$Day,times = nf))
dfspeed$Day <- as.factor(dfspeed$Day)

#display plot
cl<-0.5#0.25*max(abs(dfspeed$Pred))
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
ggplot(dfspeed,aes(x=Speed,y=Freq,fill=Pred,z=Pred))+ 
  geom_raster()+scale_x_continuous(breaks=seq(0,1,by=0.2),labels=seq(0,1,by=0.2),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(-cl,cl),oob=squish)+
  facet_wrap(~Day,nrow=6)+theme(aspect.ratio = 1)+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(aspect.ratio = 1,panel.spacing=unit(0.1,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        panel.border = element_rect(color="black",fill=NA,size=1))

