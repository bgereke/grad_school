library(mgcv)
library(ggplot2)
library(MASS)
library(itsadug)
library(RColorBrewer)
library(wesanderson)
library(directlabels)
library(scales)
library(akima)

#additive tweedie model
p1 <- plot(fit1,rug=FALSE,trans=exp,pages=1,scale=0,scheme=3,n=100)
#2D tweedie model
p2 <- plot(fit2,rug=FALSE,trans=exp,pages=1,scale=0,scheme=3,n2=100)
#additive gamma model
p3 <- plot(fit3,rug=FALSE,trans=exp,pages=1,scale=0,scheme=3,n=100)
#2D tweedie model
p4 <- plot(fit4,rug=FALSE,trans=exp,pages=1,scale=0,scheme=3,n2=100)

#plot additive tweedie model
df<-data.frame(Running_Speed=p1[[1]]$x,
               RSFit=exp(p1[[1]]$fit),
               RSLC=exp(p1[[1]]$fit-p1[[1]]$se), 
               RSUC=exp(p1[[1]]$fit+p1[[1]]$se),
               H=p1[[2]]$x,
               Hfit=exp(p1[[2]]$fit),
               HLC=exp(p1[[2]]$fit-p1[[2]]$se), 
               HUC=exp(p1[[2]]$fit+p1[[2]]$se))
df$Running_Speed <- df$Running_Speed/(2*pi)*100*pi
df$H <- 5*df$H
#Running Speed term
p<-ggplot(df,aes(x=Running_Speed,y=RSFit))
p+geom_ribbon(aes(ymin=RSLC,ymax=RSUC),color="grey1")+
  geom_line(size=1)+
  scale_x_continuous(breaks=seq(0,45,by=10),name="speed (cm/sec)",expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  geom_hline(aes(yintercept=1),linetype="dashed")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.ticks = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA))
#History term
p<-ggplot(df,aes(x=H,y=Hfit))
p+geom_ribbon(aes(ymin=HLC,ymax=HUC),color="grey1")+
  geom_line(size=1)+
  scale_x_continuous(breaks=seq(0,400,by=50),name="Lags (ms)",expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  geom_hline(aes(yintercept=1),linetype="dashed")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.ticks = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA))

#plot additive gamma model
df<-data.frame(Running_Speed=p3[[1]]$x,
               RSFit=exp(p3[[1]]$fit),
               RSLC=exp(p3[[1]]$fit-p3[[1]]$se), 
               RSUC=exp(p3[[1]]$fit+p3[[1]]$se),
               Phase=p3[[2]]$x,
               Phasefit=exp(p3[[2]]$fit),
               PhaseLC=exp(p3[[2]]$fit-p3[[2]]$se), 
               PhaseUC=exp(p3[[2]]$fit+p3[[2]]$se))
df$Running_Speed <- df$Running_Speed/(2*pi)*100*pi
#Running Speed term
p<-ggplot(df,aes(x=Running_Speed,y=RSFit))
p+geom_ribbon(aes(ymin=RSLC,ymax=RSUC),color="grey1")+
  geom_line(size=1)+
  scale_x_continuous(breaks=seq(0,45,by=10),name="speed (cm/sec)",expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  geom_hline(aes(yintercept=1),linetype="dashed")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.ticks = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA))
#Theta Phase term
p<-ggplot(df,aes(x=Phase,y=Phasefit))
p+geom_ribbon(aes(ymin=PhaseLC,ymax=PhaseUC),color="grey1")+
  geom_line(size=1)+
  scale_x_continuous(breaks=seq(-pi,pi,by=pi/2),
                     labels=c(expression(paste("-",pi)),expression(paste("-",pi,"/2")),
                              "0",expression(paste(pi,"/2")),expression(paste(pi))),
                     expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  geom_hline(aes(yintercept=1),linetype="dashed")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA))

#plot 2D tweedie model
df<-data.frame(Running_Speed=rep(p2[[1]]$x,times=100),
               Lags=rep(p2[[1]]$y,each=100),
               Fit=exp(p2[[1]]$fit),
               LC=exp(p2[[1]]$fit-p2[[1]]$se), 
               UC=exp(p2[[1]]$fit+p2[[1]]$se))
df$Running_Speed <- df$Running_Speed/(2*pi)*100*pi
df$Lags <- 5*df$Lags
#2D color plot
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
ggplot(df,aes(x=Lags,y=Running_Speed,fill=Fit,z=Fit))+
  geom_raster()+scale_x_continuous(breaks=seq(0,400,by=50),name="Lags (ms)",expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,45,by=10),name="speed (cm/sec)",expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(1-0.5,1+0.5),oob=squish)+
  stat_contour(mapping=aes(x=Lags,y=Running_Speed,z=UC),colour="blue",breaks=c(1))+
  stat_contour(mapping=aes(x=Lags,y=Running_Speed,z=LC),colour="red",breaks=c(1))+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))

#plot 2D gamma model
df<-data.frame(Running_Speed=rep(p4[[1]]$y,each=100),
               Phase=rep(p4[[1]]$x,times=100),
               Fit=exp(p4[[1]]$fit),
               LC=exp(p4[[1]]$fit-p4[[1]]$se), 
               UC=exp(p4[[1]]$fit+p4[[1]]$se))
df$Running_Speed <- df$Running_Speed/(2*pi)*100*pi
#2D color plot
myPalette <- colorRampPalette(rev(brewer.pal(10, "PuOr")), space="Lab")
ggplot(df,aes(x=Phase,y=Running_Speed,fill=Fit,z=Fit))+
  geom_raster()+ 
  scale_x_continuous(breaks=seq(-pi,pi,by=pi/2),
                     labels=c(expression(paste("-",pi)),expression(paste("-",pi,"/2")),
                     "0",expression(paste(pi,"/2")),expression(paste(pi))),
                     expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,45,by=10),name="speed (cm/sec)",expand=c(0,0))+
  scale_fill_gradientn(colours = myPalette(250),limits = c(1-0.75,1+0.75),oob=squish)+
  stat_contour(mapping=aes(x=Phase,y=Running_Speed,z=UC),colour="blue",breaks=c(1))+
  stat_contour(mapping=aes(x=Phase,y=Running_Speed,z=LC),colour="red",breaks=c(1))+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,panel.spacing=unit(0.2,"lines"),
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        strip.text = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1))
