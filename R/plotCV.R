library(mgcv)
library(RevoUtilsMath)
library(ggplot2)
library(MASS)
library(itsadug)
library(RColorBrewer)
library(wesanderson)
library(directlabels)
library(scales)
library(matrixStats)
library(parallel)

setwd("C:/Data/GAMS_full")
MYDATA<-readRDS("MYDATAPOW_full.rds")

#Set some general variables
nummice<-length(unique(MYDATA$Mouse))
nf<-50
nm<-length(unique(MYDATA$Mouse))
freq<-10^seq(log10(2),log10(100),length.out=nf)

count<-0

setwd("C:/Data/GAMS_CV")
no_cores <- detectCores() #- 1
cl<-makeCluster(no_cores, type='PSOCK')
clusterExport(cl,list("nf","MYDATA"))
clusterEvalQ(cl,{library(mgcv);NULL})

results<-parLapply(cl,1:nm,function(m)
{
  #preallocate predictions
  SSEtot<-matrix(0,nrow=1,ncol=nf)
  SSEpred<-matrix(0,nrow=1,ncol=nf)
  
  #get predictor data for this day
  mouse <- unique(MYDATA$Mouse)[m]
  testdat<-MYDATA[MYDATA$Mouse==mouse,]
  testdat$Day<-MYDATA$Day[MYDATA$Mouse!=mouse][1]
  testdat$Trial<-MYDATA$Trial[MYDATA$Mouse!=mouse][1]
  
  predmat <- matrix(nrow = length(testdat$Mouse), ncol = nf)
  
  setwd("C:/Data/GAMS_CV_rho")
  for(f in 1:nf)
  {
    #load model
    filename<-paste0("Base_TS_TP_f",f,"_m",m,"_POW.rds")
    fmod<-readRDS(filename)
    
    #get prediciton
    fit<-predict.bam(fmod,testdat,type="terms",newdata.guaranteed = TRUE)
    
    randeff1<-grep("Day",colnames(fit))
    randeff2<-grep("Trial",colnames(fit))
    randeff<-append(randeff1,randeff2)
    predmat[,f] <- rowSums(as.matrix(fit[,-randeff]))
    predmat[,f] <- rowSums(as.matrix(fit))
    
    SSEtot[f]<-SSEtot[f] + sum((testdat[,f]-mean(testdat[,f]))^2)
    SSEpred[f]<-SSEpred[f] + sum((predmat[,f]-testdat[,f])^2)
    
    #count<-count+1
    #print(paste("Completed: ",count," of ",nf*nummice),quote=FALSE)
  }
  #results<-append(results,list(predmat,SSEtot,SSEpred))
  list(predmat,SSEtot,SSEpred)
})
stopCluster(cl)

saveRDS(results,"Base_TS_TP_lomo.rds")
# registerDoSEQ()

setwd("C:/Data/GAMS_CV")
B<-readRDS("Base_lomo.rds")
BT<-readRDS("Base_T_lomo.rds")
BTS<-readRDS("Base_TS_lomo.rds")
BTSRT<-readRDS("Base_TS_RT_lomo.rds")
BTSTP<-readRDS("Base_TS_TP_lomo.rds")
BTSRTTP<-readRDS("Base_TS_RT_TP_lomo.rds")

SSEtot_B<-matrix(0,nrow=1,ncol=nf)
SSEpred_B<-matrix(0,nrow=1,ncol=nf)
SSEtot_BT<-matrix(0,nrow=1,ncol=nf)
SSEpred_BT<-matrix(0,nrow=1,ncol=nf)
SSEtot_BTS<-matrix(0,nrow=1,ncol=nf)
SSEpred_BTS<-matrix(0,nrow=1,ncol=nf)
SSEtot_BTSRT<-matrix(0,nrow=1,ncol=nf)
SSEpred_BTSRT<-matrix(0,nrow=1,ncol=nf)
SSEtot_BTSTP<-matrix(0,nrow=1,ncol=nf)
SSEpred_BTSTP<-matrix(0,nrow=1,ncol=nf)
SSEtot_BTSRTTP<-matrix(0,nrow=1,ncol=nf)
SSEpred_BTSRTTP<-matrix(0,nrow=1,ncol=nf)

for(m in 1:nm){
  SSEtot_B<-SSEtot_B+B[[m]][[2]]
  SSEpred_B<-SSEpred_B+B[[m]][[3]]
  SSEtot_BT<-SSEtot_BT+BT[[m]][[2]]
  SSEpred_BT<-SSEpred_BT+BT[[m]][[3]]
  SSEtot_BTS<-SSEtot_BTS+BTS[[m]][[2]]
  SSEpred_BTS<-SSEpred_BTS+BTS[[m]][[3]]
  SSEtot_BTSRT<-SSEtot_BTSRT+BTSRT[[m]][[2]]
  SSEpred_BTSRT<-SSEpred_BTSRT+BTSRT[[m]][[3]]
  SSEtot_BTSTP<-SSEtot_BTSTP+BTSTP[[m]][[2]]
  SSEpred_BTSTP<-SSEpred_BTSTP+BTSTP[[m]][[3]]
  SSEtot_BTSRTTP<-SSEtot_BTSRTTP+BTSRTTP[[m]][[2]]
  SSEpred_BTSRTTP<-SSEpred_BTSRTTP+BTSRTTP[[m]][[3]]
}

#plot r-squared for all

Rsq_B <- 1 - SSEpred_B/SSEtot_B
Rsq_BT <- 1 - SSEpred_BT/SSEtot_BT
Rsq_BTS <- 1 - SSEpred_BTS/SSEtot_BTS
Rsq_BTSRT <- 1 - SSEpred_BTSRT/SSEtot_BTSRT
Rsq_BTSTP <- 1 - SSEpred_BTSTP/SSEtot_BTSTP
Rsq_BTSRTTP <- 1 - SSEpred_BTSRTTP/SSEtot_BTSRTTP

rsqdf<-data.frame(Frequency=rep(freq,times=6),Rsquared=c(Rsq_B,Rsq_BT,Rsq_BTS,
                                                         Rsq_BTSRT,Rsq_BTSTP,
                                                         Rsq_BTSRTTP))
rsqdf$model<-factor(rep(c("B","BT","BTS","BTSRT","BTSTP","BTSRTTP"),each=nf),
                    levels = c("B","BT","BTS","BTSRT","BTSTP","BTSRTTP"))
# jBrewColors <- brewer.pal(n = 6, name = "Dark2",type="qual")
jBrewColors <- wes_palette(n = 5, name = "Darjeeling",type="discrete")
cols<-c("black",jBrewColors)
minrsq<-min(rsqdf$Rsquared)
maxrsq<-max(rsqdf$Rsquared)
ggplot(rsqdf,aes(x=Frequency,y=Rsquared,group=model))+
  geom_line(aes(colour=model),size=1)+
  # scale_color_brewer(palette="Dark2",type="qual")+
  scale_color_manual(values=cols)+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,0.2,by=0.025),name="R^2",expand=c(0,0),limits=c(minrsq-0.001,maxrsq+0.001))+
  geom_hline(aes(yintercept=0),linetype="dashed")+
  theme(axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.background = element_blank(),aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA,size=1))

ggsave("allcv.eps",device="eps",width=7,height=7,units="in")

#plot r-squared (around mean)
Rsq_B <- 1 - SSEpred_B/SSEtot_B
Rsq_BT <- 1 - SSEpred_BT/SSEtot_BT
Rsq_BTS <- 1 - SSEpred_BTS/SSEtot_BTS
Rsq_BTSRT <- 1 - SSEpred_BTSRT/SSEtot_BTSRT
Rsq_BTSTP <- 1 - SSEpred_BTSTP/SSEtot_BTSTP
Rsq_BTSRTTP <- 1 - SSEpred_BTSRTTP/SSEtot_BTSRTTP

Rsq_mu<-(Rsq_BT+Rsq_BTS+Rsq_BTSRT+Rsq_BTSTP+Rsq_BTSRTTP)/5

Rsq_B <- Rsq_B - Rsq_mu
Rsq_BT <- Rsq_BT - Rsq_mu
Rsq_BTS <- Rsq_BTS - Rsq_mu
Rsq_BTSRT <- Rsq_BTSRT - Rsq_mu
Rsq_BTSTP <- Rsq_BTSTP - Rsq_mu
Rsq_BTSRTTP <- Rsq_BTSRTTP - Rsq_mu

rsqdf<-data.frame(Frequency=rep(freq,times=5),Rsquared=c(Rsq_BT,Rsq_BTS,Rsq_BTSRT,
                                                         Rsq_BTSTP,Rsq_BTSRTTP))
rsqdf$model<-factor(rep(c("BT","BTS","BTSRT","BTSTP","BTSRTTP"),each=nf),
                    levels = c("BT","BTS","BTSRT","BTSTP","BTSRTTP"))

jBrewColors <- wes_palette(n = 5, name = "Darjeeling",type="discrete")
cols<-c(jBrewColors)

ggplot(rsqdf,aes(x=Frequency,y=Rsquared,group=model))+
  geom_line(aes(colour=model),size=1)+
  # scale_color_brewer(palette="Dark2",type="qual")+
  scale_color_manual(values=cols)+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(name="R^2",expand=c(0,0))+
  geom_hline(aes(yintercept=0),linetype="dashed")+
  theme(axis.text = element_text(size = 20,color="black"),
        axis.title = element_text(size = 24,color="black"),
        panel.background = element_blank(),aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA,size=1))

ggsave("dmcv.eps",device="eps",width=7.25,height=7.25,units="in")

#plot r-squared beyond base
Rsq_BT <- 1 - SSEpred_BT/SSEtot_BT - Rsq_B
Rsq_BTS <- 1 - SSEpred_BTS/SSEtot_BTS - Rsq_B
Rsq_BTSRT <- 1 - SSEpred_BTSRT/SSEtot_BTSRT - Rsq_B
Rsq_BTSTP <- 1 - SSEpred_BTSTP/SSEtot_BTSTP - Rsq_B
Rsq_BTSRTTP <- 1 - SSEpred_BTSRTTP/SSEtot_BTSRTTP - Rsq_B

rsqdf<-data.frame(Frequency=rep(freq,times=5),Rsquared=c(Rsq_BT,Rsq_BTS,Rsq_BTSRT,
                                                         Rsq_BTSTP,Rsq_BTSRTTP))
rsqdf$model<-factor(rep(c("BT","BTS","BTSRT","BTSTP","BTSRTTP"),each=nf),
                    levels = c("BT","BTS","BTSRT","BTSTP","BTSRTTP"))
ggplot(rsqdf,aes(x=Frequency,y=Rsquared,group=model))+
  geom_hline(aes(yintercept=0),linetype="dashed")+
  geom_line(aes(colour=model),size=0.25)+
  scale_color_brewer(palette="Dark2",type="qual")+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(name="R^2",expand=c(0,0))+
  geom_vline(aes(xintercept=10),linetype="dashed",col="gray")+
  geom_vline(aes(xintercept=9),linetype="dashed",col="gray")+
  geom_vline(aes(xintercept=8),linetype="dashed",col="gray")+
  geom_vline(aes(xintercept=7),linetype="dashed",col="gray")+
  geom_vline(aes(xintercept=6),linetype="dashed",col="gray")+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),aspect.ratio = 1)

#plot performance relative to full (difference in r-squared)
Rsq_B <- 1 - SSEpred_B/SSEtot_B - Rsq_BTSRTTP
Rsq_BT <- 1 - SSEpred_BT/SSEtot_BT - Rsq_BTSRTTP
Rsq_BTS <- 1 - SSEpred_BTS/SSEtot_BTS - Rsq_BTSRTTP
Rsq_BTSRT <- 1 - SSEpred_BTSRT/SSEtot_BTSRT - Rsq_BTSRTTP
Rsq_BTSTP <- 1 - SSEpred_BTSTP/SSEtot_BTSTP - Rsq_BTSRTTP

rsqdf<-data.frame(Frequency=rep(freq,times=5),Rsquared=c(Rsq_B,Rsq_BT,Rsq_BTS,Rsq_BTSRT,
                                                         Rsq_BTSTP))
rsqdf$model<-factor(rep(c("B","BT","BTS","BTSRT","BTSTP"),each=nf),
                    levels = c("B","BT","BTSRT","BTSTP","BTS"))
ggplot(rsqdf,aes(x=Frequency,y=Rsquared,group=model))+
  geom_hline(aes(yintercept=0),linetype="dashed")+
  geom_line(aes(colour=model),size=0.75)+
  geom_vline(aes(xintercept=10),linetype="dashed",col="gray")+
  geom_vline(aes(xintercept=9),linetype="dashed",col="gray")+
  geom_vline(aes(xintercept=8),linetype="dashed",col="gray")+
  geom_vline(aes(xintercept=7),linetype="dashed",col="gray")+
  geom_vline(aes(xintercept=6),linetype="dashed",col="gray")+
  scale_color_brewer(palette="Dark2",type="qual")+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(name="R^2",expand=c(0,0))+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),aspect.ratio = 1)


Rsq_BTS <- 1 - SSEpred_BTS/SSEtot_BTS - Rsq_B - Rsq_BT
Rsq_BTSRT <- 1 - SSEpred_BTSRT/SSEtot_BTSRT - Rsq_B - Rsq_BT
Rsq_BTSTP <- 1 - SSEpred_BTSTP/SSEtot_BTSTP - Rsq_B - Rsq_BT
Rsq_BTSRTTP <- 1 - SSEpred_BTSRTTP/SSEtot_BTSRTTP - Rsq_B - Rsq_BT

rsqdf<-data.frame(Frequency=rep(freq,times=4),Rsquared=c(Rsq_BTS,Rsq_BTSRT,
                                                         Rsq_BTSTP,Rsq_BTSRTTP))
rsqdf$model<-factor(rep(c("BTS","BTSRT","BTSTP","BTSRTTP"),each=nf),
                    levels = c("BTS","BTSRT","BTSTP","BTSRTTP"))
ggplot(rsqdf,aes(x=Frequency,y=Rsquared,group=model))+
  geom_hline(aes(yintercept=0),linetype="dashed")+
  geom_line(aes(colour=model),size=0.75)+
  # geom_ribbon(aes(colour=model,fill=model,ymin=low,ymax=high),alpha=0.3)+
  scale_color_brewer(palette="Dark2")+
  # scale_fill_brewer(palette="Dark2")+
  scale_x_continuous(breaks=seq(0,100,by=10),name="frequency (Hz)",expand=c(0,0))+
  scale_y_continuous(name="R^2",expand=c(0,0))+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),aspect.ratio = 1)

# #construct permutation test
# numperms <- 1000
# pSSEf <- matrix(0,nrow=numperms,ncol=nf)
# for(d in 1:numdays)
# {
#   day <- unique(MYDATA$Day)[d]
#   data <- MYDATA[MYDATA$Day==day,1:nf]
#   fdata <- data-predictions[[d]]
#   for(p in 1:numperms)
#   {
#     end <- length(data[,1])
#     ridx <- sample(seq(2,end-1),1)
#     pSSEf[p,] <- pSSEf[p,] + colSums((data[c(ridx:end,1:(ridx-1)),]-predictions[[d]])^2)
#   }
#   print(paste("Completed: ",d," of ",numdays),quote=FALSE)
# }
# 
# pRsq <- 1 - pSSEf/matrix(rep(SSEtot,numperms),numperms,byrow=TRUE)
# cutRsq <- matrix(0,nrow=1,ncol=nf)
# for(f in 1:nf){cutRsq[f]<-quantile(pRsq[,f],probs=c(0.95))}


