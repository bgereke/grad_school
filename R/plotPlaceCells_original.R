library(mgcv)
library(MASS)
library(ggplot2)

setwd("C:/Data/PCtable/")

#load data
data <- read.table("PCtable.csv",sep=",",header=TRUE)

#set appropriate variables to factors
data$cell <- ordered(data$cell)
data$day <- ordered(data$day) 
data$mouse <- ordered(data$mouse)
data$session <- ordered(data$session) 
data$com <- data$com/(2*pi)*100*pi #set units to cm
data$peakpos <- data$peakpos/(2*pi)*100*pi #set units to cm
data$firstspk <- data$firstspk/(2*pi)*100*pi #set units to cm
data$lastspk <- data$lastspk/(2*pi)*100*pi #set units to cm
data$lapwidth <- data$lapwidth/(2*pi)*100*pi #set units to cm
data$rewdist <- data$rewdist/(2*pi)*100*pi #set units to cm
data$velenter <- data$velenter/(2*pi)*100*pi #set units to cm
data$velexit <- data$velexit/(2*pi)*100*pi #set units to cm
data$speed <- (data$velenter*data$tenter + data$velexit*data$texit)/(data$tenter+data$texit)
data$rates <- data$numspks/(data$tenter+data$texit)
data$BAI<-(data$velenter-data$velexit)/(data$velenter+data$velexit)
data$mmdiff<-data$com-data$median

#z-score rates and get mean rate for each cell
cells <- unique(data$cell)
nc <- length(cells)
data$zrates <- data$rates
data$murates <- data$rates
for(c in 1:nc){
  mu <- mean(data$rates[data$cell==cells[c]])
  std <- sqrt(sum((data$rates[data$cell==cells[c]]-mu)^2))
  data$murates[data$cell==cells[c]] <- mu
  data$zrates[data$cell==cells[c]] <- (data$rates[data$cell==cells[c]] - mu)/std
  data$firstspk[data$cell==cells[c]]<-data$firstspk[data$cell==cells[c]]-mean(data$firstspk[data$cell==cells[c]])
  data$lastspk[data$cell==cells[c]]<-data$lastspk[data$cell==cells[c]]-mean(data$lastspk[data$cell==cells[c]])
  data$lapwidth[data$cell==cells[c]]<-data$lapwidth[data$cell==cells[c]]/mean(data$lapwidth[data$cell==cells[c]])
  data$fsize[data$cell==cells[c]]<-data$fsize[data$cell==cells[c]]/mean(data$fsize[data$cell==cells[c]])
  data$com[data$cell==cells[c]]<-scale(data$com[data$cell==cells[c]])
  data$median[data$cell==cells[c]]<-scale(data$median[data$cell==cells[c]])
  data$peakpos[data$cell==cells[c]]<-scale(data$peakpos[data$cell==cells[c]])
  data$skewness[data$cell==cells[c]]<-scale(data$skewness[data$cell==cells[c]])
  data$mmdiff[data$cell==cells[c]]<-scale(data$mmdiff[data$cell==cells[c]])
  
}
data <- data[data$murates>2.5,]
data <- data[data$velenter>5,]
data <- data[data$velexit>5,]
data <- data[data$rewdist>5,]

###### COM ######

#run gam 
com <- gam(com ~ session + s(time,by=session,id=1), data = data, method = "REML")

#plot(speed,rug=FALSE,scale=0)

#make prediciton data
nb <- 250
mint <- 0
maxt <- 600
tbins <- seq(mint,maxt,length.out=nb)
sessions <- unique(data$session)
ns <- length(sessions)
pdat<-with(data,data.frame(time=rep(tbins,times=ns),session=rep(sessions,each=nb)))

#simulate simultaneous intervals
nsim <- 10000
vc_com <- vcov(com)

ptype <- "iterms"
pred_com <- predict(com, pdat, type = ptype, se.fit = TRUE)

se_com <- pred_com$se.fit
fit_com <- pred_com$fit
idx1 <- which(pdat$session == "1")
idx2 <- which(pdat$session == "2")
idx3 <- which(pdat$session == "3")
want1<-grep("s\\(time\\)\\:session1",colnames(se_com))
want2<-grep("s\\(time\\)\\:session2",colnames(se_com))
want3<-grep("s\\(time\\)\\:session3",colnames(se_com))
se_com <- c(se_com[idx1,want1], se_com[idx2,want2], se_com[idx3,want3])
fit_comd <- c(fit_com[idx1,want1]-fit_com[idx2,want2], 
              fit_com[idx1,want1]-fit_com[idx3,want3], 
              fit_com[idx2,want2]-fit_com[idx3,want3])
fit_com <- c(fit_com[idx1,want1], fit_com[idx2,want2], fit_com[idx3,want3])

set.seed(42)
B_com <- mvrnorm(nsim, mu = rep(0, nrow(vc_com)), Sigma = vc_com)

LP_com <- predict(com, pdat, type = "lpmatrix")

want1<-grep("s\\(time\\)\\:session1",colnames(LP_com))
want2<-grep("s\\(time\\)\\:session2",colnames(LP_com))
want3<-grep("s\\(time\\)\\:session3",colnames(LP_com))
sim_com1 <- LP_com[idx1,want1] %*% t(B_com[,want1])
sim_com2 <- LP_com[idx2,want2] %*% t(B_com[,want2])
sim_com3 <- LP_com[idx3,want3] %*% t(B_com[,want3])
sim_com <- rbind(sim_com1, sim_com2, sim_com3)
sim_comd <- rbind(sim_com1-sim_com2, sim_com1-sim_com3, sim_com2-sim_com3)

sd_com <- abs(sweep(sim_com, 1, se_com, FUN = "/"))
se_comd <- rep(0,nrow(sd_com))
sd_comd <- sd_com
for(i in 1:nrow(sd_com)){
  se_comd[i] <- sd(sim_comd[i,])
  sd_comd[i,] <- sim_comd[i,]/se_comd[i]
  }
  
msd_com <- apply(sd_com, 2L, max)
msd_comd <- apply(sd_comd, 2L, max)

crit_com <- quantile(msd_com, prob = 0.95, type = 8)
crit_comd <- quantile(msd_comd, prob = 0.95, type = 8)

#plot fits
comdf <- data.frame(fit=fit_com, low=fit_com-crit_com*se_com, high=fit_com+crit_com*se_com)
comdf$session <- pdat$session
comdf$time <- pdat$time
# "#669966"
# "#CCFFCC"

# "gray"
# "gray95"

ggplot(comdf, aes(x = time, y = fit)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray") +
  geom_line() +
  geom_hline(aes(yintercept=0),linetype="dashed")+
  scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),name="time (min)",expand=c(0,0))+
  scale_y_continuous(name="COM (cm)",expand=c(0,0))+
  facet_wrap(~session,nrow=1)+
  theme(panel.background = element_rect(color="black",fill="white",size=1),aspect.ratio = 1,
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        panel.spacing.x=unit(0.2,"lines"),
        plot.background = element_rect(fill = "gray95",color="gray95"),
        panel.border = element_rect(colour="black",fill=NA,size=1))

ggsave("COM_Hyb.eps",device="eps",width=4,height=4,units="in")

#plot fit differences
comdf <- data.frame(fit=fit_comd, low=fit_comd-crit_comd*se_comd, high=fit_comd+crit_comd*se_comd)
comdf$session <- pdat$session
comdf$time <- pdat$time

ggplot(comdf, aes(x = time, y = fit)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray") +
  geom_line() +
  geom_hline(aes(yintercept=0),linetype="dashed")+
  scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),name="time (min)",expand=c(0,0))+
  scale_y_continuous(name="delta COM (cm)",expand=c(0,0),limits=c(-4,8))+
  facet_wrap(~session,nrow=1)+
  theme(panel.background = element_rect(color="black",fill="white",size=1),aspect.ratio = 1,
        strip.background=element_rect(color="grey50",fill="white",size=0.5),
        panel.spacing.x=unit(0.2,"lines"),
        plot.background = element_rect(fill = "gray95",color="gray95"),
        panel.border = element_rect(colour="black",fill=NA,size=1))

ggsave("dCOM_Hyb.eps",device="eps",width=4,height=4,units="in")


###### MODES ######

#run gam 
mode <- gam(MODES ~ SESSION + s(TIME,by=SESSION,id=1), data = data, method = "REML")

#plot(speed,rug=FALSE,scale=0)

#make prediciton data
nb <- 250
mint <- 0
maxt <- 600
tbins <- seq(mint,maxt,length.out=nb)
sessions <- unique(data$SESSION)
ns <- length(sessions)
pdat<-with(data,data.frame(TIME=rep(tbins,times=ns),SESSION=rep(sessions,each=nb)))

#simulate simultaneous intervals
nsim <- 10000
vc_mode <- vcov(mode)

ptype <- "iterms"
pred_mode <- predict(mode, pdat, type = ptype, se.fit = TRUE)

se_mode <- pred_mode$se.fit
fit_mode <- pred_mode$fit
idx1 <- which(pdat$SESSION == "1")
idx2 <- which(pdat$SESSION == "2")
idx3 <- which(pdat$SESSION == "3")
want1<-grep("s\\(TIME\\)\\:SESSION1",colnames(se_mode))
want2<-grep("s\\(TIME\\)\\:SESSION2",colnames(se_mode))
want3<-grep("s\\(TIME\\)\\:SESSION3",colnames(se_mode))
se_mode <- c(se_mode[idx1,want1], se_mode[idx2,want2], se_mode[idx3,want3])
fit_moded <- c(fit_mode[idx1,want1]-fit_mode[idx2,want2], 
              fit_mode[idx1,want1]-fit_mode[idx3,want3], 
              fit_mode[idx2,want2]-fit_mode[idx3,want3])
fit_mode <- c(fit_mode[idx1,want1], fit_mode[idx2,want2], fit_mode[idx3,want3])

set.seed(42)
B_mode <- mvrnorm(nsim, mu = rep(0, nrow(vc_mode)), Sigma = vc_mode)

LP_mode <- predict(mode, pdat, type = "lpmatrix")

want1<-grep("s\\(TIME\\)\\:SESSION1",colnames(LP_mode))
want2<-grep("s\\(TIME\\)\\:SESSION2",colnames(LP_mode))
want3<-grep("s\\(TIME\\)\\:SESSION3",colnames(LP_mode))
sim_mode1 <- LP_mode[idx1,want1] %*% t(B_mode[,want1])
sim_mode2 <- LP_mode[idx2,want2] %*% t(B_mode[,want2])
sim_mode3 <- LP_mode[idx3,want3] %*% t(B_mode[,want3])
sim_mode <- rbind(sim_mode1, sim_mode2, sim_mode3)
sim_moded <- rbind(sim_mode1-sim_mode2, sim_mode1-sim_mode3, sim_mode2-sim_mode3)

sd_mode <- abs(sweep(sim_mode, 1, se_mode, FUN = "/"))
se_moded <- rep(0,nrow(sd_mode))
sd_moded <- sd_mode
for(i in 1:nrow(sd_mode)){
  se_moded[i] <- sd(sim_moded[i,])
  sd_moded[i,] <- sim_moded[i,]/se_moded[i]
}

msd_mode <- apply(sd_mode, 2L, max)
msd_moded <- apply(sd_moded, 2L, max)

crit_mode <- quantile(msd_mode, prob = 0.95, type = 8)
crit_moded <- quantile(msd_moded, prob = 0.95, type = 8)

#plot fits
modedf <- data.frame(fit=fit_mode, low=fit_mode-crit_mode*se_mode, high=fit_mode+crit_mode*se_mode)
modedf$session <- pdat$SESSION
modedf$time <- pdat$TIME

ggplot(modedf, aes(x = time, y = fit)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray") +
  geom_line() +
  geom_hline(aes(yintercept=0),linetype="dashed")+
  scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),name="time (min)",expand=c(0,0))+
  scale_y_continuous(name="RI",expand=c(0,0))+
  facet_wrap(~session,nrow=1)+theme(aspect.ratio = 1)+
  theme(aspect.ratio = 1,
        panel.background = element_blank(),aspect.ratio = 1,
        panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(colour="black",fill=NA,size=1))

#plot fit differences
modedf <- data.frame(fit=fit_moded, low=fit_moded-crit_moded*se_moded, high=fit_moded+crit_moded*se_moded)
modedf$session <- pdat$SESSION
modedf$time <- pdat$TIME

ggplot(modedf, aes(x = time, y = fit)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray") +
  geom_line() +
  geom_hline(aes(yintercept=0),linetype="dashed")+
  scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),name="time (min)",expand=c(0,0))+
  scale_y_continuous(name="delta RI",expand=c(0,0))+
  facet_wrap(~session,nrow=1)+theme(aspect.ratio = 1)+
  theme(aspect.ratio = 1,
        panel.background = element_blank(),aspect.ratio = 1,
        panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(colour="black",fill=NA,size=1))


###### ZRATES ######

#run gam 
zrate <- gam(RATES ~ SESSION + s(TIME,by=SESSION,id=1), data = data, method = "REML")

#plot(speed,rug=FALSE,scale=0)

#make prediciton data
nb <- 250
mint <- 0
maxt <- 600
tbins <- seq(mint,maxt,length.out=nb)
sessions <- unique(data$SESSION)
ns <- length(sessions)
pdat<-with(data,data.frame(TIME=rep(tbins,times=ns),SESSION=rep(sessions,each=nb)))

#simulate simultaneous intervals
nsim <- 10000
vc_zrate <- vcov(zrate)

ptype <- "iterms"
pred_zrate <- predict(zrate, pdat, type = ptype, se.fit = TRUE)

se_zrate <- pred_zrate$se.fit
fit_zrate <- pred_zrate$fit
idx1 <- which(pdat$SESSION == "1")
idx2 <- which(pdat$SESSION == "2")
idx3 <- which(pdat$SESSION == "3")
want1<-grep("s\\(TIME\\)\\:SESSION1",colnames(se_zrate))
want2<-grep("s\\(TIME\\)\\:SESSION2",colnames(se_zrate))
want3<-grep("s\\(TIME\\)\\:SESSION3",colnames(se_zrate))
se_zrate <- c(se_zrate[idx1,want1], se_zrate[idx2,want2], se_zrate[idx3,want3])
fit_zrated <- c(fit_zrate[idx1,want1]-fit_zrate[idx2,want2], 
              fit_zrate[idx1,want1]-fit_zrate[idx3,want3], 
              fit_zrate[idx2,want2]-fit_zrate[idx3,want3])
fit_zrate <- c(fit_zrate[idx1,want1], fit_zrate[idx2,want2], fit_zrate[idx3,want3])

set.seed(42)
B_zrate <- mvrnorm(nsim, mu = rep(0, nrow(vc_zrate)), Sigma = vc_zrate)

LP_zrate <- predict(zrate, pdat, type = "lpmatrix")

want1<-grep("s\\(TIME\\)\\:SESSION1",colnames(LP_zrate))
want2<-grep("s\\(TIME\\)\\:SESSION2",colnames(LP_zrate))
want3<-grep("s\\(TIME\\)\\:SESSION3",colnames(LP_zrate))
sim_zrate1 <- LP_zrate[idx1,want1] %*% t(B_zrate[,want1])
sim_zrate2 <- LP_zrate[idx2,want2] %*% t(B_zrate[,want2])
sim_zrate3 <- LP_zrate[idx3,want3] %*% t(B_zrate[,want3])
sim_zrate <- rbind(sim_zrate1, sim_zrate2, sim_zrate3)
sim_zrated <- rbind(sim_zrate1-sim_zrate2, sim_zrate1-sim_zrate3, sim_zrate2-sim_zrate3)

sd_zrate <- abs(sweep(sim_zrate, 1, se_zrate, FUN = "/"))
se_zrated <- rep(0,nrow(sd_zrate))
sd_zrated <- sd_zrate
for(i in 1:nrow(sd_zrate)){
  se_zrated[i] <- sd(sim_zrated[i,])
  sd_zrated[i,] <- sim_zrated[i,]/se_zrated[i]
}

msd_zrate <- apply(sd_zrate, 2L, max)
msd_zrated <- apply(sd_zrated, 2L, max)

crit_zrate <- quantile(msd_zrate, prob = 0.95, type = 8)
crit_zrated <- quantile(msd_zrated, prob = 0.95, type = 8)

#plot fits
zratedf <- data.frame(fit=fit_zrate, low=fit_zrate-crit_zrate*se_zrate, high=fit_zrate+crit_zrate*se_zrate)
zratedf$session <- pdat$SESSION
zratedf$time <- pdat$TIME

ggplot(zratedf, aes(x = time, y = fit)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray") +
  geom_line() +
  geom_hline(aes(yintercept=0),linetype="dashed")+
  scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),name="time (min)",expand=c(0,0))+
  scale_y_continuous(name="firing rate (z-score)",expand=c(0,0))+
  facet_wrap(~session,nrow=1)+theme(aspect.ratio = 1)+
  theme(panel.background = element_blank(),aspect.ratio = 1,
        panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(colour="black",fill=NA,size=1))

#plot fit differences
zratedf <- data.frame(fit=fit_zrated, low=fit_zrated-crit_zrated*se_zrated, high=fit_zrated+crit_zrated*se_zrated)
zratedf$session <- pdat$SESSION
zratedf$time <- pdat$TIME

ggplot(zratedf, aes(x = time, y = fit)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray") +
  geom_line() +
  geom_hline(aes(yintercept=0),linetype="dashed")+
  scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),name="time (min)",expand=c(0,0))+
  scale_y_continuous(name="firing rate (delta zscore)",expand=c(0,0))+
  facet_wrap(~session,nrow=1)+theme(aspect.ratio = 1)+
  theme(aspect.ratio = 1,
        panel.background = element_blank(),aspect.ratio = 1,
        panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(colour="black",fill=NA,size=1))


###### SPEED ######

#run gam 
speed <- gam(SPEED ~ SESSION + s(TIME,by=SESSION,id=1), data = data, method = "REML")

#plot(speed,rug=FALSE,scale=0)

#make prediciton data
nb <- 250
mint <- 0
maxt <- 600
tbins <- seq(mint,maxt,length.out=nb)
sessions <- unique(data$SESSION)
ns <- length(sessions)
pdat<-with(data,data.frame(TIME=rep(tbins,times=ns),SESSION=rep(sessions,each=nb)))

#simulate simultaneous intervals
nsim <- 10000
vc_speed <- vcov(speed)

ptype <- "iterms"
pred_speed <- predict(speed, pdat, type = ptype, se.fit = TRUE)

se_speed <- pred_speed$se.fit
fit_speed <- pred_speed$fit
idx1 <- which(pdat$SESSION == "1")
idx2 <- which(pdat$SESSION == "2")
idx3 <- which(pdat$SESSION == "3")
want1<-grep("s\\(TIME\\)\\:SESSION1",colnames(se_speed))
want2<-grep("s\\(TIME\\)\\:SESSION2",colnames(se_speed))
want3<-grep("s\\(TIME\\)\\:SESSION3",colnames(se_speed))
se_speed <- c(se_speed[idx1,want1], se_speed[idx2,want2], se_speed[idx3,want3])
fit_speedd <- c(fit_speed[idx1,want1]-fit_speed[idx2,want2], 
              fit_speed[idx1,want1]-fit_speed[idx3,want3], 
              fit_speed[idx2,want2]-fit_speed[idx3,want3])
fit_speed <- c(fit_speed[idx1,want1], fit_speed[idx2,want2], fit_speed[idx3,want3])

set.seed(42)
B_speed <- mvrnorm(nsim, mu = rep(0, nrow(vc_speed)), Sigma = vc_speed)

LP_speed <- predict(speed, pdat, type = "lpmatrix")

want1<-grep("s\\(TIME\\)\\:SESSION1",colnames(LP_speed))
want2<-grep("s\\(TIME\\)\\:SESSION2",colnames(LP_speed))
want3<-grep("s\\(TIME\\)\\:SESSION3",colnames(LP_speed))
sim_speed1 <- LP_speed[idx1,want1] %*% t(B_speed[,want1])
sim_speed2 <- LP_speed[idx2,want2] %*% t(B_speed[,want2])
sim_speed3 <- LP_speed[idx3,want3] %*% t(B_speed[,want3])
sim_speed <- rbind(sim_speed1, sim_speed2, sim_speed3)
sim_speedd <- rbind(sim_speed1-sim_speed2, sim_speed1-sim_speed3, sim_speed2-sim_speed3)

sd_speed <- abs(sweep(sim_speed, 1, se_speed, FUN = "/"))
se_speedd <- rep(0,nrow(sd_speed))
sd_speedd <- sd_speed
for(i in 1:nrow(sd_speed)){
  se_speedd[i] <- sd(sim_speedd[i,])
  sd_speedd[i,] <- sim_speedd[i,]/se_speedd[i]
}

msd_speed <- apply(sd_speed, 2L, max)
msd_speedd <- apply(sd_speedd, 2L, max)

crit_speed <- quantile(msd_speed, prob = 0.95, type = 8)
crit_speedd <- quantile(msd_speedd, prob = 0.95, type = 8)

#plot fits
speeddf <- data.frame(fit=fit_speed, low=fit_speed-crit_speed*se_speed, high=fit_speed+crit_speed*se_speed)
speeddf$session <- pdat$SESSION
speeddf$time <- pdat$TIME

ggplot(speeddf, aes(x = time, y = fit)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray") +
  geom_line() +
  geom_hline(aes(yintercept=0),linetype="dashed")+
  scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),name="time (min)",expand=c(0,0))+
  scale_y_continuous(name="speed (cm/sec)",expand=c(0,0))+
  facet_wrap(~session,nrow=1)+theme(aspect.ratio = 1)+
  theme(panel.background = element_blank(),aspect.ratio = 1,
        panel.margin.x=unit(0,"lines"),panel.margin.y=unit(0,"lines"),
        panel.border = element_rect(colour="black",fill=NA,size=1))

#plot fit differences
speeddf <- data.frame(fit=fit_speedd, low=fit_speedd-crit_speedd*se_speedd, high=fit_speedd+crit_speedd*se_speedd)
speeddf$session <- pdat$SESSION
speeddf$time <- pdat$TIME

ggplot(speeddf, aes(x = time, y = fit)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray") +
  geom_line() +
  geom_hline(aes(yintercept=0),linetype="dashed")+
  scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),name="time (min)",expand=c(0,0))+
  scale_y_continuous(name="delta speed (cm/sec)",expand=c(0,0))+
  facet_wrap(~session,nrow=1)+theme(aspect.ratio = 1)+
  theme(panel.background = element_blank(),aspect.ratio = 1,
        panel.margin.x=unit(-0.02,"lines"),panel.margin.y=unit(-0.02,"lines"),
        panel.border = element_rect(colour="black",fill=NA,size=1))



ggplot(data = data,aes(x=time,y=medp)) +
  geom_density2d() +
  geom_smooth(method="gam",formula=y~s(x)) +
  facet_wrap(~session,nrow=1)





