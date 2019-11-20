library(mgcv)
library(MASS)
library(ggplot2)
library(itsadug)

setwd("E:/For Brian/PCtable/")
# setwd("G:/Rats/PCtable/")

#load data
data <- read.table("PCtable.csv",sep=",",header=TRUE)
data <- data[data$mouse>3,]

#set appropriate variables to factors and scale units
C <- 100*pi
data$cell <- factor(data$cell)
data$day <- factor(data$day) 
data$mouse <- factor(data$mouse)
data$session <- factor(data$session) 
data$MouseSession<-with(data,interaction(mouse,session))
data$com <- data$com/(2*pi)*C #set units to cm
data$median <- data$median/(2*pi)*C #set units to cm
data$peakpos <- data$peakpos/(2*pi)*C #set units to cm
data$firstspk <- data$firstspk/(2*pi)*C #set units to cm
data$lastspk <- data$lastspk/(2*pi)*C #set units to cm
data$lapwidth <- data$lapwidth/(2*pi)*C #set units to cm
data$rewdist <- data$rewdist/(2*pi)*C #set units to cm
data$velenter <- data$velenter/(2*pi)*C #set units to cm
data$velexit <- data$velexit/(2*pi)*C #set units to cm
data$speed <- (data$velenter*data$tenter + data$velexit*data$texit)/(data$tenter+data$texit)
data$rates <- data$numspks/(data$tenter+data$texit)
data$BAI<-(data$velenter-data$velexit)/(data$velenter+data$velexit)
data$mmdiff<-data$com-data$median

#normalize and filter
cells <- unique(data$cell)
nc <- length(cells)
data$zrates <- data$rates
data$murates <- data$rates
data$cnumspks <- data$numspks
for(c in 1:nc){
  mu <- mean(data$rates[data$cell==cells[c]],na.rm=TRUE)
  std <- sqrt(sum((data$rates[data$cell==cells[c]]-mu)^2))
  data$murates[data$cell==cells[c]] <- mu
  # data$zrates[data$cell==cells[c]] <- (data$rates[data$cell==cells[c]] - mu)/std
  data$zrates[data$cell==cells[c]] <- data$rates[data$cell==cells[c]]/mu
  data$firstspk[data$cell==cells[c]]<-data$firstspk[data$cell==cells[c]]-mean(data$firstspk[data$cell==cells[c]],na.rm=TRUE)
  data$lastspk[data$cell==cells[c]]<-data$lastspk[data$cell==cells[c]]-mean(data$lastspk[data$cell==cells[c]],na.rm=TRUE)
  data$com[data$cell==cells[c]]<-data$com[data$cell==cells[c]]-mean(data$com[data$cell==cells[c]],na.rm=TRUE)
  data$median[data$cell==cells[c]]<-data$median[data$cell==cells[c]]-mean(data$median[data$cell==cells[c]],na.rm=TRUE)
  data$mmdiff[data$cell==cells[c]]<-data$mmdiff[data$cell==cells[c]]-mean(data$mmdiff[data$cell==cells[c]],na.rm=TRUE)
  data$peakpos[data$cell==cells[c]]<-data$peakpos[data$cell==cells[c]]-mean(data$peakpos[data$cell==cells[c]],na.rm=TRUE)
  data$lapwidth[data$cell==cells[c]]<-data$lapwidth[data$cell==cells[c]]/mean(data$lapwidth[data$cell==cells[c]],na.rm=TRUE)
  data$fsize[data$cell==cells[c]]<-data$fsize[data$cell==cells[c]]/mean(data$fsize[data$cell==cells[c]],na.rm=TRUE)
  data$sd[data$cell==cells[c]]<-data$sd[data$cell==cells[c]]/mean(data$sd[data$cell==cells[c]],na.rm=TRUE)
  data$speed[data$cell==cells[c]]<-data$speed[data$cell==cells[c]]/mean(data$speed[data$cell==cells[c]],na.rm=TRUE)
  data$cnumspks[data$cell==cells[c]]<-sum(data$numspks[data$cell==cells[c]])
}
data <- data[data$cnumspks>300&data$cnumspks<6000,]
# data <- data[data$murates>2.5,]
data <- data[data$velenter>5,]
data <- data[data$velexit>5,]
# data <- data[data$rewdist>10,]
# data <- data[data$numspks>1,]

#run gams 
# data$session <- factor(data$session,ordered=FALSE,levels=c("1","2","3"))
com <- gam(com ~ session + s(time,by=session,id=1,k=20), data = data, method = "REML")

data1323 <- data
data1323$session <- ordered(data1323$session,levels=c("3","2","1"))
contrasts(data1323$session) <- 'contr.treatment'

data1232 <- data
data1232$session <- ordered(data1232$session,levels=c("2","3","1"))
contrasts(data1323$session) <- 'contr.treatment'

mod1323 <- gam(lastspk ~ session + s(time,k=20,id=1) + s(time,by=session,k=20,id=1), data = data1323, method = "REML")
mod1232 <- gam(lastspk ~ session + s(time,k=20,id=1) + s(time,by=session,k=20,id=1), data = data1232, method = "REML")


#make prediciton data
nb <- 250
mint <- 0
maxt <- 600
tbins <- seq(mint,maxt,length.out=nb)
data$session <- factor(data$session,ordered=FALSE,levels=c("1","2","3"))
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
  scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),name="time (sec)",expand=c(0,0))+
  scale_y_continuous(name="field width",expand=c(0,0))+
  facet_wrap(~session,nrow=1)+
  theme(aspect.ratio = 1,
        panel.spacing.x=unit(0.2,"lines"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA,size=1))
  

ggsave("fieldwidth.eps",device="eps",width=4,height=4,units="in")

#plot fit differences
comdf <- data.frame(fit=fit_comd, low=fit_comd-crit_comd*se_comd, high=fit_comd+crit_comd*se_comd)
comdf$session <- pdat$session
comdf$time <- pdat$time

ggplot(comdf, aes(x = time, y = fit)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "gray") +
  geom_line() +
  geom_hline(aes(yintercept=0),linetype="dashed")+
  scale_x_continuous(breaks=seq(0,600,by=300),labels=seq(0,10,by=5),name="time (sec)",expand=c(0,0))+
  scale_y_continuous(name="dfieldwidth",expand=c(0,0))+
  facet_wrap(~session,nrow=1)+
  theme(aspect.ratio = 1,
        panel.spacing.x=unit(0.2,"lines"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA,size=1))

ggsave("dfieldwidth.eps",device="eps",width=4,height=4,units="in")