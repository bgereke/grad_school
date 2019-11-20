setwd("C:/Users/Brian/Documents/R")

#load libraries  
library(mgcv)
library(data.table)

mydata <- fread("inputR.csv",sep=",",header=TRUE)
mod <- bam(x ~ s(phase,bs="cc",k=30) +
             ti(vel,phase,bs=c("cr","cc"),k=c(6,15))+
             ti(time,phase,bs=c("cr","cc"),k=c(6,15)),
             data = mydata)

pred <- predict.bam(mod,type="terms",newdata.guaranteed = TRUE)
fwrite(as.list(as.data.frame(pred)),file="outputR.csv")