library(mgcv)
library(ggplot2)

setwd("C:/Data/PCtable/")

#preallocate variables
nf <- 50
ns <- 12
freq<-10^seq(log10(2),log10(100),length.out=nf)
statsdf <- data.frame(Term = rep(0,times=nf*ns),
                      Frequency = rep(freq,times=ns),
                      Fval = rep(0,times=nf*ns),
                      pval = rep(0,times=nf*ns),
                      edf = rep(0,times=nf*ns),
                      refedf = rep(0,times=nf*ns))

#gather stats
for(f in 1:nf){
  
  filename <- paste0("FULLarho_f",f,"_POW.rds")
  fmod <- readRDS(filename)
  sum <- summary(fmod)$s.table
  sumdf <- as.data.frame(sum)
  
  idx <- seq(f,nf*ns,by=nf)
  
  statsdf$Term[idx] <- rownames(sum)[1:ns]
  statsdf$Fval[idx] <- sumdf$F[1:ns]
  statsdf$pval[idx] <- sumdf$"p-value"[1:ns]
  statsdf$edf[idx] <- sumdf$edf[1:ns]
  statsdf$refedf[idx] <- sumdf$"Ref.df"[1:ns]
  
  
  print(paste("Completed: ",f," of ",nf),quote=FALSE)
}

statsdf$Term = factor(statsdf$Term,ordered=TRUE,
                         levels=unique(statsdf$Term))

#make plots
ggplot(data = statsdf,aes(x=Frequency,y=pval)) +
  geom_line() + 
  scale_y_log10() +
  geom_hline(aes(yintercept=10^-5 ),linetype="dashed") +
  facet_wrap(~Term,nrow=4,scales="free")

