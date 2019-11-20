#compute simultaneous confidence intervals from simulations using line search 
#for pointwise error rate that controls familywise error rate
simCIz<-function(sims)
{

  alpha<-0.05 #target alpha
  oalpha<-1 #arbitrary starting value for observed alpha
  zlo<-3  #lower bound for line search
  zhi<-7  #upper bound for line search (may need to go higher for many comparisons)
  z_new<-5  #starting value for line search
  gr<-0.382 #split by golden ratio
  eps<-0.001 #convergence epsilon
  temp<-matrix(0,nrow = nrow(sims),ncol = 3) #matrix to hold results
  
  #z-score across smiulations
  sims_sd<-rep(0,times=nrow(sims))
  for(i in 1:nrow(sims)){
    temp[i,2]<-mean(sims[i,])
    sims_sd[i]<-sd(sims[i,])
    # sims[i,]<-(sims[i,]-temp[i,2])/sims_sd[i]
    sims[i,]<-abs(sims[i,]-temp[i,2])/sims_sd[i]
  }
  
  maxima <- apply(sims,2,max)
  z_old <- quantile(maxima,prob=0.95,type=8)
  
  # while(abs(oalpha-alpha)>eps)
  # {
  #   obs<-rep(0,times=ncol(sims))
  #   for(i in 1:ncol(sims)){obs[i]<-sum(sims[,i] > z_new)}
  #   oalpha<-sum(obs>0)/ncol(sims)
  # 
  #   print(paste0("z = ",z_new))
  #   print(paste0("Observed alpha = ",oalpha))
  # 
  #   if(oalpha>alpha){
  #     zlo<-z_new
  #     z_old<-z_new
  #     z_new<-z_old + gr*(zhi-z_old)
  # 
  #   } else{
  #     zhi<-z_new
  #     z_old<-z_new
  #     z_new<-z_old - gr*(z_old-zlo)
  #   }
  # }
  temp[,1]<--z_old*sims_sd
  temp[,3]<-z_old*sims_sd
  return(temp)
}