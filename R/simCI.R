#compute simultaneous confidence intervals from simulations using line search 
#for pointwise error rate that controls familywise error rate
simCI<-function(sims)
{
  alpha<-0.05
  oalpha<-1
  slo<-1
  shi<-10
  scale_new<-1
  gr<-0.382 #split by golden ratio
  eps<-0.001
  temp<-matrix(0,nrow = nrow(sims),ncol = 3)
  delta<-1
  
  sims<-t(scale(t(sims),center=TRUE,scale=FALSE))
  
  while(abs(oalpha-alpha)>eps && delta>.Machine$double.eps)
  {
    obs<-matrix(0,nrow = 1,ncol = ncol(sims))
    for(i in 1:nrow(sims))
    {
      temp[i,]<-quantile(scale_new*sims[i,],probs=c(alpha,0.5,1-alpha))
      obs<-obs+(sims[i,]<temp[i,1] | sims[i,]>temp[i,3])
    }
    oalpha<-sum(obs>0)/ncol(sims)

    print(paste0("Scale = ",scale_new))
    print(paste0("Observed alpha = ",oalpha))
    
    if(oalpha>alpha){
      slo<-scale_new  
      scale_old<-scale_new
      scale_new<-scale_old + gr*(shi-scale_old)
      
    } else{
      shi<-scale_new
      scale_old<-scale_new
      scale_new<-scale_old - gr*(scale_old-slo)
    }
    delta<-abs(scale_new-scale_old)
  }
  return(temp)
}