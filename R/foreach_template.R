registerDoParallel(cores=4)
for(i in 1:5) {
  a<-i
  v<-foreach(f=1:10,.combine=rbind) %dopar% {
    r<-f
  }
  print(paste0(i,"_",v))
  
}