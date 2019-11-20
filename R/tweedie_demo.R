library(mgcv)
library(ggplot2)

mu <- 3
p <- c(1.0001,1.001,1.01,1.1,1.2,1.7,1.99)
phi <- 1
nb <- 60
df <- data.frame(mu = rep(mu,times=nb*length(p)),
                 p = rep(p,each=nb),
                 phi = rep(phi,times=nb*length(p)),
                 mids = rep(0,times=nb*length(p)),
                 counts = rep(0,times=nb*length(p)))

for (i in 1:length(p)) {
  y <- rTweedie(mu=rep(mu,times=100000),p=p[i],phi=phi)
  h <- hist(y[y<10],breaks=seq(0,10,length.out=nb+1),plot=FALSE)
  df$counts[((i-1)*nb+1):(i*nb)] <- h$counts/sum(h$counts)
  df$mids[((i-1)*nb+1):(i*nb)] <- h$mids
}

ggplot(df,aes(x=mids,y=counts))+ 
  geom_col(width=10/60)+
  scale_y_continuous(expand=c(0,0))+
  facet_wrap(~p,nrow=1)+theme(aspect.ratio = 1)+
  guides(fill = guide_colorbar(title=NULL,label.theme = element_text(color="black",size=20,angle=0))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # strip.background=element_rect(color="grey50",fill="white",size=0.5),
        aspect.ratio = 1,
        panel.border = element_rect(colour="black",fill=NA))

ggsave("tweedie_demo.eps",device="eps",width=8,height=8,units="in")