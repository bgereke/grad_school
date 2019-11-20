#simple GAM simulation
x1 <- runif(1000)
x2 <- runif(1000)
x3 <- runif(1000)

y1 <- sin(2*pi*3*x1) - mean(sin(2*pi*3*x1))
y2 <- exp(x2) - mean(exp(x2))
y <- exp(5+y1+y2)
shape=2
yr <- rgamma(1000,shape=shape,scale=y/shape)
dat <- data.frame(y=yr,x1=x1,x2=x2,x3=x3)
modn <- gam(log(y)~s(x1,bs="cr",k=20)+s(x2,bs="cr",k=20)+s(x3,bs="cr",k=20),data=dat)
modg <- gam(y~s(x1,bs="cr",k=20)+s(x2,bs="cr",k=20)+s(x3,bs="cr",k=20),data=dat,family = Gamma(link="log"))


newdat <- data.frame(x1=seq(0,1,length.out=1000),
                     x2=seq(0,1,length.out=1000),
                     x3=seq(0,1,length.out=1000))
pn <- predict(modn,newdata=newdat,type="terms",se.fit=TRUE)
pg <- predict(modg,newdata=newdat,type="terms",se.fit=TRUE)
basis <- predict(modg,newdata=newdat,type="lpmatrix")


newdat$x1nfit <- pn$fit[,1]
newdat$x2nfit <- pn$fit[,2]
newdat$x3nfit <- pn$fit[,3]
newdat$x1gfit <- pg$fit[,1]
newdat$x2gfit <- pg$fit[,2]
newdat$x3gfit <- pg$fit[,3]

newdat$x1nhigh <- pn$fit[,1] + 1.96*pn$se.fit[,1]
newdat$x2nhigh <- pn$fit[,2] + 1.96*pn$se.fit[,2]
newdat$x3nhigh <- pn$fit[,3] + 1.96*pn$se.fit[,3]
newdat$x1ghigh <- pg$fit[,1] + 1.96*pg$se.fit[,1]
newdat$x2ghigh <- pg$fit[,2] + 1.96*pg$se.fit[,2]
newdat$x3ghigh <- pg$fit[,3] + 1.96*pg$se.fit[,3]

newdat$x1nlow <- pn$fit[,1] - 1.96*pn$se.fit[,1]
newdat$x2nlow  <- pn$fit[,2] - 1.96*pn$se.fit[,2]
newdat$x3nlow  <- pn$fit[,3] - 1.96*pn$se.fit[,3]
newdat$x1glow  <- pg$fit[,1] - 1.96*pg$se.fit[,1]
newdat$x2glow  <- pg$fit[,2] - 1.96*pg$se.fit[,2]
newdat$x3glow  <- pg$fit[,3] - 1.96*pg$se.fit[,3]

#make scatterplot
scatdf <- data.frame(x = rep(seq(0,1,length.out=100),times=3),
                     true = c())