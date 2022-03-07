## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(rlcv)
set.seed(12345)
x=rlnorm(300)
fit1=lcv(x.obs=x)
fit2=rlcv(x.obs=x)
c(fit1$h,fit2$h)

## ----fig.dim=c(6,4)-----------------------------------------------------------
x.new=seq(0,10,length=100)
f0=dlnorm(x.new)
f1=kde(x.obs=x,x.new=x.new,h=fit1$h)
f2=kde(x.obs=x,x.new=x.new,h=fit2$h)
# compare mean squared error
mean((f0-f1)^2)
mean((f0-f2)^2)
# density plots
matplot(x.new,cbind(f0,f1,f2),col=c("black","green","red"),type='l',xlab='x',ylab='')
legend('right',legend = c('true','LCV','RLCV'), col = c('black','green','red'), lty = 1:3, lwd = 1 , xpd = T )

## ----fig.dim=c(6,4)-----------------------------------------------------------
x=datasets::faithful
x=cbind(x[,1],x[,2])
fit1=lcv_d(x.obs=x)
fit2=rlcv_d(x.obs=x)
fit1$h
fit2$h
x1=seq(min(x[,1])*.8,max(x[,1])*1.2,length=50)
x2=seq(min(x[,2])*.8,max(x[,2])*1.2,length=50)
x11=rep(x1,each=50)
x22=rep(x2,50)
f1=kde_d(x.new=cbind(x11,x22),x.obs=x,h=fit1$h)
f2=kde_d(x.new=cbind(x11,x22),x.obs=x,h=fit2$h)
filled.contour(x1,x2,matrix(f1,50,50))
filled.contour(x1,x2,matrix(f2,50,50))

## ----fig.dim=c(6,4)-----------------------------------------------------------
library(copula)
set.seed(12345)
ncop=tCopula(.5,df=5)
n=500
u=rCopula(n,ncop)
x1=qt(u[,1],5)
x2=qt(u[,2],5)
x=cbind(x1,x2)
fit1=lcv_d(x.obs=x)
fit2=rlcv_d(x.obs=x)
fit1$h
fit2$h
# evaluation data
x1=x2=seq(-5,5,length=50)
x11=rep(x1,each=50)
x22=rep(x2,50)
f1=kde_d(x.new=cbind(x11,x22),x.obs=x,,h=fit1$h)
f2=kde_d(x.new=cbind(x11,x22),x.obs=x,,h=fit2$h)
f0=dCopula(cbind(pt(x11,5),pt(x22,5)),ncop)*dt(x11,5)*dt(x22,5)
# Mean squared errors
mean((f0-f1)^2)
mean((f0-f2)^2)
# check summary statistics
summary(cbind(f0,f1,f2))
# plot the fitted densities
filled.contour(x1,x2,matrix(f1,50,50))
filled.contour(x1,x2,matrix(f2,50,50))

