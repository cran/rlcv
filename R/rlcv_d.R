#' Robust likelihood cross validation bandwidth for multivariate kernel densities
#'
#' @param x.obs Training (observed) data (n1 by d matrix, d>=2)
#' @param x.new Evaluation data (n2 by d matrix, d>=2); default to x.obs
#'
#' @return fhat: density evaluated at x.new; h: bandwidth
#' @export
#'
#' @author Ximing Wu \email{xwu@@tamu.edu}
#' @references Wu, Ximing (2019), "Robust Likelihood Cross Validation for Kernel Density Estimation," Journal of Business and Economic Statistics, 37(4): 761-770.
#'
#' @examples
#' # old faithful data
#' x=datasets::faithful
#' x=cbind(x[,1],x[,2])
#' fit=rlcv_d(x.obs=x)
#' # evaluation data
#' x1=seq(min(x[,1])*.8,max(x[,1])*1.2,length=30)
#' x2=seq(min(x[,2])*.8,max(x[,2])*1.2,length=30)
#' x11=rep(x1,each=30)
#' x22=rep(x2,30)
#' fhat=kde_d(x.new=cbind(x11,x22),x.obs=x,h=fit$h)
#' persp(x1,x2,matrix(fhat,30,30))
rlcv_d=function(x.obs,x.new=NULL)
{
  # robust likelihood cv
  # note the bandwidth is for studentized x
  n=dim(x.obs)[1]
  d=dim(x.obs)[2]

  a=an(n,d)

  mu=apply(x.obs,2,mean)

  S=stats::cov(x.obs)
  S.eig=eigen(S)
  S.inv.root=S.eig$vectors%*%diag(1/sqrt(S.eig$values))%*%t(S.eig$vectors)
  y=x.obs
  if (is.null(x.new))
    x.new=x.obs
  for (i in 1:d)
  {
    y[,i]=y[,i]-mu[i]
    x.new[,i]=x.new[,i]-mu[i]
  }

  y=y%*%S.inv.root
  x.new=x.new%*%S.inv.root

  h0=NULL
  for (i in 1:d)
    h0=c(h0,stats::bw.SJ(y[,i]))

  fit=stats::optim(par=h0,fn=psi_d.g,x=y,a=a,lower=h0/3,upper=h0*3,method='L-BFGS-B')
  h=fit[[1]]

  f=kde_d(x.new=x.new,x.obs=y,h,stud=TRUE)
  f=f/prod(sqrt(S.eig$values))


  return(list(fhat=f,h=h))

}


an=function(n,d) gamma(d/2)/(2*pi)^(d/2)/(log(n))^(d/2-1)/n

psi_d.g=function(h,x,a)
{

  f_i=kde_d_i(x,h)

  lf_i=ifelse(f_i>=a,log(f_i),log(a)-1+f_i/a)

  f=kde_d(x.new=x,x.obs=x,h=h,stud=TRUE)
  p=mean(f>=a)
  p2=mean(f*(f<a))
  -mean(lf_i)+p+p2/2/a
}
