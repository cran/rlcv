#' Robust likelihood cross validation bandwidth for univariate densities
#'
#' @param x.obs Training (observed) data
#' @param x.new Evaluation data; default to x.obs
#'
#' @return fhat: density evaluated at x.new; h: bandwidth
#' @export
#'
#' @author Ximing Wu \email{xwu@@tamu.edu}
#' @references Wu, Ximing (2019), "Robust Likelihood Cross Validation for Kernel Density Estimation," Journal of Business and Economic Statistics, 37(4): 761-770.
#'
#' @examples
#' x=rt(200,df=5)
#' x.new=seq(-5,5,length=100)
#' fit=rlcv(x.obs=x,x.new=x.new)
#' # Mean squared errors
#' f0=dt(x.new,df=5)
#' mean((f0-fit$fhat)^2)
#'
#' matplot(x.new,cbind(f0,fit$fhat),type='l')
rlcv = function(x.obs, x.new=NULL)
{
  # x.obs: data
  # x.new: evaluation points

  # gaussian quadrature; need package 'statmod'
  g=statmod::gauss.quad(300,kind='hermite')
  w=exp(g$n^2+log(g$w))
  g$weights=w

  sdx=stats::sd(x.obs)
  mux=mean(x.obs)
  y=(x.obs-mux)/sdx

  # initial value for bandwidth
  h0=stats::bw.SJ(y)

  a=an(length(x.obs),1)

  fit=stats::optimize(psi.g,c(h0/3,h0*3),x=y,a=a,g=g,maximum=T)
  h=(fit$max)*sdx

  if (is.null(x.new))
    x.new=x.obs

  fhat=kde(x.new=x.new,x.obs=x.obs,h)

  return(list(h=h,fhat=fhat))

}

an=function(n,d) gamma(d/2)/(2*pi)^(d/2)/(log(n))^(d/2-1)/n

psi.g = function(h,x,a,g)
{
  f=kde(x.new=g$n,x.obs=x,h)
  f=ifelse(f>=a,f,f^2/2/a)

  b_psi=sum(f*g$w)

  f_i=kde_i(x,h)

  lf_i=ifelse(f_i>=a,log(f_i),log(a)-1+f_i/a)
  return(mean(lf_i)-b_psi)
}


