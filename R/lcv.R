#' Likelihood cross validation bandwidth for univariate densities
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
#' fit=lcv(x.obs=x,x.new=x.new)
#' # Mean squared errors
#' f0=dt(x.new,df=5)
#' mean((f0-fit$fhat)^2)
#'
#' matplot(x.new,cbind(f0,fit$fhat),type='l')
lcv = function(x.obs, x.new=NULL)
{
  # x.obs: data
  # x.new: evaluation points


  sdx=stats::sd(x.obs)
  mux=mean(x.obs)
  y=(x.obs-mux)/sdx

  # initial value for bandwidth
  h0=stats::bw.SJ(y)


  fit=stats::optimize(logf,c(h0/3,h0*3),x=y,maximum=T)
  h=(fit$max)*sdx

  if (is.null(x.new))
    x.new=x.obs

  fhat=kde(x.new=x.new,x.obs=x.obs,h)

  return(list(h=h,fhat=fhat))

}

logf = function(h,x)
{
  f_i=kde_i(x,h)

  lf_i=log(f_i)
  return(mean(lf_i))
}


