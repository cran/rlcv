
#' Univariate kernel density
#'
#' @param x.obs Training (observed) data (n1 vector)
#' @param x.new Evaluation data (n2 vector); default to x.obs
#' @param h Bandwidth
#'
#' @return Density evaluated at x.new
#' @export
#'
#'
#' @author Ximing Wu \email{xwu@@tamu.edu}
#' @references Wu, Ximing (2019), "Robust Likelihood Cross Validation for Kernel Density Estimation," Journal of Business and Economic Statistics, 37(4): 761-770.
#'
#' @examples
#' x=rnorm(100)
#' x.new=seq(-5,5,length=50)
#' h=1.06*sd(x)*(length(x))^(-1/5)
#' f=kde(x.new=x.new,x.obs=x,h=h)
kde <- function(x.obs,x.new=NULL,h) {
  # x.new: evaluation data
  # x.obs: observations
  if (is.null(x.new))
    x.new=x.obs
  f=rowMeans(outer(x.new,x.obs,function(p,q) stats::dnorm((p-q)/h)/h))
  return(f)
}


kde_i=function(x,h)
{ # kde, leave one out
  n=length(x)
  df=outer(x,x,'-')
  id=1:n
  id1=outer(id,id,'-')
  id1=matrix(id1,ncol=1)
  df1=matrix(df,ncol=1)
  df2=df1[id1!=0]
  df2=matrix(df2,nrow=n-1)
  f=colMeans(stats::dnorm(df2/h)/h)
}

