#' Multivariate kernel density
#'
#' @param x.obs Training (observed) data (n1 by d matrix, d>=2)
#' @param x.new Evaluation data (n2 by d matrix, d>=2); default to x.obs
#' @param h Bandwidth (d vector)
#' @param stud Indicator for whether data are studentized; default to FALSE
#'
#' @return Density evaluated at x.new
#'
#' @details For multivariate distributions, bandwidth is calculated for studentized data.
#' @export
#'
#' @author Ximing Wu \email{xwu@@tamu.edu}
#' @references Wu, Ximing (2019), "Robust Likelihood Cross Validation for Kernel Density Estimation," Journal of Business and Economic Statistics, 37(4): 761-770.
#'
#' @examples
#' x=matrix(rnorm(200),ncol=2)
#' x.new=matrix(rnorm(100),ncol=2)
#' h=c(1,1)
#' f=kde_d(x.new=x.new,x.obs=x,h=h)
kde_d=function(x.obs,x.new=NULL,h,stud=FALSE)
{
  # d-dimensional kde with product kernel
  # d>=2

  # x.new: eavluation points
  # x.obs: data
  # x.new and x.obs are assumed not studentized

  d=dim(x.obs)[2]

  if (is.null(x.new))
    x.new=x.obs

  if (stud==FALSE)
  {
    S=stats::cov(x.obs)
    S.eig=eigen(S)
    S.inv.root=S.eig$vectors%*%diag(1/sqrt(S.eig$values))%*%t(S.eig$vectors)
    x.obs=x.obs%*%S.inv.root
    x.new=x.new%*%S.inv.root
  }

  f = outer(x.new[,1],x.obs[,1],function(p,q) stats::dnorm((p-q)/h[1])/h[1])
  for (i in 2:d)
    f=f*outer(x.new[,i],x.obs[,i],function(p,q) stats::dnorm((p-q)/h[i])/h[i])

  f=rowMeans(f)

  if (stud==FALSE)
    f=f/prod(sqrt(S.eig$values))
  return(f)
}

kde_d_i=function(x,h)
{
  # leave-one-out kde
  # this is an interval function called by rlcv_d
  # x is assumed to be studentized
  n=dim(x)[1]
  d=dim(x)[2]
  id=1:n
  id=outer(id,id,'-')
  id=matrix(id,ncol=1)

  x.i=outer(x[,1],x[,1],'-')
  x.i=matrix(x.i,ncol=1)
  x.i=x.i[id!=0]
  f.i=stats::dnorm(x.i/h[1])/h[1]

  for (i in 2:d)
  {
    x.i=outer(x[,i],x[,i],'-')
    x.i=matrix(x.i,ncol=1)
    x.i=x.i[id!=0]
    f.i=f.i*stats::dnorm(x.i/h[i])/h[i]

  }

  f.i=matrix(f.i,nrow=n-1)
  f.i=colMeans(f.i)

  return(f.i)
}
