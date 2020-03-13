## definition of a mixture of Gaussian (Normal) distribution
## functions, following R's conventions: {r,p,d,q}mixnorm.

## GENERATE from a mixture of gaussians

rmixnorm <- function( size=1, p=c(1.0), mu=c(0.0), sigma=c(1.0) )
{
  # generate from sum_j p_j N(mu_j,sigma_j)
  #

  # some checks on the input
  #
  if ( abs(sum(p)-1.0)>1.0e-5 )   { stop( "p's must sum to 1" ) }
  if ( length(p)!=length(mu) )    { stop( "length(p)!=length(mu)" ) }
  if ( length(p)!=length(sigma) ) { stop( "length(p)!=length(sigma)" ) }

  J <- length(p)
  if ( J==1 ) {
    return( rnorm(size,mu[1],sigma[1]) )
  }
  z <- sample(x=J,size=size,replace=TRUE,prob=p)
  n <- tabulate(z,nbins=J)
  y <- rep( NA, size )
  for ( j in (1:J) )
  {
    y[z==j] <- rnorm( n[j], mu[j], sigma[j] )
  }
  return( y )
}
## DENSITY function
##
dmixnorm <- function( x, p=c(1.0), mu=c(0.0), sigma=c(1.0) )
{
  matrix(sapply(x, dnorm, mean=mu, sd=sigma), length(x), length(p), byrow=T) %*% p
}
## QUANTILE function
##
qmixnorm <- function( x, p=c(1.0), mu=c(0.0), sigma=c(1.0) )
{
  matrix(sapply(x, qnorm, mean=mu, sd=sigma), length(x), length(p), byrow=TRUE) %*% p
}
## PLOT function
##
plot.mixnorm <- function(n, p=c(1.0), mu=c(0.0), sigma=c(1.0), add=FALSE, type="l", return.ylim=FALSE,
                        xlim=NULL, ylim=NULL, seed=NULL, ... )
{
  if (!is.null(seed)) set.seed(seed)
  if ( length(mu)==1 && length(p)>1 )
    mu <- rep(mu,length(p))
  if ( length(sigma)==1 && length(p)>1 )
    sigma <- rep(sigma,length(p))

  Q <- sort(rmixnorm( size=n, p=p, mu=mu, sigma=sigma ))
  D <- dmixnorm(Q, p=p, mu=mu, sigma=sigma)
  xlim <- c( min(xlim[1],min(Q)), max(xlim[2],max(Q)) )
  ylim <- c( min(ylim[1],min(D)), max(ylim[2],max(D)) )
  #q <- qmixnorm((1:(n-1))/n, p=p, mu=mu, sigma=sigma)
  if ( add )
    lines( Q, D, type=type, ... )
  else
    plot( Q, D, type=type, xlim=xlim, ylim=ylim, ... )
  if (return.ylim) return( ylim )
}

###########################################
## MCLUST-related functions
###########################################

## FIT the model using mclust
##
mixnormfit <- function( x, G=1:5, model="V" )
{
  em <- summary(Mclust( x, G, model=model ), x)
  list( k=em$G, p=em$p, mu=em$mu, sd=em$sigma )
}
## PLOT mclust
##
plot.mclust <- function(mc.obj,n,...)
{
  if (class(mc.obj)!="Mclust")
    stop("Mclust object expected: ", class(mc.obj))
  mc.obj <- mc.obj$parameters
  return(plot.mixnorm(n=n,p=mc.obj$pro,mu=mc.obj$mean,sigma=sqrt(mc.obj$variance$sigmasq),...))
}
## PREDICT function
##
predict.mixnorm <- function(obj, x, uniform=FALSE)
{
  if (class(obj)!="Mclust")
    stop("Mclust object expected: ", class(obj))

  params <- obj$parameters

  if ( uniform ) {
    params$pro <- rep(1/obj$G,obj$G)
    cat("pro:"); print(params$pro)
  }
  OUT <- matrix(0.0,length(x),obj$G,dimnames=list(names(x),paste("P(C=",1:obj$G,"|x)",sep="")))
  for ( i in 1:obj$G ) {
    OUT[,i] <- cbind(params$pro[i] * dnorm(x,params$mu[i],params$variance$sigmasq[i]))
  }
  cbind(x=x,OUT,map=t(apply(OUT,1,order,decreasing=T)))
  #cbind(x=x,OUT,map=apply(OUT,1,which.max))
}
