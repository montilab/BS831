#########################################################################
##                       SUPPORT FUNCTIONS
#########################################################################

my.nlevels <- function( x )
{
  return( length(my.levels(x)) )
}
my.levels <- function( x, sort=TRUE )
{
  return( if (sort) sort(unique(as.vector(x))) else unique(as.vector(x)) )
}
# function: PRECISION
#
precision <- function( x, dir=1 )
{
  # compute a variable's precision as half the min distance btw two
  # non-equal observations
  #
  # INPUT:
  #    x - one or two dimensional array
  #  dir - compute precision by row (dir=1) or by column (dir=2)

  if ( dir!=1 && dir!=2 ) stop( "dir must be either 1 or 2: ", dir );
  
  # x can be a list/array ..
  #
  if ( is.null(dim(x)) )
  {
    x <- sort(x)
    x <- x[2:length(x)]-x[1:(length(x)-1)]
    x[x<=0] <- +Inf
    return ( min(x)/2 )   
  }
  # ..or a two-dimensional structure
  #
  else {
    return( apply(x, dir, precision) )
  }
}
# function: LOG PLUS
#
log.plus <- function( x, y=NULL )
{
  # Compute log(exp(x) + exp(y)) without ever "leaving" the log-base
  # (useful when summing over very small numbers). Either an array
  # (y==NULL) or two arguments (either scalar or vectors) are provided
  # as input
  #
  if ( is.null(y) )
  {
    if (length(x)<2) stop( "either an array or y expected" )
    i.max <- which.max(x)
    x.max <- x[i.max]; x <- x[-i.max]
    x.max + log( 1 + sum(exp(x-x.max)) )    
  }
  else if ( is.vector(y) )
  {
    log.plus(as.matrix(x),as.matrix(y))
  }
  else if ( is.matrix(x) )
  {
    if ( !is.matrix(y) ) stop( "y must be a matrix when x is a matrix" )
    if ( any(dim(x)!=dim(y)) ) stop( "x and y must have same dim" )

    Xmax <- pmax(x,y)
    Xmax + log(1 + exp(pmin(x,y)-Xmax))
  }
  else {
    log.plus( c(x,y) )
  }
}
log.minus <- function( x, y )
{
  # Compute log(exp(x) - exp(y)) without ever "leaving" the log-base
  # (useful when subtracting over very small numbers).
  #
  if (length(x)!=length(y))
    stop( "x and y must be of same length" )
  xy <- rep( -Inf, length(x) )
  idx <- x>y
  xy[idx] <- x + log( 1 - exp(y-x) )
  return( xy )
}
#########################################################################
##                          MAIN FUNCTIONS
#########################################################################


##################################################################
## class definition: NAIVE BAYES
##################################################################
#
setClass("naive.bayes",
         representation(p="vector",      # class priors
                        mu="matrix",     # class-conditional predictor means
                        sigma="matrix"), # class-conditional predictor stdev's
         validity = function(object) {
           if ( all(length(object@p)==nrow(object@mu),
                    length(object@p)==nrow(object@sigma),
                    ncol(object@mu)==ncol(object@sigma)) ) TRUE
           else paste( "Incompatible dimensions for p,mu,sigma: ",
                      length(object@p), ", [",
                      paste( dim(object@mu), collapse="," ), "], [",
                      paste( dim(object@sigma), collapse="," ), "]",
                      sep="")
         }
         )
names.naive.bayes <- function(obj) colnames(obj@mu)

##################################################################
## function: NAIVE ESTIMATE
##################################################################
##
naive.estimate <- function( dat, cls, nclass=0, priorp=0, min.sd=-1, uniform=FALSE )
{
  ## Estimate the parameters of a naive bayes model from a set of training cases
  ##
  ## INPUT
  ##       dat: M x N matrix (M rows/cases, N columns/predictors)
  ##       cls: M-vector of class labels
  ##    nclass: number of classes; if 0, count distinct values in cls
  ##    priorp: smoothing factor for class prior probability
  ##     minsd: minimum standard deviation allowed (-1: precision; 0: sample
  ##            variance; else: given value )
  ##
  ## NOTE: the smoother factor priorp is such that the class probability
  ## is computed as
  ##
  ##   P(C=c) = (N_c + priorp) / Sum_c ( N_c + priorp)
  ##
  ## therefore, if priorp=0 we have the maximum likelihood estimates,
  ## if priorp=1 we obtain the "Laplace rule of succession". The larger
  ## priorp, the more the class prior approximates a uniform
  ## distribution.

  if ( nrow(dat)!=length(cls) )
    stop( "number of cases different from number of class labels" )
  if ( is.factor(cls) ) # map class labels to 0:n-1 range
    cls <- match(cls,levels(cls))-1
  if ( any((match(cls,my.levels(cls))-1)!=cls) )
    stop("class labels must be contiguous integers starting at 0")

  M <- dim(dat)[1]                       # no. of cases
  N <- dim(dat)[2]                       # no. of predictors
  if ( nclass==0 ) {                     # no. of classes
    nclass <- length(unique(cls)) 
  }
  div <- tabulate( cls+1, nbins=nclass ) # number of cases per class
  if ( any(div==0) )
    stop("No cases for class(es) ", paste( which(div==0), collapse=", ") )

  p <- rep(1.0/nclass,nclass)
  if ( !uniform ) {
    p <- div+priorp                      # estimate class priors
    p <- p/sum(p)                        # ..
  }
  ## handling of zero variance (selection of default sd)
  ##
  if ( min.sd==-1 ) {     # set to variable precision (if more than one values)
    min.sd <- matrix( rep(precision(dat,2),nclass), nclass, N, byrow=TRUE )
    min.sd[min.sd==Inf] <-
      matrix(rep(apply(dat,2,sd),nclass), nclass, N, byrow=TRUE)[min.sd==Inf]
  }
  else if ( min.sd==0 ) { # set to overall variable variance
    min.sd <- matrix( rep(apply(dat,2,sd),nclass), nclass, N, byrow=TRUE )
  }
  else {                  # set to given variance
    if ( min.sd<0 ) stop( "min.sd must be positive" )
    min.sd <- matrix( min.sd, nclass, N )
  }
  ## estimation of model parameters
  ##
  CLS <- sapply(0:(nclass-1),function(x){as.numeric(cls==x)})
  mu <- ( t(CLS) %*% dat )
  sg <- sqrt( ( (t(CLS) %*% dat^2) - (mu^2)/div ) / (div-1) )
  mu <- mu/div
  rownames(mu) <- rownames(sg) <- paste( "class", 1:nrow(mu), sep=".")

  ## handling of zero variance (setting to default sd)
  ##
  sg[is.nan(sg)] <- min.sd[is.nan(sg)]
  sg[sg==0] <- min.sd[sg==0]
  if ( (ninf <- sum(sg==Inf))>0 ) warning( ninf, " SD's set to Inf" )

  nb <- new( "naive.bayes", p=p, mu=mu, sigma=sg )
}
##################################################################
## function: NAIVE PREDICT (IN)
##################################################################
##
naive.predict <- function( dat, model, logmath=TRUE, do.match=TRUE )
{
  if (class(model)!="naive.bayes") stop( "'naive.bayes' object expected" )

  ## check for matching between 'training' and test predictors
  if (do.match) {
    midx <- match(names(model),colnames(dat))
    if ( any(is.na(midx))) stop( "missing variables" )
    dat <- dat[,midx,drop=FALSE]
  }
  else if (any(colnames(dat)!=names(model)))
    stop( "models' variables do not match dat's columns" )
  
  return( naive.predict.in( dat, model@p, model@mu, model@sigma, logmath=logmath ) )
}
naive.predict.in <- function( x, p, mu, sigma, logmath=TRUE )
{
  ## input
  ##       x: M x N matrix (M rows/cases, N columns/predictors)
  ##       p: K-vector of class priors (K is the number of classes)
  ##       m: K x N matrix of means of conditional Gaussians
  ##       s: K x N matrix of stdev of conditional Gaussians
  ## logmath: carry out all operations in log space to avoid underflow
  ##
  if ( is.null(dim(x)) )
    x <- t(as.matrix(x))
  if ( dim(x)[2]!=dim(mu)[2] ) stop( "x's cols inconsistent with m's cols" )
  if ( dim(x)[2]!=dim(sigma)[2] ) stop( "x's cols inconsistent with s's cols" )
  
  M <- dim(x)[1]
  N <- dim(mu)[2]
  K <- length(p)
  Ik <- rep(1,K)

  joint <- matrix( NA, nrow=M, ncol=K )
  for ( j in (1:K) )
  {
    tmp <- t(apply(x,1,dnorm,mean=mu[j,],sd=sigma[j,],log=TRUE))
    if (N==1) tmp <- t(tmp)
    
    joint[,j] <- if ( logmath ) 
      log(p[j]) + apply( tmp,1,sum )
    else
      p[j] * apply( tmp,1,prod )
  }
  pred <- if ( logmath ) {
    exp( joint - apply(joint,1,log.plus) )
  }
  else {
    joint / drop(joint %*% Ik)
  }
  colnames(pred) <- paste( "P(C=", 1:length(p), ")", sep="" )
  pred
}
##################################################################
## function: NAIVE EVAL
##################################################################

naive.eval <- function( dat, cls, model, logmath=TRUE, rnd=3)
{
  nb.out <- naive.predict( dat, model=model, logmath=logmath )
  pred <- apply(nb.out,1,which.max)-1
  nb.out <- cbind(round(nb.out,rnd),
                  pred=pred,
                  true=cls,
                  error=as.numeric(pred!=cls))
  nb.out
}
