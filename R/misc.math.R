#################################################
## function: LOG PLUS
#################################################
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
#################################################
## function: LOG MINUS
#################################################
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
#################################################
## function: PRECISION
#################################################
# Compute precision of a variable as half the min distance between two
# non-equal variable observations

precision <- function( x, dir=1 )
{
  # input:
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
#################################################
## function: FAST SUM
#################################################
fast.sum <- function(x,cls=NULL,na.rm=FALSE)
{
  if (is.vector(x))
    return( fast.sum(matrix(x,1,length(x)),cls=cls) )
  if (is.null(dim(x)))
    stop( "a 2D object expected" )
  if ( any(is.na(x)) && !na.rm )
    stop( "NA values not admissible unless na.rm==T" )

  In <- rep(1,ncol(x))
  x[is.na(x)] <- 0
  return (drop(x %*% In))
}
#################################################
## function: FAST MEAN
#################################################
#' @export
fast.mean <- function( x, cls=NULL, na.rm=FALSE )
{
  if (is.vector(x))
    return( fast.mean(matrix(x,1,length(x)),cls=cls) )
  if (is.null(dim(x)))
    stop( "a 2D object expected" )
  if ( any(is.na(x)) && !na.rm )
    stop( "NA values not admissible unless na.rm==TRUE" )

  In <- rep(1,ncol(x))

  if ( is.null(cls) )
  {
    nc <- (!is.na(x)) %*% In
    x[is.na(x)] <- 0
    return ( (x %*% In)/ncol(x) )
  }
  # ELSE ..
  #
  lev <- sort(unique(cls))
  In <- sapply( 1:length(lev), function (z) as.numeric(cls==lev[z]) )
  #nc  <- apply(In,2,sum)
  #x.mn <- t( t(x %*% In)/nc )

  nc <- (!is.na(x)) %*% In
  x[is.na(x)] <- 0
  x.mn <- (x %*% In)/nc
  colnames(x.mn) <- if (is.null(levels(cls))) lev else levels(cls)

  return(x.mn)
}
#################################################
## function: FAST SD
#################################################
#' @export
fast.sd <- function( x, cls=NULL, do.sqrt=TRUE, pooled=FALSE, var.equal=TRUE, na.rm=FALSE )
{
  if (is.null(dim(x))) stop( "a 2D object expected" )
  if (pooled && is.null(cls)) stop( "cls needed when asking for pooled sd" )

  In <- rep(1,ncol(x))
  nc <- (!is.na(x)) %*% In
  if ( any(nc<2) )
    stop( "can't compute stdv w/ less than two observations" )

  if ( !is.null(cls) ) {
    lev <- sort(unique(cls))
    In <- sapply( 1:length(lev), function (z) as.numeric(cls==lev[z]) )
    nc <- (!is.na(x)) %*% In
  }
  x[is.na(x)] <- 0
  s  <- x %*% In
  s2 <- x^2 %*% In
  s2 <- (s2 - s^2/nc) / (nc-1)
  s2[s2<0] <- 0

  if (pooled)
  {
    Im <- rep(1,ncol(s2))
    s2 <- if (var.equal) {
      ( (((nc-1)*s2) %*% Im) / (length(cls)-ncol(s2)) ) * ((1/nc) %*% Im)
    }
    else
      t(t(s2)/nc) %*% rep(1,ncol(s2))
  }
  if (do.sqrt)
    s2 <- sqrt( s2 )

  if (!is.null(cls))
    colnames(s2) <- if (is.null(levels(cls))) lev else levels(cls)

  return(s2)
}
#################################################
## function: PART COR
#################################################
part.cor <- function(X,Y,Z)
{
    ## compute partial correlation between X and Y, conditioning on Z
    ##
    X1 <- lm(X~Z)$residuals
    Y1 <- lm(Y~Z)$residuals
    cor(X1,Y1)
}
#################################################
## function: CLS MEDIAN
#################################################
cls.median <- function(DAT,cls)
{
  DATctr <- sapply(unique(cls),function(z) apply(DAT[,cls==z],1,median))
  colnames(DATctr) <- levels(cls)
  DATctr
}
