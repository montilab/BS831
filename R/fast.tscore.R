## function: FAST TSCORE
##
## fast t-score calculation based on matrix operations
## (see vignettes/docs/matrixOp.Rmd for a rationale)
##
fast.tscore <- function(x, cls=NULL, y=NULL, generalized=FALSE, do.test=FALSE,
                        var.equal=FALSE, min.sd=NULL, alternative=c("two.sided","greater","less"),
                        verbose=TRUE )
{
  ## INPUT:
  ##    x - m x n1 matrix (genes by experiments, condition 1)
  ##    y - m x n2 matrix (genes by experiments, condition 2)
  ##  OR
  ##    x - m x n  matrix (genes by experiments, condition 1 & 2)
  ##  cls - n vector of class labels
  ##
  ## OUTPUT:
  ##  score - m vector (positive are upregulated for x, or for
  ##          lower label -- condition 1 -- when cls is specified)

  ## BEGIN input check
  ##
  alternative <- match.arg(alternative)

  if ( !xor(is.null(y),is.null(cls)) )
    stop( "must specify either y or cls" )

  if ( is.null(y) )
  {
    lev <- sort(unique(cls))
    if ( ncol(x)!=length(cls) )
      stop( "ncol(x) must be same as length(cls)" )
    if ( length(lev)>2 )
      stop( "cls must be binary" )
    y <- x[,cls==lev[2],drop=FALSE]
    x <- x[,cls==lev[1],drop=FALSE]
  }
  if ( nrow(x)!=nrow(y) ) stop( "x and y must be of same length\n" )
  if ( ncol(x)<3 ) stop( "x has less than 3 observations\n" )
  if ( ncol(y)<3 ) stop( "y has less than 3 observations\n" )
  ##
  ## END input check

  score <- NULL
  x.idx <- 1:ncol(x)

  VERBOSE( verbose, ifelse( robust, "\tWilcoxon test .. ", "\tt test .. " ) )

  n1 <- (ncol(x)); if (n1<2) stop( "need at least 2 obs per class" )
  n2 <- (ncol(y)); if (n2<2) stop( "need at least 2 obs per class" )
  cls <- c( rep(1,n1), rep(0,n2) ); cls <- cbind( cls, 1-cls )
  x <- cbind(x,y)

  ## FAST COMPUTATION BASED ON MATRIX MULTIPLICATION
  ##
  s  <- x %*% cls    # sum
  s2 <- x^2 %*% cls  # sum of squares
  s2[,1] <- (s2[,1] - (s[,1]^2)/n1) / (n1-1) # variance in 1st class
  s2[,2] <- (s2[,2] - (s[,2]^2)/n2) / (n2-1) # variance in 2nd class
  s[,1] <- s[,1]/n1  # mean in 1st class
  s[,2] <- s[,2]/n2  # mean in 2nd class

  ## variance thresholding (if provided)
  if ( !is.null(min.sd) ) {
    min.var <- min.sd*min.sd
    s2[s2[,1]<min.var,1] <- min.var
    s2[s2[,2]<min.var,2] <- min.var
  }
  stderr <- if (var.equal)
    sqrt( (((n1-1)*s2[,1] + (n2-1)*s2[,2])/(n1+n2-2)) * (1/n1+1/n2) )
  else
    sqrt( s2[,1]/n1 + s2[,2]/n2 )

  score <- (s[,1]-s[,2]) / stderr

  if ( do.test )
  {
    df <- if ( var.equal ) # degrees of freedom
      n1+n2-2
    else
      stderr^4 / ( (s2[,1]/n1)^2/(n1-1) + (s2[,2]/n2)^2/(n2-1)) # Welch approximation of df

    pval <- if (alternative == "less") {
      pt(score, df=df)
    }
    else if (alternative == "greater") {
      pt(score, df=df, lower.tail=FALSE)
    }
    else {
      2 * pt(-abs(score), df=df)
    }
    score <- cbind( score=score, p.value=pval )
  }
  VERBOSE( verbose, "done.\n" )
  return( score )
}
