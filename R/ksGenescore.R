## FUNCTION: KS GENESCORE
##
## simple wrapper of function ks.test to allow for plot generation and signed statistic
##
#' @export
ksGenescore <- function(
 n.x,               # length of ranked list
 y,                 # positions of geneset items in ranked list (basically, ranks)
 do.pval=TRUE,      # compute asymptotic p-value
 alternative=c("two.sided","greater","less"),
 do.plot=F,         # draw the ES plot
 bare=FALSE,        # return score & p-value only (a 2-tuple)
 cls.lev=c(0,1),    # class labels to display
 absolute=FALSE,    # takes max - min score rather than the maximum deviation from null
 plot.labels=FALSE, # hits' labels
 ...                # additional plot arguments
 )
{
  # efficient version of ks.score (should give same results as ks.test, when weight=NULL)
  #
  alternative <- match.arg(alternative)
  #DNAME <- paste( "1:", n.x, " and ", deparse(substitute(y)), sep="" )
  #METHOD <- "Two-sample Kolmogorov-Smirnov test"
  n.y <- length(y)
  if ( n.y < 1 )  stop("Not enough y data")
  if ( any(y>n.x) ) stop( "y must be <= n.x: ", max(y) )
  if ( any(y<1) ) stop( "y must be positive: ", min(y) )
  x.axis <- y.axis <- NULL

  # KS score
  #
  y <- sort(y)
  n <- n.x * n.y/(n.x + n.y)
  hit <- 1/n.y
  mis <- 1/n.x

  ## to compute score, only the y positions and their immediate preceding
  ## ..positions are needed
  ##
  Y <- sort(c(y-1,y)); Y <- Y[diff(Y)!=0]; y.match <- match(y,Y)
  D <- rep( 0, length(Y) ); D[y.match] <- (1:n.y)
  zero <- which(D==0)[-1]; D[zero] <- D[zero-1]
  
  z <- D*hit - Y*mis
  
  score <- if (absolute) max(z)-min(z) else z[which.max(abs(z))]
  names(score) <- "D"

  if (do.plot) {
      x.axis <- Y;
      y.axis <- z;
      if(Y[1]>0) {
          x.axis <- c(0,x.axis);
          y.axis <- c(0,y.axis);
      }
      if ( max(Y)<n.x ) {
          x.axis <- c(x.axis,n.x)
          y.axis <- c(y.axis,0)
      }
      plot( x.axis, y.axis, type="l",
           xlab=paste("up-regulated for class ", cls.lev[2], " (KS>0) vs ",
               "up-regulated for class ", cls.lev[1], " (KS<0)", sep="" ),
           ylab="gene hits",...)
      abline(h=0)
      abline(v=n.x/2,lty=3)
      axis(1,at=y,labels=plot.labels,tcl=0.25,las=2)
      i.max <- which.max(abs(y.axis))
      points( x.axis[i.max], y.axis[i.max], pch=20, col="red")
      text(x.axis[i.max]+n.x/20,y.axis[i.max],round(y.axis[i.max],2))
  }
  if ( !do.pval ) {
      return(score)
  }
  ## ELSE compute p-value as in function ks.test but return signed statistic
  ##
  tmp <- suppressWarnings(ks.test(1:n.x,y,alternative=alternative))
  tmp$statistic <- score # use the signed statistic
  return( if (bare) c(tmp$statistic, tmp$p.value) else tmp )
}  
