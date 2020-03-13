#' subset an expression dataset by extracting genes based on their level of variation 
#'
#' @param dat an expressionSet object
#' @param score type of variation measure to use (can be 'mad', 'sd', or 'cv')
#' @param dir direction of the gene ranking ('top': highest varying genes; 'bottom': lowest varying genes)
#' @param transform transform the data before measuring variation
#' @param ngenes number of top (or bottom) genes to extract
#' @param do.plot plot the center vs. scale plot with selected genes highlighted
#' @param do.log indicates whether to log the "x", "y", or "xy" axes, or "" none
#'
#' @return the filtered expressionSet object
#' 
#' @examples
#'
#' # select the top 100 genes by median absolute deviation (mad) and also generate plot
#'
#' # Use example data, for data set information: ?eSet.brca.100
#' data(eSet.brca.100)
#'
#' variationFilter(eSet.brca.100,score='mad',ngenes=50,do.plot=TRUE)
#' 
#' @export
VERBOSE <- function( v, ... )
{
  if ( v ) cat( ... )
}
#' @export 

########################################################################
#                              MAIN                                    #
########################################################################
variationFilter <- function(dat,
                            score=c("mad","sd","cv"),
                            dir=c("top","bottom"),
                            transform=c("none","log2","exp2","log","exp"),
                            ngenes=NULL,
                            min.score=NULL,
                            min.log=1,
                            rnd=4,
                            do.plot=FALSE,
                            pch=".",
                            lgnd.coord=1,
                            do.log="",
                            qnt.lev=0.5,
                            min.qnt=-Inf,
                            no.affx=FALSE,
                            verbose=TRUE)
{
  if (is.null(ngenes) && is.null(min.score) )
    stop( "must specify either ngenes or min.score" )
  if (!is.null(ngenes) && !is.null(min.score) )
    stop( "cannot specify both ngenes and min.score" )
  if (!is.null(ngenes) && ngenes>nrow(dat) )
    stop( "ngenes is too large" )
  if (min.log<=0)
    stop( "min.log must be positive: ", min.log )
  if ( class(dat)!='ExpressionSet' )
    stop( "ExpressionSet object expected: ",class(dat) )

  transform <- match.arg(transform)
  dir <- match.arg(dir)
  score <- match.arg(score)
  score.fun <- match.fun(score)

  if (transform=="log2" || transform=="log") { # threshold before log-transformation
    VERBOSE( verbose, "Thresholding before log-transforming .. " )
    Biobase::exprs(dat)[Biobase::exprs(dat)<min.log] <- min.log
    VERBOSE( verbose, "done.\n" )
  }
  Biobase::exprs(dat) <- switch(transform,
                       none=Biobase::exprs(dat),
                       log2=round(log2(Biobase::exprs(dat)),rnd),
                       exp2=round(2^(Biobase::exprs(dat)),rnd),
                       log=round(log(Biobase::exprs(dat)),rnd),
                       exp=round(exp(Biobase::exprs(dat)),rnd))

  if (no.affx) {
    if ( length(rm.idx <- grep("AFFX-",featureNames(dat)))>0 ) {
      VERBOSE(verbose,"Removing 'AFFX-' probes ..")
      dat <- dat[-rm.idx,,drop=FALSE]
      VERBOSE(verbose," done,", length(rm.idx),"removed.\n")
    }
  }
  ctr <- if (score=="mad") apply( Biobase::exprs(dat), 1, median ) else rowMeans( Biobase::exprs(dat) )
  SC <- SC1 <- apply( Biobase::exprs(dat), 1, score.fun )

  if (min.qnt>0)
  {
    VERBOSE(verbose, "Filtering out genes w/ ",round(100*qnt.lev,2), "-percentile < ", min.qnt, " .. ",sep="" )
    QNT <- apply(Biobase::exprs(dat),1,quantile,probs=qnt.lev)
    if ( sum(QNT>=min.qnt)<2 )
      stop( "filtering by min.qnt returns less than 2 genes (try decreasing min.qnt)" )
    dat <- dat[QNT>=min.qnt,,drop=FALSE]
    VERBOSE(verbose, "done,", nrow(dat), "genes left.\n")

    if ( !is.null(ngenes) && nrow(dat)<=ngenes ) {
      VERBOSE(verbose,"Number of genes left is less than required, no further filtering necessary")
      return(dat)
    }
    SC1 <- SC[QNT>=min.qnt]
  }
  VERBOSE(verbose, "Variation filtering based on", score, ".. " )
  
  VERBOSE(verbose, "done.\n" )
  
  idx <- NULL
  if ( is.null(ngenes) ) {
    VERBOSE(verbose, "Selecting genes with", score, "<=", min.score, ".. " )
    idx <- if(dir=="top") SC1>=min.score else SC1<=min.score
    if (sum(idx)==0)
      stop( "no genes passed the filtering criteria" )
  }
  else {
    VERBOSE(verbose, "Selecting top", ngenes, "by", score, ".. " )
    if (dir=="top") SC1 <- -SC1
    idx <- order(SC1)[1:ngenes]
  }
  dat <- dat[idx,,drop=FALSE]
  VERBOSE(verbose, "done,", nrow(dat), "genes selected.\n" )

  if (do.plot) {
    VERBOSE(verbose, "Creating scatter plot .. ")
    if (is.null(do.log)) {
      do.log <- if (transform=="none" || transform=="exp2" || transform=="exp" )
        "xy"
      else
        ""
    }
    SC <- abs(SC)
    plot( ctr, SC, pch=pch, col="gray", xlab=if (score=="mad") "median" else "mean", ylab=score, log=do.log)
    plot.idx <- match( featureNames(dat),names(SC) )
    points( ctr[plot.idx], SC[plot.idx],pch=pch,col="red")
    lx <- min(ctr); ly <- max(SC); xjust <- 0; yjust <- 1
    if (lgnd.coord==2) {
      lx <- max(ctr); xjust <- 1
    }
    else if (lgnd.coord==3) {
      lx <- max(ctr); ly <- min(SC); xjust <- 1; yjust <- 0
    }
    else if (lgnd.coord==4) {
      ly <- min(SC); yjust <- 0
    }
    else if (lgnd.coord!=1)
      stop( "lgnd.coord must be btw 1 and 4" )
    
    legend(lx, ly, xjust=xjust, yjust=yjust,
           legend=c("all","passing filter"), col=c("gray","red"), pch=20)
    VERBOSE(verbose, "done.\n")
  }
  dat
}
