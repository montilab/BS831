#' @export
VERBOSE <- function( v, ... )
{
  if ( v ) cat( ... )
}
assert <- function( expr, msg=NULL )
{
  if ( !eval(parse(text=expr)) )
    stop( if (is.null(msg)) paste("'",expr,"' failed",sep="") else msg )
}
my.nlevels <- function( x )
{
  return( length(my.levels(x)) )
}
my.levels <- function( x, sort=T )
{
  #if ( !is.numeric(x) ) stop( "x must be numeric" )
  return( if (sort) sort(unique(as.vector(x)),na.last=T) else unique(as.vector(x)) )
}
my.tabulate <- function( x,nbins=length(unique(x)) )
{
  tabulate(match(x,unique(x)),nbins=nbins)
}
rowMax <- function( MX )
{
  if ( !(class(MX) %in% c('matrix','data.frame')) )
    stop( "class(MX): ", class(MX) )
  apply(MX,1,max)
}
rowMin <- function( MX )
{
  if ( !(class(MX) %in% c('matrix','data.frame')) )
    stop( "class(MX): ", class(MX) )
  apply(MX,1,min)
}
which.rowMax <- function( MX )
{
  if ( !(class(MX) %in% c('matrix','data.frame')) )
    stop( "class(MX): ", class(MX) )
  apply(MX,1,which.max)
}
which.rowMin <- function( MX )
{
  if ( !(class(MX) %in% c('matrix','data.frame')) )
    stop( "class(MX): ", class(MX) )
  apply(MX,1,which.min)
}
cv <- function(x,mad=F)
{
  if ( is.matrix(x) )
    apply( x, 2, cv, mad=mad )
  else if ( is.vector(x) ) {
    m <- abs( if ( mad ) median(x) else mean(x) )
    if ( m==0 ) stop( "mean(x) is 0" )
    s <- if ( mad ) mad(x) else sd(x)
    s/m
  }
  else if (is.data.frame(x))
    sapply(x, cv, mad=mad )
  else cv(as.vector(x),mad=mad)
}
cvm <- function(x,mad=F)
{
  if ( is.matrix(x) ) {
    apply( x, 2, cvm, mad=mad )
  }
  else if ( is.vector(x) )
  {
    if ( mad ) {
      m <- abs(median(x))
      median( abs(x-m)/m )
    }
    else {
      m <- abs( mean(x) )
      if ( m==0 ) stop( "mean(x) is 0" )
      s <- sd(x)
      s/m
    }
  }
  else if (is.data.frame(x)) {
    sapply(x, cvm, mad=mad )
  }
  else {
    cvm(as.vector(x),mad=mad)
  }
}
lookup.fun <- function( x, X, Y )
{
  i <- sum(X<=x)+1
  #if ( i>length(Y) ) stop("index out of bound")
  return( Y[i] )
}
var2string <- function( x )
{
  deparse(substitute(x))
}
strip.path <- function( filen, path=F )
{
  splits <- unlist(strsplit(filen,"\\/"))
  if ( path )
    paste(splits[-length(splits)],collapse="/")
  else
    rev(splits)[1]
}
make.dir <- function( dname )
{
  if ( is.na(file.info(dname)$isdir) ) system(paste('mkdir',dname))
}
#' @export
file.ext <- function( filen, w.sep=F )
{
  ext <- rev(unlist(strsplit( filen, "\\.")))[1]
  if ( w.sep )
    ext <- paste(".",ext,sep="")
  ext
}
file.stub <- function( filen, w.sep=F )
{
  stub <- paste(rev(rev(unlist(strsplit( filen, "\\.")))[-1]),collapse=".")
  if ( w.sep )
    stub <- paste(stub,".",sep="")
  stub
}
file.test <- function( filen, mode=4, string=NULL, nostop=F )
{
  action <- c("exist","execute","write","","read")[mode+1]
  if ( file.access(filen,mode)[1]!=0 ) {
    out.string <- if (is.null(string)) paste("cannot ",action," file '",filen,"'",sep="") else string
    if (nostop)
      warning( out.string )
    else
      stop( out.string )
    return( F )
  }
  else
    return( T )
}
dir.test <- function( dirname )
{
  dirinfo <- file.info(dirname)
  if ( is.na(dirinfo$size) ) {
    warning( "directory '", dirname, "' doesn't exist", sep="" )
    return( F )
  }
  if ( !dirinfo$isdir ) {
    warning( "'", dirname, "' is not a directory", sep="" )
    return( F )
  }
  T
}
#' @export
mmatch <- function( x, y )
{
  #unlist(lapply( x, function(z) which(y==z) ))
  unlist(lapply( x, function(z) { tmp <- which(y==z); names(tmp) <- rep(z,length(tmp)); tmp}))
}
#' @export
match.nona <- function( a1, a2, na.rm=FALSE, suppressWarnings=FALSE )
{
  if ( any(is.na(idx <- match(a1,a2))) )
    if ( na.rm ) {
      if ( !suppressWarnings )
        warning("mismatch | missing in a2 (",sum(is.na(idx)),"): ",
                paste(a1[is.na(idx)],collapse=", ") )
    }
    else
      stop( "mismatch | missing in a2: ", paste(a1[is.na(idx)],collapse=", ") )
  idx[!is.na(idx)]
}
#' @export
matchIndex <- function( key, names, ignore.case=FALSE )
{
  if ( ignore.case ) {
    key <- toupper(key); names <- toupper(names)
  }
  if ( sum(key==names,na.rm = TRUE)>1 ) {
    warning('multiple matches, returning first only:',key)
  }
  idx <- match(key,names)
  if (is.na(idx)) stop( "index not found: ", key )
  return( idx )
}
#' @export
plot.norm <- function( n=1000, mean=0, sd=1, add=F, lty="solid", ...)
{
  p <- (1:(n-1))/n
  if ( add ) {
    lines( qnorm(p,mean=mean,sd=sd), dnorm( qnorm(p,mean=mean,sd=sd),mean=mean,sd=sd ),type="l", lty=lty,...)
  }
  else {
    plot( qnorm(p,mean=mean,sd=sd), dnorm( qnorm(p,mean=mean,sd=sd),mean=mean,sd=sd ),type="l", lty=lty, ...)
  }
}
#' @export
plot.distn <- function( qfun, dfun, n, lty=1, add=FALSE, xlab=NULL, ylab=NULL, main=NULL, ... )
{
  p <- (1:(n-1))/n
  q <- qfun(p,...)

  if ( add )
    lines( q, dfun(q,...), lty=lty )
  else
    plot( q, dfun(q,...), type="l", lty=lty, xlab=xlab,ylab=ylab,main=main )
}
#' @export
cumineq <- function( prm, obs, dir=1, debug=F )
{
  ## INPUT:
  ##  - prm    n-sized array
  ##  - obs    n-sized array
  ## WHAT:
  ##  (support function for pval2fdr)
  ##  for each entry in obs, count how many entries in prm
  ##  are <= (dir=1) or >= (dir=2) than that entry
  ##
  p.ord <- order(if ( dir==1 ) prm else -prm)
  o.ord <- order(if ( dir==1 ) obs else -obs)
  o.rnk <- rank(if ( dir==1 ) obs else -obs)

  ## sort entries
  ##
  prm <- prm[p.ord]
  obs <- obs[o.ord]

  u.obs <- unique(obs)
  cup <- c(prm,u.obs)
  cup <- cup[order(if (dir==1) cup else -cup)]
  fp <- length(cup)+1-match(obs,rev(cup))-match(obs,u.obs)

  ## return values sorted according to original order (o.rnk)
  ##
  return ( if (debug)
             cbind( prm[o.rnk], obs[o.rnk], fp[o.rnk] )
           else
             fp[o.rnk] )
}
#' @export
pval2fdr <- function( p, monotone=T, nh=length(p), na.rm=F )
{
  if (length(p)==1)
    return(p)

  if (nh<length(p)) warning( "nh<length(p)" )

  na.idx <- is.na(p)
  if ( !na.rm && any(na.idx) ) {
    stop( "use na.rm=T to handle NA's")
  }
  if ( any(na.idx) ) {
    p <- p[!na.idx]
    nh <- nh-sum(na.idx)
  }
  p.ord <- order(p)
  p1 <- p[p.ord]
  fdr <- p1 * nh / cumineq(p1,p1)
  if (monotone) for ( i in (length(p)-1):1 ) fdr[i] <- min(fdr[i],fdr[i+1])
  fdr[fdr>1] <- 1
  fdr <- fdr[rank(p)]
  names(fdr) <- names(p)

  if ( any(na.idx) ) {
    tmp <- rep(NA,length(na.idx))
    tmp[!na.idx] <- fdr
    fdr <- tmp
  }
  fdr
}
lmchoose <- function( x )
{
  if (any(x<0)) stop( "all terms must be non-negative" )
  if (any(x-round(x)!=0)) stop( "all terms must be integer" )
  lfactorial(sum(x))-sum(lfactorial(x))
}
mchoose <- function( x )
{
  exp(lmchoose(x))
}
list.intersect <- function( X )
{
  if (!is.list(X)) stop( "list expected: ", class(X) )
  if (length(X)<=1) return(X[[1]])
  INT <- X[[1]]
  for ( i in 2:length(X) ) {
    INT <- intersect(X[[i]])
  }
  INT
}
upper.case <- function( X )
{
#  sapply(X,function(z) gsub("(.*)","\\U\\1",z,perl=T))
  toupper(X)
}
#' @export
## LOAD VAR
##
load.var <- function(file,verbose=F)
{
  ## load a variable stored in an R-binary file

  varlist <- c(ls(),"varlist")
  load(file=file)
  if ( length(setDiff <- setdiff(ls(),varlist))>1 ) stop( "loading more than one variable" )
  VERBOSE(verbose,"(variable loaded: ",setdiff(ls(),varlist),")\n",sep="")
  return( eval(parse(text=setDiff)) )
}
robust.load <- function( file, envir=parent.env(environment()), max.try=100, verbose=T )
{
  ntry <- 0
  if ( file.access(file)!=0 ) {
    stop("cannot read '", file, "'")
  }
  while( ntry<max.try ) {
    out <- try(load(file=file,envir=envir),silent=T)
    if (class(out)!="try-error")
      return(T)
    ntry <- ntry+1
    Sys.sleep(5)
  }
  warning("robust load attempts exhausted, trying one last time (likely to fail)")
  load(file=file,envir=envir)
}
rjust <- function(x,n)
{
  format(x,width=n,justify="right")
}
ljust <- function(x,n)
{
  format(x,width=n,justify="left")
}
segment.overlap <- function(seg1,seg2,do.plot=F)
{
  ## is 1st segment within, left-overlap, right-overlap, containing 2nd segment?
  ##
  ovlp.inn <- seg1[1]>=seg2[1] && seg1[2]<=seg2[2]
  ovlp.lft <- seg1[1]<=seg2[1] && seg1[2]>=seg2[1] && seg1[2]<=seg2[2]
  ovlp.rgt <- seg1[1]>=seg2[1] && seg1[1]<=seg2[2] && seg1[2]>=seg2[2]
  ovlp.out <- seg1[1]<=seg2[1] && seg1[2]>=seg2[2]

  if ( do.plot ) {
    xlim <- c(min(seg1,seg2),max(seg1,seg2))
    plot( seg1, c(1,1), pch="|", type="b",col="red",xlim=xlim)
    points( seg2, c(1.1,1.1), pch="|", type="b",col="blue")
    legend("topleft",c("segment 1","segment 2"),col=c("red","blue"),pch="|",lty=1)
  }
  return( c(inn=ovlp.inn,lft=ovlp.lft,rgt=ovlp.rgt,out=ovlp.out) )
}
## REPMAT
##
## replicate a matrix or data.frame row- or column-wise

repmat <- function(MX,MARGIN=2,n,each=FALSE)
{
  m <- if ( MARGIN==1 ) ncol(MX) else nrow(MX)
  idx <- if ( each ) rep(1:m,each=n) else rep(1:m,times=n)

  if ( MARGIN==1 )
    return(MX[,idx])
  else if ( MARGIN==2 )
    return(MX[idx,])
  else
    stop("MARGIN:",MARGIN)
}
which.2D <- function( x )
{
  xr <- apply( x, 1, which )
  xr <- unlist(sapply(1:nrow(x),
                      function(z){ if ((n <- length(xr[[z]]))>0) t(cbind(rep(z,n),xr[[z]]))}))
  matrix(xr,length(xr)/2,2,byrow=T)
}
read.csv.delim <- function(file,header=T,sep=",",stringsAsFactors=F,check.names=F,blank.lines.skip=T,...)
{
  read.csv(file,header=header,sep=sep,stringsAsFactors=stringsAsFactors,check.names=check.names,
           blank.lines.skip=blank.lines.skip,...)
}
#' @export
read.tab.delim <- function(file,header=T,sep="\t",stringsAsFactors=F,check.names=F,blank.lines.skip=T,...)
{
  read.delim(file,header=header,sep=sep,stringsAsFactors=stringsAsFactors,check.names=check.names,
             blank.lines.skip=blank.lines.skip,...)
}
#' @export
read.delim.wsave <- function(file,do.save=T,force.read=F,verbose=T,ext=".RData",...)
{
  binfile <- paste(file.stub(file),ext,sep="")
  binfo <- file.info(binfile)
  finfo <- file.info(file)

  #if ( !force.read && file.access(binfile,mode=4)[1]==0 )
  if ( !force.read && !is.na(binfo$size) && binfo$mtime>finfo$mtime )
  {
    VERBOSE(verbose, "Loading binary version '", binfile, "' ..",sep="")
    dat <- load.var(binfile)
    VERBOSE(verbose, "done.\n")
    return(dat)
  }
  ## ELSE ..
  ##
  VERBOSE(verbose,"Reading '",file,"' ..", sep="")
  dat <- if (file.ext(file)=="csv") read.csv.delim(file,...) else read.tab.delim(file,...)
  VERBOSE(verbose,"done.\n")

  if ( do.save ) {
    VERBOSE(verbose,"Saving binary version to '", binfile, "' ..", sep="")
    save(dat,file=binfile)
    VERBOSE(verbose, "done.\n")
  }
  dat
}
my.write.matrix <- function ( x, file = "", sep = "\t",
                              col.names=T,
                              append=F,
                              row.names=F,
                              justify = c( "none", "left", "right"),
                              pval=NULL,
                              newline="\n",
                              names=NULL )
{
  justify = match.arg( justify )
  x <- as.matrix(x)
  p <- ncol(x)
  cnames <- colnames(x)
  rnames <- rownames(x)

  if ( !is.null(pval) ) {
    x[,pval] <- format.pval( as.numeric(x[,pval]) )
  }
  if ( col.names && !is.null(cnames) )
  {
    x <- format(rbind( cnames, x ), justify=justify)
  }
  if ( row.names )
  {
    p <- p+1
    if ( col.names && !is.null(cnames) ) {
      rnames <- if (is.null(names)) c("",rnames) else c(names, rnames)
    }
    x <- cbind( format(rnames,justify=justify), x )
  }
  cat( t(x), file=file, sep=c(rep(sep,p - 1), newline), append=append )
}
matrix2list <- function( MX, MARGIN=2 )
{
  if (MARGIN!=1 && MARGIN!=2) {
    stop( "MARGIN must be 1 or 2: ", MARGIN )
  }
  LS <- vector("list",length=if(MARGIN==2) ncol(MX) else nrow(MX))

  if ( MARGIN==2 ) {
    names(LS) <- colnames(MX)
    for ( i in 1:ncol(MX) ) LS[[i]] <- MX[,i]
  }
  else {
    names(LS) <- rownames(MX)
    for ( i in 1:nrow(MX) ) LS[[i]] <- MX[i,]
  }
  LS
}
interleave <- function(n,m)
{
  ## generate indeces to interleave m vectors of size n
  ## before interleaving:
  ##    1,2,3,..n,n+1,n+2,2n,..,(m-1)*n+1,..,m*n
  ## after interleaving:
  ##    1,n+1,2n+1,..,(m-1)n+1, 2,n+2,..,(m-1)n+2,..,n, n+n,2n+n,..,(m-1)n+n
  ##
  inc <- rep(n*(0:(m-1)),n)
  rep(1:n,each=m)+inc
}
## UNIQUEFY
##
## make all names unique by appending a count

uniquefy <- function( NAMES,sep="." )
{
  if ( length(NAMES)>length(unique(NAMES)) )
  {
    rpl.idx <- unique(NAMES)[which(tabulate(match(NAMES,unique(NAMES)))>1)]
    for ( rpl in rpl.idx ) {
      idx <- NAMES==rpl
      NAMES[idx] <- paste(NAMES[idx],1:sum(idx),sep=sep)
    }
  }
  NAMES
}
## N TABULATE
##
## count the number of occurrences of each unique entry

ntabulate <- function( X, do.tab=F )
{
  tmp <- tabulate(match(X,unique(X)))
  names(tmp) <- unique(X)

  if ( do.tab ) {
    tmp <- data.frame(ID=names(tmp),count=tmp)
    rownames(tmp) <- 1:nrow(tmp)
  }
  tmp
}
list2table <- function( obj, fill=NA )
{
  if ( !is.list(obj) )
    stop( "input object must be a list" )

  mx <- matrix( fill, nrow=max(sapply(obj,length)), ncol=length(obj), dimnames=list(NULL,names(obj)) )
  for ( i in 1:length(obj) ) {
    if ( length(obj[[i]])<1 ) next
    mx[1:length(obj[[i]]),i] <- obj[[i]]
  }
  mx
}
table2list <- function( tab, header=TRUE, fill=NA )
{
  LS <- vector('list',ncol(tab))
  if ( header ) names(LS) <- colnames(tab)
  for ( i in 1:ncol(tab) ) LS[[i]] <- tab[tab[,i]!=fill,i]
  LS
}
## write.table wrapper that adds a column header to row names column, if specified
##
#' @export
my.write.table <- function(x, file="", append=F, sep="\t", row.names=T, col.names=T, names=NULL, dec=".",
                           quote=F, eol="\n", na="NA", qmethod=c("escape", "double"), fileEncoding="",newline=0)
{
  qmethod <- match.arg(qmethod)

  if ( !is.null(names) ) {
    cat(names,sep,sep="",file=file,append=append)
    append <- TRUE
  }
  suppressWarnings(write.table(x,file=file,append=append,sep=sep,row.names=row.names,col.names=col.names,dec=dec,
                               quote=quote,eol=eol,na=na,qmethod=qmethod,fileEncoding=fileEncoding))
  if (newline>0) {
    cat( paste(rep("\n",newline),collapse=""), file=file, append=TRUE )
  }
}
list.append <- function( LIST, item, name=NULL )
{
  length(LIST) <- length(LIST)+1
  LIST[[length(LIST)]] <- item
  names(LIST)[length(LIST)] <- name
  LIST
}
#' @export
colGradient <- function( cols, length, cmax=255 )
{
  ## e.g., to create a white-to-red gradient with 10 levels
  ##
  ##   colGradient(cols=c('white','red'),length=10)
  ##
  ## or, to create a blue-to-white-to-red gradients with 9 colors (4 blue's, white, 4 red's)
  ##
  ##   colGradient(cols=c('blue','white','red'),length=9)
  ##
  ramp <- colorRamp(cols)
  rgb( ramp(seq(0,1,length=length)), maxColorValue=cmax )
}
genName <- function( stub=NULL )
{
  while ( TRUE ) {
    fname <- paste(stub,format(Sys.time(), "%Y%m%d_%H%M%S"),sep="")
    if ( file.access(fname)!=0 ) break
  }
  fname
}
## map any type of vector (e.g., vector of strings) to a vector of
## consecutive numbers (starting at base)
##
#' @export
genIndex <- function( X, base=0, add.levels=TRUE, do.sort=FALSE )
{
  idx <- match( X, if (do.sort) sort(unique(X)) else unique(X) ) + (base-1)
  if ( add.levels )
    levels(idx) <- if (is.null(levels(X))) unique(idx) else levels(X)
  idx
}
## depends on library(xlsx)
##
## write a list of matrices (or data.frames) to a multi-sheet excell workbook
##
#' @export
multi.write.xlsx <- function
(
 LIST,       # named list of data.frames; list names will be used as sheet names
 file,       # output file
 xlsx2=TRUE, # efficient (xlsx2) vs. inefficient (xlsx) with large data frames
 ...
 )
{
  for ( i in 1:length(LIST) )
  {
    if (xlsx2)
      write.xlsx2(LIST[[i]],sheetName=names(LIST)[i],file=file,append=i>1,...)
    else
      write.xlsx(LIST[[i]],sheetName=names(LIST)[i],file=file,append=i>1,...)
  }
}
## reads a multi-sheet excel workbook into a list of data.frames
## presumes the names of the sheets are known (snames)
##
#multi.read.xlsx <- function( file, snames )
#{
#  OUT <- lapply(snames, function(X) read.xlsx2(file,sheetName=X,check.names=FALSE))
#  names(OUT) <- snames
#  OUT
#}
## hierarchical clustering with optimal leaf ordering (needs package 'cba')
## require(cba)
##
hcopt <- function(d, HC=NULL, method = "ward.D", members = NULL)
{
  if ( is.null(HC) ) {
    HC <- hclust(d,method=method,members=members)
  }
  ## using optimal ordering only with more than two items
  if ( length(HC$order)>2 ) {
      ORD <- order.optimal(d,merge=HC$merge)
      HC$merge <- ORD$merge
      HC$order <- ORD$order
  }
  HC
}
## remove unused levels from factors and data frames (usually as a consequence of factor subsetting)
##
pruneFactor <- function( X )
{
  if ( !is.factor(X) ) return(X)
  factor(X,levels=levels(X)[levels(X) %in% unique(X)])
}
pruneFactors <- function( DF )
{
  if ( !is.data.frame(DF) ) stop( "!is.data.frame(DF)" )
  for ( i in 1:ncol(DF) )
    DF[,i] <- pruneFactor(DF[,i])
  DF
}
many2one <- function( cls ) many2one.cls(cls) # just for backward compatibility
many2one.cls <- function( cls )
{
  if ( is.null(dim(cls)) || ncol(cls)==1 ) {
    cls2 <- as.numeric(match(drop(cls),sort(unique(drop(cls)))))
    levels(cls2) <- levels(cls)
    return(cls2)
  }
  # ELSE ..
  #
  levs <- sort(apply(unique(cls),1,paste,collapse="."))
  cls.new <- as.numeric( match( apply(cls,1,paste,collapse="."), levs) )
  levels(cls.new) <- levs
  cls.new
}
