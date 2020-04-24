##################################################################
# DISTANCE functions
##################################################################
euclid.dst <- function( x, y )
{
  sqrt(sum((x-y)^2))
}
spearman.dst <- function( x, y )
{
  1 - rank.cor(x,y)
}
pearson.dst <- function( x, y )
{
  1 - cor(x,y)
}
cosine.dst <- function( x, y )
{
  1 - drop( (x %*% y)/sqrt(x%*%x * y%*%y) )
}
##################################################################
# function: KNN ESTIMATE
##################################################################

knn.estimate <- function( dat, cls,
                          method=c("cosine.dst","pearson.dst","spearman.dst","euclid.dst"),
                          weight=c("distance","rank","none"),
                          k=3 )
{
  if (is.null(colnames(dat)))
    stop( "dat expected to have column names" )
  if ( nrow(dat)!=length(cls) )
    stop( "nrow(dat) must be same as length(cls)" )
  if ( !is.factor(cls) )
      stop( "cls must be a factor" )
  if ( k<1 )
    stop( "k must be > 0" )
  if ( k>ncol(dat) )
    stop( "k cannot be > ncol(dat)" )

  ## do nothing, just store data and labels
  method <- match.arg(method)
  weight <- match.arg(weight)
  list( dat=dat, cls=cls, method=method, weight=weight, k=k )
}
##################################################################
# function: KNN PREDICT
##################################################################

knn.predict <- function
(
  dat,  # sample-by-variables dataset
  model # list with elements (dat,cls,method,weight,k). See above.
)
{
  ## BEGIN checks on inputs
  ##
  if (is.null(colnames(dat)))
      stop( "dat expected to have column names" )
  if ( any(is.na( idx  <- match( colnames(model$dat), colnames(dat) ))) )
      stop( "missing variables in dat" )
  if ( !is.factor(model$cls) )
      stop( "model$cls must be a factor" )
  dist.fun <- match.fun(model$method)
  ##
  ## END checks

  ## map test data to feature space (i.e., as many columns as the predictive features)
  dat <- dat[,idx,drop=FALSE]

  ## some useful entities referred to often
  ncls <- length(levels(model$cls)) # number of classes
  levs <- levels(model$cls)         # extract class labels
  
  ## PWD is a [NxM] matrix of pairwise distances between each of N
  ## test samples and each of M training samples. That is, the i-th
  ## row in PWD is the vector of distances between the i-th test
  ## sample and each of the training samples.
  ##
  PWD <- t(apply(dat, 1, function(z) apply(model$dat, 1, dist.fun, z)))

  ## CLS and WHT are intermediate objects to store class assignments
  ## and weights of the training samples sorted in their descending
  ## order of distance from each of the test samples (got that?)
  ##
  CLS <- data.frame( NA, nrow(dat), model$k )
  WHT <- matrix(  1, nrow(dat), model$k )

  ## classify each of the test samples
  ##
  for ( i in 1:nrow(dat)  )
  {
    ## sort trainig samples by their distance from the test sample
    ord <- order(PWD[i,])

    ## pick the class assignments of the K closest training samples
    CLS[i,] <- model$cls[ord[1:model$k]]

    ## compute weights for the closest K training samples
    
    if ( model$weight=="distance" && model$method!="euclid" ) # correlation-based weigth
      WHT[i,] <- 1-PWD[i,ord[1:model$k]]
    else if ( model$weight=="rank" )                          # rank correlation weigth
      WHT[i,] <- 1/rank(PWD[i,ord[1:model$k]])
    else if ( model$weight=="distance" && model$method=="euclid" ) # euclidean distance weight
      WHT[i,] <- 1/sapply(PWD[i,ord[1:model$k]],max,1.0e-5)
  }
  CLS01 <- matrix(0,nrow(CLS),ncol(CLS))
  OUT <- matrix( NA, nrow(dat), ncls )

  #In <- rep(1,model$k)
  #Div <- WHT %*% In
  ## Compute (weighted) votes for each sample and for each class

  Div <- rowSums(WHT)
  for ( i in 1:ncls )
  {
    CLS01[,] <- 0; CLS01[CLS==levs[i]] <- 1
    OUT[,i] <- rowSums(WHT * CLS01) / Div
    #OUT[,i] <- ((WHT * CLS01) %*% In) / Div
  }
  colnames(OUT) <- paste("P(C=",levels(model$cls),")",sep="")
  rownames(OUT) <- rownames(dat)
  OUT
}
knn.summary <- function( results, cls, rnd=2 )
{
  as.data.frame(cbind(round(results,rnd),Pred=levels(cls)[apply(results,1,which.max)]))
}
