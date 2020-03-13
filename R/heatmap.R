suppressPackageStartupMessages(require(heatmap.plus))

ramp.wr <- colorRamp(c( "white","red"))
palette.wr <- rgb( ramp.wr(seq(0, 1, length = 12)), max = 255)
ramp.br <- colorRamp(c( "blue","white","red"))
palette.br <- rgb( ramp.br(seq(0, 1, length = 14)), max = 255)
old.palette.br <- c("#0000FF", "#4040FF", "#7070FF", "#8888FF",
                    "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0",
                    "#FF7080", "#FF5A5A", "#FF4040", "#FF0000",
                    "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D")

my.heatmap <- function
(dat,
 col=palette.br,
 rng=c(NA,NA),
 Rowv=NA,
 Colv=NA,
 jpg=NULL,
 do.x11=(names(dev.cur())=="null device" || names(dev.cur())=="X11") && do.palette,
 negative=F,
 do.rev=F,
 color.code=FALSE,
 do.palette=F,
 pal.names=NULL,
 ...
 )
{
  if (is.null(col)) col <- heat.colors(12)
  rng[is.na(rng)] <- c(min(dat,na.rm=T),max(dat,na.rm=T))[is.na(rng)]
  pal <- seq(rng[1], rng[2], (rng[2]-rng[1])/(length(col)-1) )
  labCol <- round(pal,2); if (negative) labCol <- rev(labCol)

  #heatmap( if (negative) max(dat,na.rm=T)-dat else dat, col=col, Rowv=Rowv, Colv=Colv, ... )
  heatmap.plus( if (negative) max(dat,na.rm=T)-dat else dat, col=col, Rowv=Rowv, Colv=Colv, ... )

  if ( color.code ) {
    CSCn <- apply(CSC01,2,function(z) match(z,unique(z)))
    CSCn[,-1] <- t(t(CSCn[,-1,drop=FALSE])+cumsum(apply(CSCn[,-ncol(CSCn),drop=FALSE],2,max)))
    CSCncol <- unlist(apply(CSC01,2,unique))
  }
  if (do.x11) {
    x11()
  }
  else if ( !is.null(jpg) && do.palette ) {
    dev.off()
    jpeg(jpg)
  }
  if (do.palette)
    heatmap( rbind(pal,pal), labRow=NA, col=col, Rowv=NA, Colv=NA,cexCol=1.00,scale="n",
             labCol=if(is.null(pal.names)) labCol else pal.names, margins=c(40,5) )
  if ( !is.null(jpg) && do.palette ) dev.off()
}
sd.map <- function(dat,max.sd=3,n.sd=1,robust=FALSE)
{
  SD <- apply(dat,1,if(robust) mad else sd)
  MN <- apply(dat,1,if(robust) median else mean)
  TMP <- 
    sapply(1:nrow(dat),function(i)
             cut(dat[i,],
                 breaks=unique(c(-Inf,MN[i]+seq(from=-(max.sd*SD[i]-SD[i]/(2*n.sd)),to=-(SD[i]/(2*n.sd)),by=SD[i]/n.sd),
                   MN[i]+seq(from=SD[i]/(2*n.sd),to=max.sd*SD[i]-SD[i]/(2*n.sd),by=SD[i]/n.sd),+Inf)),
                 include.lowest=T,labels=FALSE))
  dat01 <- dat
  dat01[,] <- t(TMP)
  dat01
}
value.cap <- function(dat,qnt=.95)
{
  cap <- quantile(as.vector(dat),probs=qnt)
  dat[dat>cap] <- cap
  dat
}
clustColors <- function(clust,
                        ncuts,
                        COL=c('magenta','darkgreen','blue','orange','gray','black'),
                        do.matrix=TRUE)
{
    CC <- COL[cutree(clust,k=ncuts)];
    if ( do.matrix ) {
      CC <- cbind(CC,CC)
      colnames(CC) <- c('','cluster')
    }
    return(CC)
}
 
