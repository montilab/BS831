## This script takes a signature, consisting of a set of gene
## identifiers, and uses ASSIGN to rank expression profiles based on
## the given signature 'activitiy'
#'
#' @import ASSIGN
#' @import Biobase
#' @import RColorBrewer
#' @import gplots
#' @export
run_assign <- function(
    ES,        # ExpressionSet object 
    eSig,      # Expression signature
    oDir,      # directory where ASSIGN output will be stored
    iter=3000, # Number of iterations of the MCMC chain run
    beta=TRUE  # Whether or not to use mixture beta when running assign
)
{
  if ( (noverlap <- sum(fData(ES)$gene_symbol %in% eSig))<1 )
    stop( "insufficient signature genes in input data:", noverlap )
  
  cat("Running with signature:                  ", substitute(eSig), "\n")
  cat("Number of genes in signature:            ", length(eSig), "\n")
  cat("Overlap between signature and ESet genes:", noverlap,"\n")
  cat("Run ASSIGN with mixture modelling beta:  ", beta," ..\n")
  cat("Saving output to directory:              ", eval(oDir),"\n")
  
  ## Create output directory
  dir.create(oDir,showWarnings=F)

  ## extract expression data limited to signature genes
  mat <- Biobase::exprs(ES)[fData(ES)$gene_symbol %in% eSig,]
  rownames(mat) <- fData(ES)$gene_symbol[fData(ES)$gene_symbol %in% eSig]
  gene_list <- rownames(mat)
  
  ## Take log transform (if only normalized RNASeq data), adding pseudocount to avoid NA's
  mat <- log2(mat+1)
  
  ## Initialize assign
  processed.data <- assign.preprocess(trainingData=NULL,
                                      testData=mat,
                                      trainingLabel=NULL,
                                      geneList=list(gene_list),
                                      n_sigGene=NA,
                                      theta0=0.05,
                                      theta1=0.9)
  
  ## Run mcmc
  mcmc.chain <- assign.mcmc(Y=processed.data$testData_sub, 
                            Bg = processed.data$B_vector, 
                            X=processed.data$S_matrix, 
                            Delta_prior_p = processed.data$Pi_matrix, 
                            iter = iter, 
                            adaptive_B=TRUE, 
                            adaptive_S=TRUE, 
                            mixture_beta=beta)

  ## generate summary
  mcmc.pos.mean <- assign.summary(test=mcmc.chain, burn_in=1000, 
                                  iter=iter, adaptive_B=TRUE, 
                                  adaptive_S=TRUE,mixture_beta=beta)

  ## save summary to output
  assign.output(processed.data=processed.data, 
                mcmc.pos.mean.testData=mcmc.pos.mean, 
                trainingData=NULL, 
                testData=mat, 
                trainingLabel=NULL, 
                testLabel=NULL, 
                geneList=list(gene_list), 
                adaptive_B=TRUE, 
                adaptive_S=TRUE, 
                mixture_beta=beta, 
                outputDir=oDir)
  
  output.data <- list(processed.data = processed.data,
                      mcmc.pos.mean.testData = mcmc.pos.mean)
  save(output.data, file = paste(oDir,"output.rda",sep='/'))

  cat("Done!\n")
}
########################################################################
##                           SUPPORT FUNCTIONS                         #
########################################################################

## function: SCALE 01
##
scale01 <- function(x, low = min(x), high = max(x)) {
  x <- (x - low)/(high - low)
  x
}
## function: SIG FRAME
##
sig_frame <- function(L) {
  
  pad.na <- function(x,len) {
    c(x,rep("",len-length(x)))
  }
  maxlen <- max(sapply(L,length))
  
  d <- do.call(data.frame,lapply(L,pad.na,len=maxlen))
  colnames(d) <- names(L)
  return(d)
}
## function: FETCH GENES
##
fetch_genes <- function(ESet,
                        output.data,
                        gene_list)
{
  gene_scores <- as.vector(output.data$mcmc.pos.mean.testData$S_pos)
  names(gene_scores) <- names(output.data$processed.data$Pi_matrix)
  
  pos_score_genes <- names(gene_scores)[gene_scores > 0]
  neg_score_genes <- names(gene_scores)[gene_scores < 0]
  
  cat ("Percentage genes depending on positivity/negativity of scores:\n\n")
  cat("Positive: ",round((length(pos_score_genes)/length(gene_scores))*100,2),"\n")
  cat("Negative: ",round((length(neg_score_genes)/length(gene_scores))*100,2),"\n\n")
  return(list("pos.genes"=pos_score_genes,"neg.genes"=neg_score_genes)) 
}
## function: ASSIGN HEATMAP
##
assign_heatmap <- function(x,
                         gene_scores,
                         sample_scores,
                         sample_probs,
                         posteriors,
                         inclusion_cutoff,
                         gene_cols=rep('green',nrow(x)),
                         drop_outliers=T,
                         fancy_order=T,
                         colSideCols=NULL,
                         main=''){
  
  
  #ordering the heatmap
  gene_order <- order(gene_scores,decreasing=T)
  sample_order <- order(sample_scores,decreasing=T)   
  x <- x[gene_order,sample_order]
  
  if (fancy_order){
    fancy_score <- colSums(sweep(x,1,gene_scores,'*'))
    fancy_order <- order(fancy_score,decreasing=T)
    x <- x[,fancy_order]
    sample_scores <- sample_scores[fancy_order]
    sample_probs <- sample_probs[fancy_order] 
  }
  rm <- rowMeans(x, na.rm = T)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd, na.rm = T)
  x <- sweep(x, 1, sx, "/")
  
  if (drop_outliers){
    x[x>4]   <- 4
    x[x<(-4)] <-  -4
  }
  
  ## defining the canvas
  lhei <- c(1.5, 4)
  lwid <- c(1.5, 4)
  lmat <- rbind(4:3, 2:1)
  
  if (!is.null(colSideCols)) {
    lmat <- rbind(lmat[1,] + 1, c(0, 1), lmat[2, ] +  1)
    lhei <- c(lhei[1], 0.2, lhei[2])
  }
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!is.null(colSideCols)) {
    par(mar = c(0.1, 1.1, 0.5, 2))
    image(cbind(1:ncol(x)), col = colSideCols[sample_order], axes = FALSE)
  }
  ## heatmap
  col <- brewer.pal(11,"RdBu")[11:1]
  breaks <- seq(min(x, na.rm = T), max(x, na.rm = T), length = length(col)+1)
  par(mar=c(2,1.1,0.5,2))
  image(t(x[nrow(x):1,]), 
        axes = FALSE, 
        xlab = '', 
        ylab = "Signature", 
        col = col,
        breaks=breaks)
  mtext("Samples",1,cex=1.5,line=0.5)
  mtext("Signature",4,cex=1.5,line=0.5)
  
  gene_cols[posteriors<inclusion_cutoff | gene_scores <0] <- 'grey'   
  bar_ord <- order(gene_scores)
  par(mar=c(0.9,4,0,0))
  barplot(-1*abs(gene_scores[bar_ord]),
          horiz=T,
          border=NA,
          axes=F,
          col=gene_cols[bar_ord],
          space=0)
  mtext('Gene scores',2,line=-1.5)
  
  ## top histogram
  par(mar=c(0,0,7,1.2))
  barplot(sort(sample_probs,decreasing=T),
          border=NA,
          col='grey',
          axes=F,
          main='',
          space=0)
  barplot(sort(sample_scores,decreasing=T),
          border=NA,
          col='purple',
          space=0,
          axes=F,
          add=T)
  axis(2,at=seq(0,1,0.2),tick=T,line=-0.8)
  axis(1,labels=F)
  legend('topright',
         legend=c('Score', 'Activation Prob.'),
         col=c('purple','grey'),
         pch=15,
         cex=1,
         pt.cex=1)
  
  title(main, cex.main = 2.5)
  
  ## key
  par(mar = c(5, 2, 2, 3), cex = 0.75)
  tmpbreaks <- breaks
  min.raw <- min(x, na.rm = TRUE)
  max.raw <- max(x, na.rm = TRUE)
  z <- seq(min.raw, max.raw, length = length(col))
  image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, xaxt = "n", yaxt = "n")
  par(usr = c(0, 1, 0, 1))
  lv <- pretty(breaks)
  xv <- scale01(as.numeric(lv), min.raw, max.raw)
  axis(1, at = xv, labels = lv)
  mtext(side = 1, "Row Z-Score", line = 2)
  
  
  h <- hist(x, plot = FALSE, breaks = breaks)
  hx <- scale01(breaks, min.raw, max.raw)
  hy <- c(h$counts, h$counts[length(h$counts)])
  lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
        col = "cyan")
  axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
  title("Color Key\nand Histogram")
  par(cex = 0.5)
}
## function: PLOT ALL
##
plotAll <- function(eSet,output.data,title,inclusion_cutoff=0.75,colSideCols=NULL,gene_list){
  
  assay.data <- log2(Biobase::exprs(eSet)+1)
  #print(head(assay.data))
  assay.data <- assay.data[fData(eSet)$gene_symbol %in% gene_list,]
  rownames(assay.data) <- fData(eSet)$gene_symbol[fData(eSet)$gene_symbol %in% gene_list]
  #print(head(assay.data))
  gene_list <- rownames(assay.data)
  
  #all
  hc01.col <- hclust(dist(t(assay.data)),method="ward.D2")
  hc01.row <- hclust(as.dist(1-cor(t(assay.data))),method="ward.D")
  
  if (!is.null(colSideCols)){
    heatmap.2(assay.data,
              scale='row',
              trace='none',
              main=title,
              Colv=as.dendrogram(hc01.col),
              Rowv=as.dendrogram(hc01.row),
              col=brewer.pal(11,"RdBu")[11:1],
              ColSideColors=colSideCols)
  }else{
    heatmap.2(assay.data,
              scale='row',
              trace='none',
              main=title,
              Colv=as.dendrogram(hc01.col),
              Rowv=as.dendrogram(hc01.row),
              col=brewer.pal(11,"RdBu")[11:1])   
  }
  
  gene_scores <- output.data$mcmc.pos.mean.testData$S_pos[match(rownames(assay.data),rownames(output.data$processed.data$Pi_matrix))]
  posteriors <- output.data$mcmc.pos.mean.testData$Delta_pos[match(rownames(assay.data),rownames(output.data$processed.data$Pi_matrix))]
  sample_scores <- output.data$mcmc.pos.mean.testData$kappa_pos
  sample_probs <- output.data$mcmc.pos.mean.testData$gamma_pos
  
  assign_heatmap(assay.data,
                 gene_scores,
                 sample_scores,
                 sample_probs,
                 posteriors,
                 inclusion_cutoff,
                 colSideCols=colSideCols,
                 main=title)
}

#Read arguements from command line
#args <- commandArgs(TRUE)
#ES <- args[1]
#sig <- args[2]
#out_dir <- args[3]

#Run ASSIGN
#run_assign(ES=ES,eSig=sig,oDir=out_dir,iter=3000)
