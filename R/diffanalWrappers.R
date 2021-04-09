#################################################
## WRAPPER for running DESeq2
#################################################
##
#' @import Biobase
#' @import DESeq2
#' @export
run_deseq <- function(eset, class_id, control, treatment)
{
  ## require(DESeq2)
  control_inds <- which(pData(eset)[, class_id] == control)
  treatment_inds <- which(pData(eset)[, class_id] == treatment)
  eset.compare <- eset[, c(control_inds, treatment_inds)]

  ## make deseq2 compliant dataset
  colData <- data.frame(condition=as.character(pData(eset.compare)[, class_id]))
  dds <- DESeqDataSetFromMatrix(exprs(eset.compare), colData, formula( ~ condition))

  ## set reference to control, otherwise default is alphabetical order
  dds$condition <- factor(dds$condition, levels=c(control,treatment))

  ## run deseq2
  ## 3 steps:
  ##   1. estimate size factors
  ##   2. estimate dispersion
  ##   3. negative binomial GLM fitting and wald test
  dds_res <- DESeq(dds)
  res <- results(dds_res)
  res$dispersion <- dispersions(dds_res)
  return(res)
}
#################################################
## WRAPPER for running edgeR
#################################################
##
#' @import Biobase
#' @import edgeR
#' @export
run_edgeR <- function(eset, class_id, control, treatment)
{
  ##library(edgeR)
  control_inds <- which(pData(eset)[, class_id] == control)
  treatment_inds <- which(pData(eset)[, class_id] == treatment)

  ## make edgeR compliant dataset
  eset.compare <- eset[, c(control_inds, treatment_inds)]
  condition <- as.character(pData(eset.compare)[, class_id])

  ## run edgeR
  y <- DGEList(counts=exprs(eset.compare), group = condition)
  y <- calcNormFactors(y)
  y <- estimateGLMCommonDisp(y)
  y <- estimateGLMTrendedDisp(y)
  y <- estimateGLMTagwiseDisp(y)
  et <- exactTest(y)
  res <- topTags(et, n = nrow(eset.compare),  sort.by = "none")
  return(res)
}
#################################################
## WRAPPER for running limma
#################################################
##
## (assumes data is already log2 normalized)
##
#' @import Biobase
#' @import limma
#' @export
run_limma <- function(eset, class_id, control, treatment)
{
  control_inds <- which(pData(eset)[, class_id] == control)
  treatment_inds <- which(pData(eset)[, class_id] == treatment)

  eset.compare <- eset[, c(control_inds, treatment_inds)]
  condition <- as.character(pData(eset.compare)[, class_id])
  colData <- data.frame(condition=as.character(pData(eset.compare)[, class_id]))

  design <- model.matrix(~ 0 + factor(condition))
  colnames(design) <- levels( factor(condition))
  fit <- lmFit(eset.compare, design)
  command_str <- paste("makeContrasts(",
                       "(", treatment , "-", control, ")",
                       ",levels = design)", sep = "")

  contrast.matrix <- eval(parse(text=command_str))
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  res <- topTable(fit2, coef=1, adjust="BH", sort.by = "none", number=Inf)
  return(res)
}
#################################################
## optional: reattach empirical measurements/gene annotation to DE results
#################################################
##
#' @import Biobase
#' @export
summarize_results <- function(res, eset, class_id, control, treatment){
  control_inds <- which(pData(eset)[, class_id] == control)
  treatment_inds <- which(pData(eset)[, class_id] == treatment)
  eset.control <- eset[, control_inds]
  eset.treatment <- eset[, treatment_inds]

  eset.ordered <- match(rownames(res), rownames(eset))
  res <- cbind(res, fData(eset)[eset.ordered,])
  res$rowmeans.control <- rowMeans(exprs(eset.control)[eset.ordered,])
  res$rowmeans.treatment <- rowMeans(exprs(eset.treatment)[eset.ordered,])

  res$log2fc <- log2(res$rowmeans.treatment/res$rowmeans.control)
  index <- (res$rowmeans.treatment <= res$rowmeans.control)
  res$log2fc[index] <- -abs(res$log2fc[index])
  return(res)
}
##################################################################
## DESEQ DIFFANAL (w/ optional covariates)
##################################################################
#' @import Biobase
#' @import DESeq2
#' @export
deseqDiffanal <- function( dat, pheno, covariate=NULL, verbose=TRUE )
{
  ## BEGIN input checks
  if ( !is(dat,"ExpressionSet"))
    stop( "dat expected to be an ExpressionSet: ", class(dat) )
  if ( !(pheno %in% colnames(pData(dat))) )
    stop( "unrecognized pheno:", pheno )
  if ( !is.null(covariate) && !all(covariate %in% colnames(pData(dat))) )
    stop( "unrecognized covariate:", paste(covariate,sep=", ") )
  ## END input checks

  ## Define the regression formula.
  ## The phenotype is always last (so that DESeq2::results pulls it by default)
  formula <- {
    if (is.null(covariate))
      ## ~ pheno
      eval(paste("~",parse(text=pheno)))
    else
      ## ~ cov_1 + ... + cov_n + pheno
      eval(paste("~",paste(sapply(covariate,function(Z) parse(text=Z)),collapse=" + "),
                 "+",parse(text=pheno)))
  }
  VERBOSE(verbose, "formula: ", formula, "\n")

  ## Properly format the dataset (and remove phenotype NA's if any)
  datI <- dat[,!is.na(pData(dat)[,pheno])]
  VERBOSE( verbose, paste0(paste(table(pData(datI)[,pheno]),collapse=" vs. "), "\n") )
  dds <-DESeq2::DESeqDataSetFromMatrix(countData = exprs(datI),
                                       colData = pData(datI),
                                       design = eval(parse(text=formula)))

  ## Perform differential analysis
  dds <-DESeq2::DESeq(dds, parallel=TRUE)
  dds_res <- DESeq2::results(dds, parallel = TRUE)
  dds_res$dispersion <- DESeq2::dispersions(dds)

  ## return results
  dds_res
}
