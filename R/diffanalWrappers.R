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
  cntIdx <- which(pData(eset)[, class_id] == control)
  trtIdx <- which(pData(eset)[, class_id] == treatment)
  eset.compare <- eset[, c(cntIdx, trtIdx)]

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
  cntIdx <- which(pData(eset)[, class_id] == control)
  trtIdx <- which(pData(eset)[, class_id] == treatment)

  ## make edgeR compliant dataset
  eset.compare <- eset[, c(cntIdx, trtIdx)]
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
  cntIdx <- which(pData(eset)[, class_id] == control)
  trtIdx <- which(pData(eset)[, class_id] == treatment)

  eset.compare <- eset[, c(cntIdx, trtIdx)]
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
## LIMMA DIFFANAL (w/ covariates)
#################################################
#' @import Biobase
#' @import limma
#' @export
limmaDiffanal <- function(
    eset,
    pheno,
    control,
    treatment,
    covariate=NULL,
    sort.by=NULL,
    verbose=TRUE
)
{
    ## BEGIN input checks
    if ( !is(eset,"ExpressionSet"))
        stop( "eset expected to be an ExpressionSet: ", class(eset) )
    if ( !(pheno %in% colnames(pData(eset))) )
        stop( "unrecognized pheno:", pheno )
    if ( !is.null(covariate) && !all(covariate %in% colnames(pData(eset))) )
        stop( "unrecognized covariate:", paste(covariate,sep=", ") )
    if ( !is.factor(pData(eset)[,pheno]) )
        stop( "pheno must be a factor: ", class(pData(eset)[,pheno]) )
    ## END input checks

    ## Define the regression formula.
    formula <- {
        if (is.null(covariate))
            ## ~ pheno
            paste("~",parse(text=pheno))
        else
            ## ~ cov_1 + ... + cov_n + pheno
            paste("~",paste(sapply(covariate,function(Z) parse(text=Z)),collapse=" + "),
                  "+",parse(text=pheno))
    }
    VERBOSE(verbose, "formula: ", formula, "\n")

    ## Properly format the dataset (and remove phenotype NA's if any)
    esetI <- eset[,pData(eset)[,pheno] %in% c(control,treatment)]
    pData(esetI)[,pheno] <- droplevels(factor(pData(esetI)[,pheno],levels=c(control,treatment)))
    VERBOSE( verbose, paste0(paste(table(pData(esetI)[,pheno]),collapse=" vs. "), "\n") )

    ## Run limma
    design <- model.matrix( eval(parse(text=formula)), data=pData(esetI) )
    fit <- lmFit(esetI, design)
    fit2 <- eBayes(fit)
    fit2.table <- topTable(fit2, coef=paste0(pheno,treatment), adjust="BH", number=Inf, sort.by="none")
    return(fit2.table)
}
#################################################
## optional: reattach empirical measurements/gene annotation to DE results
#################################################
##
#' @import Biobase
#' @export
summarize_results <- function(res, eset, class_id, control, treatment){
  cntIdx <- which(pData(eset)[, class_id] == control)
  trtIdx <- which(pData(eset)[, class_id] == treatment)
  eset.control <- eset[, cntIdx]
  eset.treatment <- eset[, trtIdx]

  eset.ordered <- match(rownames(res), rownames(eset))
  res <- cbind(res, fData(eset)[eset.ordered,])
  res$rowmeans.control <- rowMeans(exprs(eset.control)[eset.ordered,])
  res$rowmeans.treatment <- rowMeans(exprs(eset.treatment)[eset.ordered,])

  res$log2fc <- log2(res$rowmeans.treatment/res$rowmeans.control)
  index <- (res$rowmeans.treatment <= res$rowmeans.control)
  res$log2fc[index] <- -abs(res$log2fc[index])
  return(res)
}
