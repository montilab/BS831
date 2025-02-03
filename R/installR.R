## You can source this function as follows:
## source("https://raw.githubusercontent.com/montilab/BS831/master/R/installR.R")

## various packages (more than you probably need)
#' @import devtools
#' @import BiocManager
#' @export
installR <- function(
    add = FALSE,           # additional packages not needed by BS831
    install_shiny = FALSE  # ditto
)
{
    ## core CRAN packages
    core_cran_pkgs <- c(
        "caret",
        "cba",
        "circlize",
        "dplyr",
        "eulerr",
        "ggdendro",
        "ggplot2",
        "ggpmisc",
        "gridExtra",
        "gtable",
        "heatmap.plus",
        "magrittr",
        "manhattanly",
        "mclust",
        "openxlsx",
        "pROC",
        "pheatmap",
        "plotly",
        "RColorBrewer",
        "reshape2",
        "scales",
        "tidyr",
        "tsne",
        "umap",
        "VennDiagram",
        "vennr",
        "venn")
    ## additional optional CRAN packages
    add_cran_pkgs <- c(
        "box",
        "Cairo",
        "combinat",
        "data.table",
        "dendextend",
        "devtools",
        "dynamicTreeCut",
        "e1071",
        "ff",
        "gee",
        "ggcorplot",
        "gitlabr",
        "glmnet",
        "gplots",
        "hierarchicalSets",
        "leiden",
        "lme4",
        "lmerTest",
        "markdown",
        "matrixStats",
        "mediation",
        "MEGENA",
        "misc3d",
        "missForest",
        "msigdbr",
        "pamr",
        "pbapply",
        "pkgdown",
        "PMCMRplus", # jonckheereTest
        "psych",     # matrix correlation w/ p-values and q-values
        "Rmisc",
        "randomForest",
        "randomForestSRC",
        "reactable",
        "rgl",
        "rjson",
        "rmarkdown",
        "rmeta",
        "roxygen2",
        "qdapTools",
        "scCustomize",
        "shinyjs",
        "shinythemes",
        "SILGGM",
        "simr",
        "statmod",
        "survminer",
        "tidygraph",
        "tidyverse",
        "vip"
    )
    ## Bioconductor packages
    core_bioC_pkgs <- c(
        "affy",
        "ASSIGN",
        "Biobase",
        "biomaRt",
        "ComplexHeatmap",
        "ConsensusClusterPlus",
        "DESeq2",
        "GEOquery",
        "GSVA",
        "edgeR",
        "limma",
        "multtest")
    ## additional optional BioC packages
    add_bioC_pkgs <- c(
        "BioinformaticsFMRP/TCGAbiolinksGUI.data",
        "BioinformaticsFMRP/TCGAbiolinks",
        "bumphunter",
        "dittoSeq"
        "EnhancedVolcano",
        "EnsDb.Hsapiens.v86",
        "flashClust",
        "GO.db",
        "gprofiler2",
        "impute",
        "MultiAssayExperiment",
        "pdInfoBuilder",
        "RNASeqPower",
        "ROC",
        "sesameData",
        "sesame",
        "Seurat",
        "IlluminaHumanMethylation450kanno.ilmn12.hg19"
    )
    ## install CRAN packages
    install.packages(core_cran_pkgs,repos="http://cran.r-project.org")
    if (add) {
        install.packages(add_cran_pkgs,repos="http://cran.r-project.org")
    }
    ## install Bioconductor packages
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos="http://cran.r-project.org")
    BiocManager::install()                ## install basic distribution
    BiocManager::install(core_bioC_pkgs)  ## additional packages
    if (add) {
        BiocManager::install(add_bioC_pkgs)
    }
    ## packages better installed directly from github
    devtools::install_github("montilab/hypeR")
    devtools::install_github("montilab/vennr")
    #devtools::install_github("montilab/ConAn")
    #devtools::install_github("montilab/BS831")
}
