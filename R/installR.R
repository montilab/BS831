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
        "ggdendro",
        "ggplot2",
        "gridExtra",
        "gtable",
        "heatmap.plus",
        "magrittr",
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
        "vennr")
    ## additional optional CRAN packages
    add_cran_pkgs <- c(
        "combinat",
        "data.table",
        "dendextend",
        "devtools",
        "dynamicTreeCut",
        "e1071",
        "ff",
        "ggcorplot",
        "gitlabr",
        "glmnet",
        "gplots",
        "markdown",
        "matrixStats",
        "MEGENA",
        "msigdbr",
        "pamr",
        "pkgdown",
        "Rmisc",
        "randomForest",
        "reactable",
        "rgl",
        "rjson",
        "rmarkdown",
        "roxygen2",
        "statmod",
        "shinyjs",
        "shinythemes",
        "survminer",
        "tidygraph",
        "tidyverse"
    )
    ## Bioconductor packages
    core_bioC_pkgs <- c(
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
        "pdInfoBuilder",
        "RNASeqPower",
        "ROC",
        "TCGAbiolinks"
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
    #devtools::install_github("montilab/BS831")
}
