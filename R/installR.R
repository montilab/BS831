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
        "ggdendro",
        "ggplot2",
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
        "gee",
        "ggcorplot",
        "gitlabr",
        "glmnet",
        "gplots",
        "hierarchicalSets",
        "lme4",
        "lmerTest",
        "markdown",
        "matrixStats",
        "MEGENA",
        "misc3d",
        "msigdbr",
        "pamr",
        "pbapply",
        "pkgdown",
        "psych",
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
        "shinyjs",
        "shinythemes",
        "SILGGM",
        "statmod",
        "survminer",
        "tidygraph",
        "tidyverse"
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
        'flashClust',
        "GO.db",
        "gprofiler2",
        "impute",
        "MultiAssayExperiment",
        "pdInfoBuilder",
        "RNASeqPower",
        "ROC",
        "sesameData",
        "sesame",
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
