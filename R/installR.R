## various packages (more than you probably need)
#' @import devtools
#' @import BiocManager
#' @export
installR <- function() {
     pkgs <- c("cba",
               "caret",
               "combinat",
               "data.table",
               "dendextend",
               "devtools",
               "dynamicTreeCut",
               "e1071",
               "ff",
               "ggdendro",
               "ggplot2",
               "glmnet",
               "gplots",
               "gridExtra",
               "heatmap.plus",
               "mclust",
               "msigdbr",
               "openxlsx",
               "pamr",
               "pheatmap",
               "pkgdown",
               "pROC",
               "Rmisc",
               "randomForest",
               "rgl",
               "rjson",
               "rmarkdown",
               "roxygen2",
               "statmod",
               "tsne",
               "umap",
               "VennDiagram"
               )
     install.packages(pkgs,repos="http://cran.r-project.org")

     ## packages better installed directly from github
     devtools::install_github("montilab/hypeR")
     devtools::install_github("montilab/vennr")

     ## install bioconductor and needed packages
     if (!requireNamespace("BiocManager", quietly = TRUE))
         install.packages("BiocManager",repos="http://cran.r-project.org")
     BiocManager::install() ## install basic distribution
                            ## install packages of interest
     BiocManager::install(c("GEOquery","Biobase","biomaRt","ROC","ConsensusClusterPlus","ASSIGN",
                            "limma","DESeq2","edgeR","GSVA","RNASeqPower","pdInfoBuilder","ComplexHeatmap"))
}
