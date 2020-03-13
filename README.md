# BS831
Course Materials for Genomics Data Mining

## Documentation

Please visit <https://montilab.github.io/BS831/>

## Requirements

## Installation

Install the the package from Github.

```r
devtools::install_github("montilab/BS831")
```

Developer notes on rebuilding documentation.
```r
library(pkgdown)

# Use lazy=FALSE to rebuild entire site
build_site(pkg=".", lazy=TRUE)

# Only rebuild changed articles
build_articles(pkg=".", lazy=TRUE)

# Rebuild one article
build_article(name, pkg=".", lazy=FALSE)
```