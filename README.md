# BS831
Course Material for **_Genomics Data Mining and Statistics_**

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

## Course Description

<p align="justify">
The goal of this course is for the participants to develop a good
understanding and hands-on skills in the design and analysis of 
'omics' data from microarray and high-throughput sequencing
experiments, including data collection and management, statistical
techniques for the identification of genes that have differential
expression in different biological conditions, development of
prognostic and diagnostic models for molecular classification, and the
identification of new disease taxonomies based on their molecular
profile. These topics will be covered using real examples, extensively
documented hands-on's (see <a
href="https://montilab.github.io/BS831/">Markdowns' menu</a>), class
discussion and critical readings. Principles of reproducible research
will be emphasized, and participants will become proficient in the use
of the statistical language <a
href="https://cran.r-project.org/">R</a> (an advanced beginners’
knowledge of the language is expected) and associated packages
(including <a href="https://bioconductor.org/">Bioconductor</a>), and
in the use of <a href="https://rmarkdown.rstudio.com/">R markdown</a>
(and/or <a href="https://jupyter.org/">electronic notebooks</a>) for
the redaction of analysis reports.
</p>

<p align="justify">
The course is organized in seven modules covering: <b>1)</b> Introduction to
Genomics Analysis; <b>2)</b> Data Preprocessing and Quality Control; <b>3)</b>
Comparative Experiments based on Microarrays and Linear Models (LM);
<b>4)</b> Comparative Experiments based on RNA-sequencing and Generalized Linear
Models (GLM); <b>5)</b> Comparative Experiments based on Differential Enrichment
Analysis; <b>6)</b> Classification; and <b>7)</b> Clustering and Class Discovery.
</p>
