
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SpaMTP <img src="man/figures/logo.png" align="right" height="100" alt="" />

<!-- badges: start -->


## New *R*-based User-Friendly Spatial Metabolomic, Transcriptomic, and Proteomic Data Analysis Tool

<br>

<!-- badges: end -->


SpaMTP is an *R* based wrapper package for [*Seurat*](https://satijalab.org/seurat/) for the analysis of spatial metabolomic data. This user-freindly package contains various function for the analysis, integration and visalisation of multi-modal datasets, in particularly focusing on metabolomic/transcriptomics integration.
This package includes various functions to preform pre-processing,
spatial visualisation, various down-stream biological centered analyses,
data integration and data export of Spatial Metabolomic data. In
addition, this package has the ability to use both Cardinal (Spatial
Metabolomic based software) and Seurat (Spatial Transcriptomic based
software) functions.

<br>

## Installation

You can install the current version of SpaMTP from
[GitHub](https://github.com/) with:

``` r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("BiomedicalMachineLearning/SpaMTP")
```

For tutorials and more information please visit the [SpaMTP website](https://genomicsmachinelearning.github.io/SpaMTP/)
