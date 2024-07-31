
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SpaMTP !UNDER DEVELOPMENT! <img src="man/figures/logo.png" align="right" height="100" alt="" />

<!-- badges: start -->
<br>

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

<img src="inst/figures/SpaMTPSumFig.png" height="400" alt="" style="background-color: white;" />

## Another development version of the github is below
https://github.com/BCRL-tylu/SpaMTP

## Installation

You can install the current version of SpaMTP from
[GitHub](https://github.com/) with:

``` r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("BiomedicalMachineLearning/SpaMTP")
```

## Installing with conda

Below is an example of how to set up SpaMTP using a conda environment

``` console
conda create -n SpaMTP -c conda-forge r-base=4.3.3 r-essentials r-devtools
conda activate SpaMTP
```

Within the SpaMTP environment open *R* and install ***SpaMTP***
packages

``` r
devtools::install_github("BiomedicalMachineLearning/SpaMTP")
```

### Possible Installiation Errors

#### Cairo failed to install

``` console
checking if R was compiled with the RConn patch... no
checking for cairo.h... no
configure: error: Cannot find cairo.h! Please install cairo (http://www.cairographics.org/) and/or set CAIRO_CFLAGS/LIBS correspondingly.
ERROR: configuration failed for package ‘Cairo’
```

check these libraries exist:

``` console
conda install conda-forge::r-cairo
```

This should resolve any issues. Rerun:

``` r
BiocManager::install("scater")
```

#### Cardinal Failed to install

``` console
ERROR: dependency ‘EBImage’ is not available for package ‘Cardinal’

BiocManager::install("EBImage")

fftwtools.c:28:9: fatal error: fftw3.h: No such file or directory
   28 | #include<fftw3.h>
      |         ^~~~~~~~~
```

If EBImage failed to installed it ist most likely due to an issue with the installation of *fftwtools*

Try:

``` console
conda config --add channels conda-forge
conda config --set channel_priority strict

conda install r-fftwtools
```

Else try:

``` console
conda install bioconda::r-fftwtools
conda install conda-forge::fftw
```

Try in *R*:

```r

if (!require("BiocManager", quietly = TRUE)) #Check if BiocManager is installed
    install.packages("BiocManager")

BiocManager::install("fftwtools")

```

If none of the above methods resolve instiallation issues, install *Cardinal* directly through *conda*:

``` console
conda install bioconda::bioconductor-cardinal  
```

#### *rgoslin* Failed to install

If rgoslin fails to install please head to their [github](https://github.com/lifs-tools/rgoslin). Alternative, you can try installing through BiocManager or mamba shown below:

```r
BiocManager::install("rgoslin")
```
or
```console
mamba install bioconductor-rgoslin
```

For other issue please flag on github under *Issues*



## SpaMTP Vignette

To come:

``` r
## To Come ...
```
