% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialAnalysis.R
\name{RunMetabolicPCA}
\alias{RunMetabolicPCA}
\title{Generates PCA analysis results for a SpaMTP Seurat Object}
\usage{
RunMetabolicPCA(
  SpaMTP,
  npcs = 30,
  variance_explained_threshold = 0.9,
  assay = "SPM",
  slot = "counts",
  show_variance_plot = FALSE,
  bin_resolution = NULL,
  resolution_units = "ppm",
  bin_method = "sum",
  reduction.name = "pca",
  verbose = TRUE
)
}
\arguments{
\item{SpaMTP}{SpaMTP Seurat class object that contains spatial metabolic information.}

\item{npcs}{is an integer value to indicated preferred number of PCs to retain (default = 30).}

\item{variance_explained_threshold}{Numeric value defining the explained variance threshold (default = 0.9).}

\item{assay}{Character string defining the SpaMTP assay to extract intensity values from (default = "SPM").}

\item{slot}{Character string defining the assay slot containing the intensity values (default = "counts").}

\item{show_variance_plot}{Boolean indicating weather to display the variance plot output by this analysis (default = FALSE).}

\item{reduction.name}{Character string indicating the name associated with the PCA results stored in the output SpaMTP Seurat object (default = "pca").}

\item{verbose}{Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).}

\item{resampling_factor}{is a numerical value > 0, indicate how you want to resample the size of original matrix (default = 1).}

\item{byrow}{is a boolean to indicates whether each column of the matrix is built byrow or bycol (default = FALSE).}
}
\value{
SpaMTP object with pca results stored in the
}
\description{
Generates PCA analysis results for a SpaMTP Seurat Object
}
\examples{
# HELPER FUNCTION
}
