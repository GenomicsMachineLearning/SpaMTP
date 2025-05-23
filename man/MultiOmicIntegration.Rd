% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiOmicAnalysis.R
\name{MultiOmicIntegration}
\alias{MultiOmicIntegration}
\title{Mult-Omic data integration}
\usage{
MultiOmicIntegration(
  multiomic.data,
  weight.list = NULL,
  reduction.list = list("spt.pca", "spm.pca"),
  dims.list = list(1:30, 1:30),
  return.intermediate = FALSE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{multiomic.data}{SpaMTP dataset contain Spatial Transcriptomics and Metabolomic datasets in two different assays}

\item{weight.list}{List containing the relative weightings for each modality, matching the reduction order. If NULL, weights will be automatically calculated else, two values must add to 1 (default = NULL).}

\item{reduction.list}{List containing character strings defining the reduction to use for each modality, in the order matching weight.list if applicable (default = list("spt.pca", "spm.pca")).}

\item{dims.list}{List containing the numeric range of principle component dimension to include for each modality (default = list(1:30,1:30)).}

\item{return.intermediate}{Boolean value indicating whether to store intermediate results in misc slot of SpaMTP Seurat class object (default = FALSE).}

\item{verbose}{Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).}

\item{...}{Additional arguments that can be parsed through Seurat's FindMultModalNeighbors function. For possible inputs please visit: https://www.rdocumentation.org/packages/Seurat/versions/5.0.3/topics/FindMultiModalNeighbors.}
}
\value{
SpaMTP Seurat class object containing a weighted nearest neighbours graph which integrates Metabolic and Transcriptomic modalities. This graph can be used for clustering.
}
\description{
This function performs multi-omic integration of Spatial Metabolomics and Spatial Transcriptomics data using Seurat's Weighted Nearest Neighbours function.
}
\examples{
# SpaMTP.obj <- MultiOmicIntegration(SpaMTP.obj, weight.list = list(0.5, 0.5), reduction.list =  list("spt.pca", "spm.pca"), dims.list = list(1:30, 1:30))
}
