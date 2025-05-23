% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LoadData.R
\name{LoadSM}
\alias{LoadSM}
\title{Loads spatial metabolic data into a SpaMTP Seurat Object}
\usage{
LoadSM(
  name,
  path,
  mass.range = NULL,
  resolution = 10,
  units = "ppm",
  verbose = TRUE,
  assay = "Spatial",
  bin_package = "SpaMTP",
  multi.run = FALSE,
  ...
)
}
\arguments{
\item{name}{Character string of the object name. This should match the filename.}

\item{path}{Character string defining the directory path of the file. This should not include the file name.}

\item{mass.range}{Vector of numeric values indicating the mass range to use for the imported data (default = NULL).}

\item{resolution}{Numeric value defining the the accuracy to which the m/z values will be binned after reading. This value can be in either "ppm" or "mz" depending on the units type specified (default = 10).}

\item{units}{Character string defining the resolution value unit type, either c("ppm", "mz") (default = "ppm")}

\item{verbose}{Boolean indicating whether to show informative processing messages. If TRUE the message will be show, else the message will be suppressed (default = TRUE)}

\item{assay}{Character string describing the name of the new assay which stores the imported data (default = "Spatial").}

\item{bin_package}{Character string defining the package used to bin the imported data. Options are either "SpaMTP" or "Cardinal" (default = "SpaMTP").}

\item{multi.run}{Boolean indicating if there are multiple runs within the imported data. If \code{TRUE}, an index will be added to the pixel names per run, and an individual FOV will be generated per run in the Seurat Object (default = FALSE).}

\item{...}{Additional arguments passed to the \code{readMSIData} function.}
}
\value{
A new SpaMTP Seurat object contain the imported spatial metabolic intensity values
}
\description{
This function loads raw spatial metabolic data in a .imzML and .ibd format and generates a SpaMTP Seurat Object.
This function adapts the \code{readImzML} function implmented in Cardinal to correctly import large data.
}
\examples{
# data <-LoadSM(name = "run1", folder = "/Documents/SpaMTP_test_data/", mass.range = c(160,1500), resolution = 10, assay = "Spatial")
}
