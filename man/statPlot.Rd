% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMPreprocessing.R
\name{statPlot}
\alias{statPlot}
\title{Helper function for QC plots by generating intensity count data}
\usage{
statPlot(
  seurat.obj,
  group.by = NULL,
  assay = "Spatial",
  slot = "counts",
  bottom.cutoff = NULL,
  top.cutoff = NULL,
  log.data = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{seurat.obj}{Seruat object containing the intensity data.}

\item{group.by}{Character string specifying the meta.data column to group by (default = NULL).}

\item{assay}{Character string defining the name of the Seurat Object assay to pull the corresponding intensity data from (default = "Spatial").}

\item{slot}{Character string defining the name of the slot within the Seurat Object assay to pull the corresponding intensity data from (default = "counts").}

\item{bottom.cutoff}{Numeric value defining the percent of data to exclude for the lower end of the distribution. A bottom.cutoff = 0.05 will remove the bottom 5\% of data point (default = NULL).}

\item{top.cutoff}{Numeric value defining the percent of data to exclude for the upper end of the distribution. A top.cutoff = 0.05 will remove the top 5\% of data point (default = NULL).}

\item{log.data}{Boolean value indicating whether to log transform the y-axis values (default = FALSE).}

\item{verbose}{Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = FALSE).}
}
\value{
A data.frame containing the relative transformed and sum counts required for various QC plots
}
\description{
Helper function for QC plots by generating intensity count data
}
\examples{
# df <- statPlot(SeuratObj, group.by = "sample", bottom.cutoff = 0.05, top.cutoff = 0.05, log.data = TRUE)
}
