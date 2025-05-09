% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MZAnnotation.R
\name{FindDuplicateAnnotations}
\alias{FindDuplicateAnnotations}
\title{Finds if any metabolite is duplicated across multiple m/z values.}
\usage{
FindDuplicateAnnotations(data, assay = "Spatial")
}
\arguments{
\item{data}{Seurat Spatial Metabolomic Object containing annotated m/z values.}

\item{assay}{Character string defining the Seurat assay that contains the annotated metadata corresponding to the m/z values (default = "Spatial").}
}
\value{
Vector of character strings describing metabolites that are assigned to multiple m/z values
}
\description{
Finds if any metabolite is duplicated across multiple m/z values.
}
\examples{
# FindDuplicateAnnotations(SeuratObj)
}
