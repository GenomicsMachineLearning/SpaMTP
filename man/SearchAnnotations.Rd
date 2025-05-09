% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MZAnnotation.R
\name{SearchAnnotations}
\alias{SearchAnnotations}
\title{Find Annotation}
\usage{
SearchAnnotations(
  data,
  metabolite,
  assay = "Spatial",
  search.exact = FALSE,
  column.name = "all_IsomerNames"
)
}
\arguments{
\item{data}{Seurat Spatial Metabolomic Object containing annotated m/z values.}

\item{metabolite}{Character string of metabolite search term.}

\item{assay}{Character string defining the Seurat assay that contains the annotated metadata corresponding to the m/z values (default = "Spatial").}

\item{search.exact}{Boolean value defining if to only return m/z values which contain the exact match to the metabolite search term (default = FALSE).}

\item{column.name}{Character string defining the column name where the annotations are stored in the slot meta.data (default = "all_IsomerNames").}
}
\value{
A Data.Frame containing the peak metadata corresponding to the metabolite search term provided
}
\description{
Searches through annotated m/z values to return all which contain the metabolite search term provided
}
\examples{
# SearchAnnotations(SeuratObj, "Glucose", search.exact = TRUE)
}
