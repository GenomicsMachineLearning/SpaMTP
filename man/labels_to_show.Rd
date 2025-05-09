% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MZAnnotation.R
\name{labels_to_show}
\alias{labels_to_show}
\title{Filters the annotation list to only include the first n number of annotations per m/z}
\usage{
labels_to_show(annotation_column, n = 3)
}
\arguments{
\item{annotation_column}{Vector of the meta.data column containing the m/z annotations.}

\item{n}{Numeric value defining the number of annotations to keep (default = 3).}
}
\value{
Vector containing the first n number of annotations
}
\description{
Filters the annotation list to only include the first n number of annotations per m/z
}
\examples{
# labels_to_show(`SeuratObject[["Spatial"]]@meta.data$annotations`, n = 3)
}
