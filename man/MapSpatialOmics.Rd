% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialMultiomicAlignment.R
\name{MapSpatialOmics}
\alias{MapSpatialOmics}
\title{Maps Spatial Metabolomic (MALDI) data to corresponding Spatial Transcriptomics data and coordinates.}
\usage{
MapSpatialOmics(
  SM.data,
  ST.data,
  ST.hires = FALSE,
  SM.assay = "Spatial",
  ST.assay = "Spatial",
  SM.fov = "fov",
  ST.image = "slice1",
  ST.scale.factor = "hires",
  SM.pixel.width = NULL,
  overlap.threshold = 0.2,
  annotations = TRUE,
  add.metadata = TRUE,
  merge.unique.metadata = TRUE,
  map.data = FALSE,
  new_SPT.assay = "SPT",
  new_SPM.assay = "SPM",
  verbose = FALSE
)
}
\arguments{
\item{SM.data}{A SpaMTP Seurat object representing the Spatial Metabolomics data.}

\item{ST.data}{A Seurat object representing the Spatial Transcriptomics data.}

\item{ST.hires}{Boolean string defining if the ST data is at a higher resolution compared to the SM pixel data. For example, generally Visium data will be lower res whereas Xenium/single-cell resolution spatial data will be a higher resolution (default = FALSE).}

\item{SM.assay}{Character string defining the Seurat assay that contains the annotated counts and metadata corresponding to the m/z values (default = "Spatial").}

\item{ST.assay}{Character string specifying the current assay to use to extract transcriptional data from (default = "Spatial").}

\item{SM.fov}{Character string of the image fov associated with the spatial metabolomic data (default = "fov").}

\item{ST.image}{Character string matching the image name associated with the ST data such as 'fov' or 'slice1' object (default = "slice1").}

\item{ST.scale.factor}{Character string defining the image resolution associated with the Visium image pixel data. If \code{NULL} the full-res coordinates will be used and no scaling will be performed. Note: This parameter is only required for aligning lowres ST data (default = "hires").}

\item{SM.pixel.width}{Numeric value defining the width of each SM pixel. If set to \code{NULL}, the median pixel width will be calculated based on the distance between each pixel (default = NULL).}

\item{overlap.threshold}{Numeric value defining the overlap proportion threshold for a SM pixel to be associated with a ST spot. For example, if res_increase = 0.2 then pixels that have at least 20\% area overlap with the respective visium spot will be assigned a match. Note: This parameter is only required for aligning lowres ST data (default = 0.2).}

\item{annotations}{Boolean value indicating if the Spatial Metabolomics (MALDI) Seurat object contains annotations assigned to m/z values (default = TRUE).}

\item{add.metadata}{Boolean defining whether to add the current metadata stored in the SM object to the new mapped multi-omic SpaMTP object (default = TRUE)}

\item{merge.unique.metadata}{Boolean indicating whether to summaries duplicated metadata terms to only store unique values in the metadata. Note: This parameter is only required for aligning lowres ST data, and \code{add.metadata} must be set to \code{TRUE} for this functionality to be implemented (default = TRUE).}

\item{map.data}{Boolean indicating whether to map normalised/additional data stored in the \code{data} slot of the SpaMTP assay. Note: this process is computationally expensive with large datasets (default = FALSE).}

\item{new_SPT.assay}{Character string defining the assay name of the new overlaid SpaMTP Seurat object containing all updated transcriptomics data (default = "SPT").}

\item{new_SPM.assay}{Character string defining the assay name of the new overlaid SpaMTP Seurat object containing all updated metabolomic data (default = "SPM").}

\item{verbose}{Boolean value indicating whether to print informative progression update messages and progress bars (default = TRUE).}
}
\value{
A SpaMTP Seurat object with Spatial Metabolomic data mapped to equivalent Spatial Transcripomics coordinates (Visium spots/Xenium cells).
}
\description{
Maps Spatial Metabolomic (MALDI) data to corresponding Spatial Transcriptomics data and coordinates.
}
\examples{

## Mapping MALDI data to equivalent Visium spots
# MapSpatialOmics(VisiumObj, SeuratObj, ST.scale.factor = "hires", SM.assay = "Spatial", ST.assay = "Spatial")

#' ## Mapping MALDI data to equivalent Xenium cells
# MapSpatialOmics(VisiumObj, SeuratObj, SM.assay = "Spatial", ST.assay = "Xenium")
}
