% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PathwayNetworks.R
\name{PathwayNetworkPlots}
\alias{PathwayNetworkPlots}
\title{Constructs an interactive network for exploring spatial metabolomics and transcriptomics data.}
\usage{
PathwayNetworkPlots(
  SpaMTP,
  ident,
  regpathway,
  DE.list,
  selected_pathways = NULL,
  path = getwd(),
  SM_slot = "counts",
  ST_slot = "counts",
  colour_palette = NULL,
  SM_assay = "SPM",
  ST_assay = "SPT",
  analyte_types = c("genes", "metabolites"),
  image = "slice1",
  verbose = T
)
}
\arguments{
\item{SpaMTP}{A \code{SpaMTP} Seurat object containing spatial metabolomics (SM) and/or spatial transcriptomics (ST) data. If SM data is included, it must be annotated using the \code{SpaMTP::AnnotateSM()} function.}

\item{ident}{A character string specifying the cluster identifier used to group regions, corresponding to a column name in the \code{SpaMTP@meta.data} slot.}

\item{regpathway}{A dataframe output from the \code{SpaMTP::FindRegionalPathways()} function, containing identified regional pathways.}

\item{DE.list}{A list containing differential expression results from the \code{FindAllMarkers()} function, with items matching the order of the \code{analyte_types} argument.}

\item{selected_pathways}{A character vector specifying the names or IDs of pathways used to construct the network (e.g., \code{c("Amino acid metabolism", "WP1902", "Aspartate and asparagine metabolism")}). This argument is not case-sensitive.}

\item{path}{The directory to save the output. If not provided, the default is the current working directory.}

\item{SM_slot}{The slot name containing the spatial metabolomics assay matrix (default: \code{"counts"}).}

\item{ST_slot}{The slot name containing the spatial transcriptomics assay matrix (default: \code{"counts"}).}

\item{colour_palette}{The color palette used to plot the spatial image in the output HTML file. Default: \code{grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100)}.}

\item{SM_assay}{A character string specifying the assay name for spatial metabolomics data in \code{SpaMTP} to extract intensity values (default: \code{"SPM"}).}

\item{ST_assay}{A character string specifying the assay name for spatial transcriptomics data in \code{SpaMTP} to extract RNA count values (default: \code{"SPT"}).}

\item{analyte_types}{A subset of \code{c("genes", "metabolites")}. Can be \code{c("genes")}, \code{c("metabolites")}, or both.}

\item{image}{Character string specifying which image stored within the SpaMTP object to use for plotting (default = "slice1").}

\item{verbose}{A logical value indicating whether to display detailed messages during execution (default: \code{FALSE}).}
}
\value{
An interactive HTML file visualizing the network structure of the specified pathways.
}
\description{
Constructs an interactive network for exploring spatial metabolomics and transcriptomics data.
}
\examples{
#PathwayNetworkPlots(SpaMTP, ident = "Custom_ident", regpathway = regpathway, DE.list = DE.list, selected_pathways = selected_pathways)
}
