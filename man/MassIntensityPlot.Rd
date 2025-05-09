% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlottingSpaMTP.R
\name{MassIntensityPlot}
\alias{MassIntensityPlot}
\title{Plot mass intensity spectra}
\usage{
MassIntensityPlot(
  data,
  group.by = NULL,
  split.by = NULL,
  cols = NULL,
  assay = "Spatial",
  slot = "counts",
  label.annotations = FALSE,
  annotation.column = "all_IsomerNames",
  mz.labels = NULL,
  metabolite.labels = NULL,
  xlab = "m/z",
  ylab = "intensity",
  mass.range = NULL,
  y.lim = NULL,
  labelCex = 5,
  labelAdj = -1,
  labelOffset = 0,
  labelCol = "#eb4034",
  nlabels.to.show = NULL
)
}
\arguments{
\item{data}{Seurat object containing data to be plot.}

\item{group.by}{Character string defining the name of the meta.data column to group the data by. Results from each group will be overlayed on the one plot (default = NULL).}

\item{split.by}{Character string defining the name of the meta.data column to group and split the data by. Results from each group will be plotted individually (default = NULL).}

\item{cols}{Vector of character strings defining the colours to be used to identify each group (default = NULL).}

\item{assay}{Character string defining the relative Seurat Object assay to pull the required intensity data from (default = "Spatial").}

\item{slot}{Character string defining the relative slot from the Seurat Object assay to pull the required intensity data from (default = "counts").}

\item{label.annotations}{Boolean value defining whether to plot metabolite annotations of the supplied mz.labels or metabolite.labels on the plot (default = FALSE).}

\item{annotation.column}{Character string defining the name of the feature meta.data column which contains the stored m/z annotations (default = "all_IsomerNames").}

\item{mz.labels}{Vector of character strings defining the m/z values to display on the plot (default = NULL).}

\item{metabolite.labels}{Vector of character strings defining the metabolite names to display on the plot (default = NULL).}

\item{xlab}{Character string describing the x-axis title (default = "m/z").}

\item{ylab}{Character string describing the y-axis title (default = "intensity").}

\item{mass.range}{Vector of numeric values defining the range of m/z values to include (default = NULL).}

\item{y.lim}{Vector of numeric values defining the range of intensity values to include (default = NULL).}

\item{labelCex}{Numeric values for character expansion factor. Seen graphics::text() for more details (default = 0.7).}

\item{labelAdj}{One or two values in \verb{[0,1]} which specify the x and y adjustments for the label. Seen graphics::text() for more details (default = NULL).}

\item{labelOffset}{Value that controls the distance of the text label from the specified coordinates. Seen graphics::text() for more details (default = 0).}

\item{labelCol}{Character string defining the colour of the annotation labels (default = "#eb4034").}

\item{nlabels.to.show}{Numeric value defining the number of annotations to show per m/z (default = NULL).}
}
\value{
A mass spectrometry plot displaying mean intensity values
}
\description{
This function plots mean mass spectra intensity values for a given SpaMTP Seurat Object, and groups/splits by categories if supplied.
}
\examples{
## Plot mean of whole tissue section
# MassIntensityPlot(SeuratObj)

## Plot ssc segmentation groups on the same plot
# MassIntensityPlot(SeuratObj, group.by = "ssc")

## Plot mean of each ssc segmentation of separate plot with mz annotations
# MassIntensityPlot(SeuratObj, split.by= "ssc", mz.labels = c(329.166), plot.layout = c(5,2))

## Plot mean of each ssc segmentation of separate plot with metabolite annotations
# MassIntensityPlot(SeuratObj, split.by= "ssc", mz.labels = c(329.166), label.annotations = TRUE)
}
