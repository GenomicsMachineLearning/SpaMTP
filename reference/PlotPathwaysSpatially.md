# Plots the spatial expression profile of a feature set corresponding to specified pathways

This function has been adapted from fgsea package. For more detail about
function inputs please visit their
[documentation](https://github.com/alserglab/fgsea/blob/master/R/geseca-plot.R#L266C1-L266C31).

## Usage

``` r
PlotPathwaysSpatially(
  pathways,
  object,
  images,
  title = NULL,
  image.alpha = 1,
  assay = SeuratObject::DefaultAssay(object),
  slot = "scale.data",
  colors = c("darkblue", "lightgrey", "darkred"),
  guide = "colourbar",
  crop = TRUE,
  min.cutoff = NA,
  max.cutoff = NA,
  ncol = NULL,
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  image.scale = "lowres",
  shape = 21,
  stroke = NA,
  interactive = FALSE,
  information = NULL,
  image.labels = NULL
)
```

## Arguments

- pathways:

  Character vector of pathway names to plot. These pathways must be
  those stored in the `chempathway` object. To check call
  `unique(chempathway$pathwayName)`.

- object:

  SpaMTP Seurat object containing the spatial data and coordinates for
  plotting.

- images:

  Character vector specifying which images to use for plotting from the
  SpaMTP Seurat Object.

- title:

  Optional title for the plot. If set to `NULL` the title will be the
  pathway name (default = NULL).

- image.alpha:

  Numeric value between 0 and 1 controlling the transparency of the
  underlying image (default = 1).

- assay:

  Character string defining the name of assay to use for data
  extraction. This `assay` should contain RAMP_IDs as feature names
  (default = `SeuratObject::DefaultAssay(object)`).

- slot:

  Character string defining the name of slot to extract data from
  (default = "scale.data").

- colors:

  Vector of colors for the gradient (default = c("darkblue",
  "lightgrey", "darkred")).

- guide:

  Character string stating the type of legend to display (default =
  "colourbar").

- crop:

  Boolean logical indicating whether to crop images (default = TRUE).

- min.cutoff:

  Numeric value defining the minimum cutoff value for color scale. If
  `NA` this value will be determined based on the minimum value of the
  dataset (default = NA).

- max.cutoff:

  Numeric value defining the maximum cutoff value for color scale. If
  `NA` this value will be determined based on the maximum value of the
  dataset (default = NA).

- ncol:

  Numeric value defining Number of columns for arranging plots. Note
  this is useful if multiple images have been supplied/containing within
  the one SpaMTP Seurat Object (default = NULL).

- pt.size.factor:

  Numeric value defining the spot point size for plotting (default =
  1.6).

- alpha:

  Numeric vector controlling the transparency of points (default =
  c(1,1).

- image.scale:

  Character string defining the image resolution to use. This can be
  either "hires" or "lowres" (default = "lowres").

- shape:

  Integer value defining the shape of points. If shape = 21, circles
  will be plotted (default = 21).

- stroke:

  Numeric value stating the stroke width for points (default = NA).

- interactive:

  Boolean logical indicating whether to create interactive plots. NOTE:
  this functionality is implemented through
  [`Seurat::SpatialFeaturePlot`](https://satijalab.org/seurat/reference/SpatialPlot.html)
  (default = FALSE).

- information:

  An optional data.frame or matrix of extra information to be displayed
  on hover. NOTE: this functionality is implemented through
  [`Seurat::SpatialFeaturePlot`](https://satijalab.org/seurat/reference/SpatialPlot.html)
  (default = NULL).

- image.labels:

  Character vector specifying optional labels for multiple images
  (default = NULL).

## Value

A ggplot object visualizing the pathway score spatially.

## Examples

``` r
#PlotPathwaysSpatially(c("Glycolysis", "Acylcarnitine 3-Butenylcarnitine", "ABC transporters"), spamtp_obj)
```
