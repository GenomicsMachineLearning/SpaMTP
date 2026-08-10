# Plot expression of metabolites in spatially from a Spatial SpaMTP Objects.

This is for plotting values of annotated metabolites stored in a SpaMTP
Seurat Object containing an image(i.e. H&E image). This function
inherits off Seurat::ImageFeaturePlot(). Look here for more detailed
documentation about inputs.

## Usage

``` r
SpatialMZAnnotationPlot(
  object,
  metabolites,
  plusminus = NULL,
  images = NULL,
  crop = TRUE,
  assay = "Spatial",
  slot = "counts",
  keep.scale = "feature",
  min.cutoff = NA,
  max.cutoff = NA,
  ncol = NULL,
  combine = TRUE,
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  image.alpha = 1,
  stroke = 0.25,
  interactive = FALSE,
  information = NULL,
  column.name = "all_IsomerNames",
  plot.exact = TRUE,
  verbose = TRUE
)
```

## Arguments

- object:

  Seurat Spatial Metabolomic Object to Visualise.

- metabolites:

  Vector of metabolite names to plot (e.g. c("Glucose", "Glutamine")).
  The Seurat Object provided must contain annotations in the respective
  assay metadata.

- plusminus:

  Numeric value defining the range/threshold either side of the target
  peak/peaks to be binned together for plotting (default = NULL).

- images:

  Character string of the name of the image to plot (default = NULL).

- crop:

  Boolean value indicating if to crop the plot to focus on only points
  being plotted (default = TRUE).

- assay:

  Character string indicating which Seurat object assay to pull data
  form (default = "Spatial").

- slot:

  Character string indicating the assay slot to use to pull expression
  values form (default = "counts").

- keep.scale:

  Character string describing how to handle the color scale across
  multiple plots. Check Seurat::SpatialFeaturePlot() for all options
  (default = "feature").

- min.cutoff:

  Vector of numeric value describing the minimum cutoff values for each
  m/z feature (default = NA).

- max.cutoff:

  Vector of numeric value describing the maximum cutoff values for each
  m/z feature (default = NA).

- ncol:

  Integer defining the number of columns if plotting multiple plots
  (default = NULL).

- combine:

  Boolean value stating if to combine plots into a single patchworked
  ggplot object (default = TRUE).

- pt.size.factor:

  Numeric value defining the point size for spots when plotting (default
  = 1.6).

- alpha:

  Numeric value between 0 and 1 defining the spot alpha (default = 1).

- image.alpha:

  Numeric value between 0 and 1 defining the image alpha (default = 1).

- stroke:

  Numeric value describing the width of the border around the spot
  (default. =0.25).

- interactive:

  Boolean value of if to launch an interactive SpatialDimPlot or
  SpatialFeaturePlot session, see Seurat::ISpatialDimPlot() or
  Seurat::ISpatialFeaturePlot() for more details (default = FALSE).

- information:

  An optional dataframe or matrix of extra infomation to be displayed on
  hover (default = NULL).

- column.name:

  Character string defining the column name where the annotations are
  stored in the slot meta.data (default = "all_IsomerNames").

- plot.exact:

  Boolean value describing if to only plot exact matches to the
  metabolite search terms, else will plot all metabolites which contain
  serach word in name (default = TRUE).

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

## Value

A ggplot showing the Spatial representation of expression data of
specified m/z values for Spatial data with H&E Image

## Examples

``` r
# SpatialMZAnnotationPlot(SeuratObj, mzs = c("Glucose", "Glutamine"))
# SpatialMZAnnotationPlot(SeuratObj, mzs = c("Glucose", "Glutamine"), plusminus = 0.05)
```
