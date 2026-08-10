# Generates a 3D spatial feature plot from a SpaMTP object

This function generates an interactive plotly 3D spatial plot of 2
specified features. These features can include a combination of
metabolites and/or genes expression values, or numerical meta.data
values. Features must be numeric values. Please see
` SpaMTP Additional Features` vignette for instructions on how to plot
catagorical data, or \> 2 features.

## Usage

``` r
Plot3DFeature(
  data,
  features,
  assays = c("SPT", "SPM"),
  slots = "counts",
  between.layer.height = 100,
  names = NULL,
  size = 3,
  col.palette = "Reds",
  x.axis.label = "x",
  y.axis.label = "y",
  z.axis.label = "z",
  show.x.ticks = FALSE,
  show.y.ticks = FALSE,
  show.z.ticks = FALSE,
  show.image = NULL,
  plot.height = 800,
  plot.width = 1500,
  image.sf = "lowres",
  downscale.image = NULL
)
```

## Arguments

- data:

  A SpaMTP Seurat Object.

- features:

  A character vector specifying the features to be plotted.

- assays:

  A character vector specifying the Seurat Object assays to be used for
  plotting (default = c("SPT", "SPM")).

- slots:

  A character vector specifying the assay slots to be used for plotting
  (deafult = "counts").

- between.layer.height:

  A numeric value specifying the height between layers (default = 100).

- names:

  A character vector specifying custom names for the features. If NULL
  will use feature names (default = NULL).

- size:

  A numeric value specifying the size of markers in the plot (default =
  3).

- col.palette:

  Character list defining the colour palettes to use for plotting. If
  only 1 is supplied than it will be used for both variables. Possible
  palettes to use are:
  Blackbody,Bluered,Blues,Cividis,Earth,Electric,Greens,Greys,Hot,Jet,Picnic,Portland,Rainbow,RdBu,Reds,Viridis,YlGnBu,YlOrRd
  (default = "Reds").

- x.axis.label:

  A character string specifying the label for the x-axis (default =
  "x").

- y.axis.label:

  A character string specifying the label for the y-axis (default =
  "y").

- z.axis.label:

  A character string specifying the label for the z-axis (default =
  "z").

- show.x.ticks:

  A logical value specifying whether to show ticks on the x-axis
  (default = FALSE).

- show.y.ticks:

  A logical value specifying whether to show ticks on the y-axis
  (default = FALSE).

- show.z.ticks:

  A logical value specifying whether to show ticks on the z-axis
  (default = FALSE).

- show.image:

  Character string specifying the image name to plot. If NULL then no
  image is plot (default = NULL).

- plot.height:

  Numeric value defining the height of the returned plot (default =
  800).

- plot.width:

  Numeric value defining the width of the returned plot (default =
  1500).

- image.sf:

  Character string defining the image scalefactor to use (default =
  "lowres").

- downscale.image:

  Numeric value defining the amount to downscale the image by. This
  parameter is only necessary when show.image != NULL. Downscaling the
  image is not necessary, however when rendering plotly to HTML
  excessively large images can cause issues whereby setting a
  downscaling value can improve this (default = NULL).

## Value

A 3D Plotly plot

## Examples

``` r
# Plot3DFeature(data = my_data, features = c("gene1", "gene2"), assays = c("SPT", "SPM"), show.image = "slice1")
```
