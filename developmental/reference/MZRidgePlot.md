# Generates a ridge plot of spatial metabolic intensity data

Generates a ridge plot of spatial metabolic intensity data

## Usage

``` r
MZRidgePlot(
  seurat.obj,
  group.by = NULL,
  mzs = NULL,
  assay = "Spatial",
  slot = "counts",
  title = "RidgePlot",
  x.lab = "intensity",
  y.lab = "var",
  bottom.cutoff = NULL,
  top.cutoff = NULL,
  bins = 1000,
  log.data = FALSE,
  cols = NULL,
  verbose = FALSE
)
```

## Arguments

- seurat.obj:

  Seurat object containing the metabolomic intensity data.

- group.by:

  Character string specifying the meta.data column to group by (default
  = NULL).

- mzs:

  Vector of characters defining which features (m/z's) to label on the
  plot. If `NULL` no features will be labeled (default = NULL).

- assay:

  Character string defining the name of the Seurat Object assay to pull
  the corresponding intensity data from (default = "Spatial").

- slot:

  Character string defining the name of the slot within the Seurat
  Object assay to pull the corresponding intensity data from (default =
  "counts").

- title:

  Character string of the plot title (default = "RidgePlot").

- x.lab:

  Character string of the x-axis label (default = "var").

- y.lab:

  Character string of the y-axis label (default = "intensity").

- bottom.cutoff:

  Numeric value defining the percent of data to exclude for the lower
  end of the distribution. A bottom.cutoff = 0.05 will remove the bottom
  5% of data point (default = NULL).

- top.cutoff:

  Numeric value defining the percent of data to exclude for the upper
  end of the distribution. A top.cutoff = 0.05 will remove the top 5% of
  data point (default = NULL).

- bins:

  number of bins to group

- log.data:

  Boolean value indicating whether to log transform the y-axis values
  (default = FALSE).

- cols:

  Vector of strings defining the colours to use for plotting. This
  vector should match the length of unique groups (default = NULL).

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = FALSE).

## Examples

``` r
# MZRidgePlot(SeuratObj, group.by = "sample")
```
