# Helper function for QC plots by generating intensity count data

Helper function for QC plots by generating intensity count data

## Usage

``` r
statPlot(
  seurat.obj,
  group.by = NULL,
  assay = "Spatial",
  slot = "counts",
  bottom.cutoff = NULL,
  top.cutoff = NULL,
  log.data = FALSE,
  verbose = FALSE
)
```

## Arguments

- seurat.obj:

  Seruat object containing the intensity data.

- group.by:

  Character string specifying the meta.data column to group by (default
  = NULL).

- assay:

  Character string defining the name of the Seurat Object assay to pull
  the corresponding intensity data from (default = "Spatial").

- slot:

  Character string defining the name of the slot within the Seurat
  Object assay to pull the corresponding intensity data from (default =
  "counts").

- bottom.cutoff:

  Numeric value defining the percent of data to exclude for the lower
  end of the distribution. A bottom.cutoff = 0.05 will remove the bottom
  5% of data point (default = NULL).

- top.cutoff:

  Numeric value defining the percent of data to exclude for the upper
  end of the distribution. A top.cutoff = 0.05 will remove the top 5% of
  data point (default = NULL).

- log.data:

  Boolean value indicating whether to log transform the y-axis values
  (default = FALSE).

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = FALSE).

## Value

A data.frame containing the relative transformed and sum counts required
for various QC plots

## Examples

``` r
# df <- statPlot(SeuratObj, group.by = "sample", bottom.cutoff = 0.05, top.cutoff = 0.05, log.data = TRUE)
```
