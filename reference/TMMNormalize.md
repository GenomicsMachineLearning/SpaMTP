# Performs TMM normalization between categories based on a specified ident

This function is mainly used for normalising a merged SpaMTP Seurat
Object containing multiple samples.

## Usage

``` r
TMMNormalize(
  combined.obj,
  ident,
  refIdent,
  normalisation.type = "CPM",
  CPM.scale.factor = 1e+06,
  assay = "Spatial",
  slot = "counts",
  verbose = FALSE
)
```

## Arguments

- combined.obj:

  Seurat object that contains groups being normalized.

- ident:

  Character string defining the column name or ident group to normalize
  between.

- refIdent:

  Character string specifying one class/group type to use as a reference
  for TMM normalisation.

- normalisation.type:

  Character string defining the normalization method to run. Options are
  either c("CPM", "TIC", "LogNormalize") which represent counts per
  million (CPM), Total Ion Current (TIC) normalization or Log
  Normalization, respectively (default = "CPM").

- CPM.scale.factor:

  Numeric value that sets the scale factor for pixel/spot level
  normalization. Following normalization the total intensity value
  across each pixel will equal this value (default = 1e6).

- assay:

  Character string defining the name of the Seurat Object assay to pull
  the corresponding intensity data from (default = "Spatial").

- slot:

  Character string defining the name of the slot within the Seurat
  Object assay to pull the corresponding intensity data from (default =
  "counts").

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = FALSE).

## Value

Seurat object with count values normalised and corrected for between
categories

## Examples

``` r
# norm.data <- TMMNormalize(SeuratObj, ident = "samples", refIdent = "sample1", normalisation.type = "CPM")
```
