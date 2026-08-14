# Converts a SpaMTP Seurat object to a Cardinal object, including annotations and metadata

Converts a SpaMTP Seurat object to a Cardinal object, including
annotations and metadata

## Usage

``` r
ConvertSeuratToCardinal(
  data,
  assay = "Spatial",
  slot = "counts",
  run_col = NULL,
  feature.metadata = FALSE,
  verbose = TRUE
)
```

## Arguments

- data:

  Seurat object being converted.

- assay:

  Character string defining the Seurat Object assay name to pull
  intensity count data from (default = "Spatial").

- slot:

  Character string defining which slot from the Seurat Object assay to
  gather intensity data from (default = "counts").

- run_col:

  Character string describing the Seurat meta.data column where the run
  identities are stored (default = NULL).

- feature.metadata:

  Boolean value of whether the Seurat Object contains annotations stored
  in the feature metadata slot of the specified assay (default = FALSE).

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

## Value

A Cardinal object containing intensity values and feature metadata
(annotations)

## Examples

``` r
# cardinal.obj <- ConvertSeuratToCardinal(SeuratObject, feature.metadata = TRUE)
```
