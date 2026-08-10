# Converts a Cardinal Object into a SpaMTP Seurat Object

Converts a Cardinal Object into a SpaMTP Seurat Object

## Usage

``` r
CardinalToSeurat(
  data,
  multi.run = FALSE,
  seurat.coord = NULL,
  assay = "Spatial",
  verbose = TRUE
)
```

## Arguments

- data:

  A Cardinal Object that is being converted into a Seurat Object.

- multi.run:

  Boolean indicating if there are multiple runs within the imported
  data. If `TRUE`, an index will be added to the pixel names per run,
  and an individual FOV will be generated per run in the Seurat Object
  (default = FALSE).

- seurat.coord:

  A Data.Frame containing two columns titled 'X_new' and 'Y_new'
  specifying the pixel coordinates of each data point. This is only
  required if mapping Spatial Metabolic data with a H&E image was done
  externally, and the SM coordinates need to change to align correctly
  to the ST data. Else, set to `NULL` (default = NULL).

- assay:

  Character string containing the name of the assay (default =
  "Spatial").

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

## Value

A Seurat Object containing the mz count data of the supplied Cardinal
Object

## Examples

``` r
# CardinalToSeurat(CardinalObj, run_name = "run_1", seurat.coord = NULL)
```
