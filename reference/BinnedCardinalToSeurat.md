# Converts a SpaMTP binned Cardinal Object into a SpaMTP Seurat Object

This function is used by
[`LoadSM()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/LoadSM.md)
when `bin_package` is set to 'SpaMTP'.

## Usage

``` r
BinnedCardinalToSeurat(
  data,
  mtx,
  multi.run = FALSE,
  assay = "Spatial",
  verbose = TRUE
)
```

## Arguments

- data:

  A Cardinal Object that is being converted into a Seurat Object.

- mtx:

  Matrix object containing

- multi.run:

  Boolean indicating if there are multiple runs within the imported
  data. If `TRUE`, an index will be added to the pixel names per run,
  and an individual FOV will be generated per run in the Seurat Object
  (default = FALSE).

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
# CardinalToSeurat(CardinalObj, run_name = "run_1")
```
