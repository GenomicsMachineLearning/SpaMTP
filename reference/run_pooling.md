# Pools SpaMTP Seurat object into random pools for pseudo-bulking.

Runs pooling of a SpaMTP dataset to generate pseudo-replicates for each
unique identity provided. This function is used by
[`FindAllDEMs()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/FindAllDEMs.md).

## Usage

``` r
run_pooling(data.filt, idents, n, assay, slot, seed = 1234, verbose = TRUE)
```

## Arguments

- data.filt:

  A Seurat Object containing count values for pooling.

- idents:

  A character string defining the idents column to pool the data
  against.

- n:

  An integer defining the amount of pseudo-replicates to generate for
  each sample (default = 3).

- assay:

  Character string defining the assay where the mz count data and
  annotations are stored (default = "Spatial").

- slot:

  Character string defining the assay storage slot to pull the relative
  mz intensity values from (default = "counts").

- seed:

  Numeric value used to set the seed for reproducible randomisation
  (default = 1234).

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

## Value

A SinglCellExpereiment object which contains pooled (n)-pseudo-replicate
counts data based on the Seurat Object input

## Examples

``` r
# run_pooling <- list(seuratObj, idents = "sample", n = 3, assay = "Spatial", slot = "counts")
```
