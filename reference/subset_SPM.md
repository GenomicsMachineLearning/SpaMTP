# Subsets SpaMTP Seurat Object containing FOVs

Intermediate solution to
[`subset()`](https://rdrr.io/r/base/subset.html): subset FOVs/centroids
if selected cells are NOT found in each FOV NOTE: some code parts and
args are taken from SeuratObject

## Usage

``` r
subset_SPM(
  object = NULL,
  subset = NULL,
  cells = NULL,
  idents = NULL,
  features = NULL,
  Update.slots = TRUE,
  Update.object = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- object:

  An S4 object or A `FOV` object

- subset:

  Logical expression indicating features/variables to keep

- cells:

  A vector of cells to keep; if `NULL`, defaults to all cells

- idents:

  A vector of identity classes to keep

- features:

  A vector of feature names or indices to keep

- Update.slots:

  If to update slots of an object

- Update.object:

  If to update final object, default to TRUE.

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = FALSE).

- ...:

  Arguments passed to [`subset()`](https://rdrr.io/r/base/subset.html)
  and other methods

## Value

A subset Seurat object

## Details

Function params/args:

## Examples

``` r
# sub <- subset_obt(seurat.obj, idents = "Sample1")
```
