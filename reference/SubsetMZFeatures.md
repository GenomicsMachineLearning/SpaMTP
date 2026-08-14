# Subset a SpaMTP Seurat Spatial Metabolomic object by a list of m/z's

Subset a SpaMTP Seurat Spatial Metabolomic object by a list of m/z's

## Usage

``` r
SubsetMZFeatures(data, features, assay = "Spatial")
```

## Arguments

- data:

  A Seurat Spatial Metabolomic Object for subsetting.

- features:

  A list of character strings defining the features/mz values to subset
  against.

- assay:

  A character string identifying the Seurat assay which contains the
  count data being subset.

## Value

A subset Seurat object containing only m/z values that were specified

## Examples

``` r
# SubsetMZFeatures(SeuratObj, c("mz-160","mz-170","mz-180"))
```
