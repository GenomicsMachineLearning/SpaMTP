# Get top spatially variable metabolites

This function returns the names of the top n number of spatially
variable features.

## Usage

``` r
GetSpatiallyVariableMetabolites(object, assay = "SPM", n = 10)
```

## Arguments

- object:

  SpaMTP Seurat class object containing the intensity values for each
  m/z

- assay:

  Character string indicating which Seurat object assay to pull data
  form (default = "SPM").

- n:

  Numeric value defining the top number of metabolites to return
  (default = 10).

## Value

Vector containing m/z feature names corresponding to the top n number of
spatially variable metabolites

## Examples

``` r
# features <- GetSpatiallyVariableMetabolites(SpaMTP, n = 6)
```
