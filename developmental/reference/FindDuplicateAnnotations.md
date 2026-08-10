# Finds if any metabolite is duplicated across multiple m/z values.

Finds if any metabolite is duplicated across multiple m/z values.

## Usage

``` r
FindDuplicateAnnotations(data, assay = "Spatial")
```

## Arguments

- data:

  Seurat Spatial Metabolomic Object containing annotated m/z values.

- assay:

  Character string defining the Seurat assay that contains the annotated
  metadata corresponding to the m/z values (default = "Spatial").

## Value

Vector of character strings describing metabolites that are assigned to
multiple m/z values

## Examples

``` r
# FindDuplicateAnnotations(SeuratObj)
```
