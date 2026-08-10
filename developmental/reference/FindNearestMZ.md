# Finds the nearest m/z peak to a given value in a SpaMTP Object

Finds the nearest m/z peak to a given value in a SpaMTP Object

## Usage

``` r
FindNearestMZ(data, target_mz, assay = NULL)
```

## Arguments

- data:

  Seurat Spatial Metabolomic object containing mz values

- target_mz:

  Numeric value defining the target m/z peak

- assay:

  Character string indicating which Seurat object assay to pull data
  form. If set to NULL will use current default assay based on Seurat's
  `DefaultAssay()` function (default = NULL).

## Value

String of the closest m/z value within the given dataset

## Examples

``` r
# FindNearestMZ(SeuratObj, target_mz = 400.01)
```
