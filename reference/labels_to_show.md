# Filters the annotation list to only include the first n number of annotations per m/z

Filters the annotation list to only include the first n number of
annotations per m/z

## Usage

``` r
labels_to_show(annotation_column, n = 3)
```

## Arguments

- annotation_column:

  Vector of the meta.data column containing the m/z annotations.

- n:

  Numeric value defining the number of annotations to keep (default =
  3).

## Value

Vector containing the first n number of annotations

## Examples

``` r
# labels_to_show(`SeuratObject[["Spatial"]]@meta.data$annotations`, n = 3)
```
