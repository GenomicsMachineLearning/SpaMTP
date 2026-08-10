# Assign custom annotations to m/z values

Adds custom metabolite annotations to respective m/z values (ideal for
specific matrices such as FMP10).

## Usage

``` r
AddCustomMZAnnotations(
  data,
  annotations,
  assay = "Spatial",
  return.only.annotated = FALSE,
  mass.threshold = 0.05,
  annotation.column = "all_IsomerNames"
)
```

## Arguments

- data:

  SpaMTP Seurat object containing m/z intensity values.

- annotations:

  data.frame containing two columns named 'annotation' and 'mass'. These
  columns should contain the custom metabolite annotation and the
  relative m/z mass respectively.

- assay:

  Character string defining the Seurat object assay to store the
  respective annotations in the feature meta.data dataframe (default =
  "Spatial").

- return.only.annotated:

  Boolean defining whether to return a SpaMTP Seurat object containing
  only successfully annotated m/z values (default = FALSE).

- mass.threshold:

  Numeric value defining the acceptable threshold (plus-minus) between
  the custom annotations and the actual m/z values contained within the
  SpaMTP object (default = 0.05).

- annotation.column:

  Character string defining the feature meta.data column name that will
  contain the assigned annotations (default = "all_IsomerNames").

## Value

SpaMTP Seurat object containing the custom annotations stored in the
feature metadata dataframe.

## Examples

``` r
# annotated_data <- AddCustomMZAnnotations(SpaMTP.obj, annotation.df)
```
