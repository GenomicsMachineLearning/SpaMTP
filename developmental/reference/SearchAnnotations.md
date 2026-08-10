# Find Annotation

Searches through annotated m/z values to return all which contain the
metabolite search term provided

## Usage

``` r
SearchAnnotations(
  data,
  metabolite,
  assay = "Spatial",
  search.exact = FALSE,
  column.name = "all_IsomerNames"
)
```

## Arguments

- data:

  Seurat Spatial Metabolomic Object containing annotated m/z values.

- metabolite:

  Character string of metabolite search term.

- assay:

  Character string defining the Seurat assay that contains the annotated
  metadata corresponding to the m/z values (default = "Spatial").

- search.exact:

  Boolean value defining if to only return m/z values which contain the
  exact match to the metabolite search term (default = FALSE).

- column.name:

  Character string defining the column name where the annotations are
  stored in the slot meta.data (default = "all_IsomerNames").

## Value

A Data.Frame containing the peak metadata corresponding to the
metabolite search term provided

## Examples

``` r
# SearchAnnotations(SeuratObj, "Glucose", search.exact = TRUE)
```
