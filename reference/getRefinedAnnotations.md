# Refines and reduces m/z annotations

Used to subset dataset to only include annotations that have n number of
entries. For example some peaks can have multiple annotations. Peaks
which have above n number of annotations assigned will be removed from
Seurat Object.

## Usage

``` r
getRefinedAnnotations(obj, assay = "Spatial", n = 1)
```

## Arguments

- obj:

  Seurat object needing annotation refinement. This object must have
  annotations present in 'obj`[[assay]]@meta.data`'

- assay:

  Character string defining the Seurat object assay where the annotation
  data is stored (default = "Spatial").

- n:

  Integer defining the number of entries an annotation can have
  assigned. Any higher counts will be removed (default = 1).

## Value

Refined Seurat object that only contains annotated mz values that have n
number of annotations assigned (per mz value)

## Examples

``` r
# HMDB_db <- load("data/HMDB_1_names.rds")
# AnnotatedSeuratObj <- AnnotateSeuratMALDI(SeuratObj, HMDB_db)

# getRefinedAnnotations(AnnotatedSeuratObj, n = 2)
```
