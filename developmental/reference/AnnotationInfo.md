# Inspect the metabolite annotation provenance stored in a SpaMTP object

Reports whether downstream pathway functions will use the current
indexed, scored RaMP annotation pipeline or a legacy mass-match result.

## Usage

``` r
AnnotationInfo(SpaMTP)
```

## Arguments

- SpaMTP:

  A SpaMTP Seurat object.

## Value

A named list containing annotation schema, engine, RaMP version,
provenance, and candidate count where available.
