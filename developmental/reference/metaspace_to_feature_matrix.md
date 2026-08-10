# Convert METASPACE Images to Feature Matrix

Transforms the list of 2D ion image matrices from `get_metaspace` into a
single matrix suitable for spatial analysis packages.

## Usage

``` r
metaspace_to_feature_matrix(metaspace_data, transform = FALSE, verbose = TRUE)
```

## Arguments

- metaspace_data:

  A list object returned by `get_metaspace`.

- transform:

  Logical. If TRUE, transposes each image matrix before flattening
  (e.g., to align coordinate systems) (default = FALSE).

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

## Value

A numeric matrix where rows are molecular features (m/z) and columns are
spatial pixels (named 'x_y').

## Examples

``` r
# metaspace_to_feature_matrix(mtx, transform = TRUE)
```
