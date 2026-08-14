# Load METASPACE Data and Create Seurat Object

Downloads METASPACE data, converts it to a feature matrix, and packages
it into a Seurat object with spatial metadata.

## Usage

``` r
Load_METASPACE(
  dataset_id,
  fdr = 0.1,
  database = c("HMDB", "v4"),
  api_key = NULL,
  relative = FALSE,
  transform = FALSE,
  verbose = TRUE
)
```

## Arguments

- dataset_id:

  Character string of the METASPACE dataset ID.

- fdr:

  Numeric, the FDR level threshold (default = 0.1).

- database:

  Character vector, specifying the database name and version (default =
  c("HMDB", "v4")).

- api_key:

  Optional character string, your METASPACE API key (default = NULL).

- relative:

  Logical, if TRUE, normalizes pixel intensities (0-1 range) (default =
  FALSE).

- transform:

  Boolean specifiying whether to transform/flip the image (default =
  FALSE).

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

## Value

A Seurat object initialized with the METASPACE data in the 'Spatial'
assay.

## Examples

``` r
# Load_METASPACE(dataset_id = "2020-12-07_03h16m14s", fdr = 0.1, database= c("HMDB", "v4"), relative = TRUE)
```
