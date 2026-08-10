# Core METASPACE Data Retrieval

The main user-facing function to download annotations and ion images for
a METASPACE dataset.

## Usage

``` r
get_metaspace(
  dataset_id,
  fdr = 0.1,
  database = c("HMDB", "v4"),
  include_images = TRUE,
  api_key = NULL,
  isotope_idx = 1,
  relative = TRUE
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

- include_images:

  Logical, if TRUE, downloads all ion image matrices (default = TRUE).

- api_key:

  Optional character string, your METASPACE API key (default = NULL).

- isotope_idx:

  Integer, index of the isotope image to download for all annotations (1
  for main peak) (default = 1).

- relative:

  Logical, if TRUE, normalizes pixel intensities (0-1 range) (default =
  TRUE).

## Value

A list with two elements: `annotations` (data.frame) and `images` (list
of 2D matrices).

## Examples

``` r
# get_metaspace(dataset_id = "2020-12-07_03h16m14s", fdr = 0.1, database= c("HMDB", "v4"), relative = TRUE)
```
