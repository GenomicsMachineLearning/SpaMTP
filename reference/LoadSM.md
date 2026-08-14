# Loads spatial metabolic data into a SpaMTP Seurat Object

This function loads raw spatial metabolic data in a .imzML and .ibd
format and generates a SpaMTP Seurat Object. This function adapts the
`readImzML` function implmented in Cardinal to correctly import large
data.

## Usage

``` r
LoadSM(
  file,
  mass.range = NULL,
  resolution = NA,
  units = "ppm",
  verbose = TRUE,
  assay = "Spatial",
  multi.run = FALSE
)
```

## Arguments

- file:

  Character string defining the directory path of the file and the file
  name.

- mass.range:

  Vector of numeric values indicating the mass range to use for the
  imported data (default = NULL).

- resolution:

  Numeric value defining the the accuracy to which the m/z values will
  be binned after reading. This value can be in either "ppm" or "mz"
  depending on the units type specified (default = 10).

- units:

  Character string defining the resolution value unit type, either
  c("ppm", "mz") (default = "ppm")

- verbose:

  Boolean indicating whether to show informative processing messages. If
  TRUE the message will be show, else the message will be suppressed
  (default = TRUE)

- assay:

  Character string describing the name of the new assay which stores the
  imported data (default = "Spatial").

- multi.run:

  Boolean indicating if there are multiple runs within the imported
  data. If `TRUE`, an index will be added to the pixel names per run,
  and an individual FOV will be generated per run in the Seurat Object
  (default = FALSE).

## Value

A new SpaMTP Seurat object contain the imported spatial metabolic
intensity values

## Examples

``` r
# data <-LoadSM(name = "run1", folder = "/Documents/SpaMTP_test_data/", mass.range = c(160,1500), resolution = 10, assay = "Spatial")
```
