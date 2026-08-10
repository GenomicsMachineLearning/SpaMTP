# Read Spatial Metabolomics matrix file (.csv format)

This function reads in SM data stored in a table format. This table must
contain two rows named 'x' and 'y' storing the respective pixel
coordinates. All other columns should be the relative m/z values
containing the intensity values of each pixel.

## Usage

``` r
ReadSM_mtx(
  mtx.file,
  assay = "Spatial",
  verbose = TRUE,
  feature.start.column = 1,
  mz.prefix = NULL,
  project.name = "SpaMTP"
)
```

## Arguments

- mtx.file:

  Character string defining the path of the spatial metabolomic image
  matrix .csv file.

- assay:

  Character string of the Seurat object assay name to store the relative
  intensity data (default = "Spatial").

- verbose:

  Boolean indicating whether to show the message. If TRUE, the message
  will be shown; else, it will be suppressed (default = TRUE).

- feature.start.column:

  Numeric value defining the start index containing the x, y, and m/z
  value columns within the table (default = 1).

- mz.prefix:

  Character string matching the prefix string in front of each m/z name
  (default = NULL).

- project.name:

  Character string defining the name of the sample to be assigned as
  `orig.ident` (default = "SpaMTP").

## Value

A SpaMTP Seurat class object containing the intensity values in the
counts slot of the designated assay.

## Details

**NOTE:** The input file must be in a format similar to the table below:

    A data.frame: 5 × 5
       x   y   mz1  mz2  mz3
    1  0   1    0    0   11
    2  0   2    0    0    0
    3  0   3    0    0    0
    4  0   4   20    0    0
    5  0   5    0    0    0

- The first two columns (`x`, `y`) contain the respective spatial
  coordinates.

- The subsequent columns contain the m/z values and their intensities
  for each spatial pixel.

## Examples

``` r
# msi_data <- ReadSM_mtx("~/Documents/msi_mtx.csv")
```
