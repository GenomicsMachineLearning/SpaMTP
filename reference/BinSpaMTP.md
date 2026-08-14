# Bin SpaMTP Object

Bins m/z peaks of a specified SpaMTP object to reduce dimensionality and
remove noise.

## Usage

``` r
BinSpaMTP(
  data,
  resolution,
  units = "ppm",
  assay = "Spatial",
  slot = "counts",
  method = c("sum"),
  return.only.mtx = FALSE
)
```

## Arguments

- data:

  SpaMTP Seurat Object containing the intensity data to be binned

- resolution:

  Numeric value indicating the relative bin size to use when binning the
  intensity data.

- units:

  Character string defining the relative units of the provided
  resolution size. Values can be either 'ppm' or 'mz' (default = 'ppm').

- assay:

  Character string defining the SpaMTP assay to extract intensity values
  from (default = "Spatial").

- slot:

  Character string defining the assay slot containing the intensity
  values (default = "counts").

- method:

  A character vector specifying the binning method. Input values must be
  either `"sum"`, `"mean"`, `"max"` or `"min"` (default = "sum").

- return.only.mtx:

  Boolean value indicating whether to return only the binned intensity
  matrix. If `FALSE`, a SpaMTP object will be returned with the binned
  values stored in an assay called `binned` (default = FALSE).

## Value

Either a binned intensity matrix or a SpaMTP object contatining the
binned intensity matrix stored in the `binned` assay.

## Examples

``` r
#BinSpaMTP(spamtp.obj, resolution = 10, units = "ppm", return.only.mtx = TRUE)
```
