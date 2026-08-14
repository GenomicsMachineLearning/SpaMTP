# Spectral binning of intensity values stored in a Matrix object, converted from matter.

Spectral binning of intensity values stored in a Matrix object,
converted from matter.

## Usage

``` r
spectral_binning(
  matrix,
  ref,
  index,
  method = c("sum", "mean", "max", "min"),
  tolerance
)
```

## Arguments

- matrix:

  matter matrix object containing the intensity values to be binned.

- ref:

  A vector of reference mass-to-charge (m/z) values for binning. If left
  unspecified, mass range will be used to generate reference peaks.

- index:

  Character string specifying the column name for the m/z values
  (default = "mz").

- method:

  A character vector specifying the binning method. Options include
  `"sum"`, `"mean"`, `"max"`, `"min"`. If not specified default method
  used is "sum".

- tolerance:

  Numeric value specifying the tolerance for m/z matching (default =
  `NA`).

## Value

Matrix object containing the binned intensity values matching the
provided reference list.

## Examples

``` r
#Helper function for binning data in Matrix format
```
