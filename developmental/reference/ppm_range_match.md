# Check whether m/z values match within a ppm tolerance

Calculates the absolute ppm difference between observed and reference
m/z values and returns whether each comparison satisfies the tolerance.

## Usage

``` r
ppm_range_match(observed_mz, reference_mz, ppm)
```

## Arguments

- observed_mz:

  Numeric value defining the observed mz value.

- reference_mz:

  Numeric value defining the reference mz value.

- ppm:

  Numeric value defining the maximum acceptable ppm_error/threshold for
  searching.

## Value

Boolean value indicating if a match is found (TRUE) or not (FALSE)

## Examples

``` r
### Helper Function ###
```
