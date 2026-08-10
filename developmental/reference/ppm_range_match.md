# Calculates the ppm range and check if mz values are within the range

- Returns TRUE if match is found and false if no match.

Calculates the ppm range and check if mz values are within the range

- Returns TRUE if match is found and false if no match.

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
