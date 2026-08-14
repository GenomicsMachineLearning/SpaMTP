# List MALDI matrix and on-tissue derivatization profiles

The registry separates conventional MALDI matrices, reactive matrices,
and derivatization reagents. A registered profile does not necessarily
imply that a universal mass-shift rule is available: some reactions
require a coupling reagent, produce several products, or depend strongly
on analyte structure and acquisition conditions.

## Usage

``` r
MALDIMatrixProfiles(matrix = NULL)
```

## Arguments

- matrix:

  Optional matrix/reagent name or alias. When `NULL`, return the
  complete registry.

## Value

A data frame describing supported matrix profiles and their current
automatic-rule status.
