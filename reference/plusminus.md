# Identifies all mz peaks within a plus-minus range of the target_mz

Identifies all mz peaks within a plus-minus range of the target_mz

## Usage

``` r
plusminus(data, target_mz, plus_minus)
```

## Arguments

- data:

  Seurat Spatial Metabolomic object containing mz values

- target_mz:

  Numeric value defining the target m/z peak

- plus_minus:

  Numeric value defining the range/threshold either side of the target
  peak to also be identified

## Value

Vector of m/z values within plus-minus range of the target mz value

## Examples

``` r
# plusminus(SeuratObj, target_mz = 400.01, plus_minus = 0.005)
```
