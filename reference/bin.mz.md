# Bins multiple m/z values into one.

This function is a helper function for `BinMetabolites`.

## Usage

``` r
bin.mz(
  data,
  mz_list,
  assay = "Spatial",
  slot = "counts",
  stored.in.metadata = FALSE
)
```

## Arguments

- data:

  Seurat Spatial Metabolomic object containing mz values

- mz_list:

  Vector of m/z names to be binned into one value

- assay:

  Character string indicating which Seurat object assay to pull data
  form (default = "Spatial").

- slot:

  Character string indicating the assay slot to use to pull expression
  values form (default = "counts").

- stored.in.metadata:

  Boolean value indicating if the mz_list should be searched in the
  metadata DataFrame. If FALSE searches in Seurat object assay slot, if
  TRUE searches in metadata slot (default = FALSE).

## Value

Vector of binned mz counts for each spot/pixel

## Examples

``` r
# mz_values <- plusminus(SeuratObj, 448.2, 0.05)
# bin.mz(SeuratObj, mz_values)
```
