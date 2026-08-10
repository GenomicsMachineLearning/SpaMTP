# Sums the intensity values of multiple m/z values into one

This function inputs any vector of m/z values and combines their
intensity values per pixel. This 'binned' data is stored in the SpaMTP
Object's `@meta.data` slot.

## Usage

``` r
BinMetabolites(
  data,
  mzs,
  assay = "Spatial",
  slot = "data",
  bin_name = "Binned_Metabolites"
)
```

## Arguments

- data:

  SpaMTP Seurat class object containing m/z intensities.

- mzs:

  Vector of m/z names to be binned together

- assay:

  Character string indicating which Seurat object assay to pull data
  form (default = "Spatial").

- slot:

  Character string indicating the assay slot to use to pull expression
  values form (default = "counts").

- bin_name:

  Character string defining the name of the meta.data column that stores
  the data (default = "Binned_Metabolites").

## Value

Binned intensity value stored in barcode meta.data slot

## Examples

``` r
# SpaMPT.obj <- BinMetabolites(SpaMPT.obj, mz = c('mz-740.471557617188','mz-784.528564453125','mz-897.603637695312'), bin_name = "Lipids")
```
