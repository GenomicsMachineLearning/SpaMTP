# Gets values from a single metadata column for a respective m/z value.

Gets values from a single metadata column for a respective m/z value.

## Usage

``` r
GetMZMetadata(
  obj,
  mz,
  assay = "Spatial",
  metadata.column = "all_IsomerNames",
  separate = TRUE
)
```

## Arguments

- obj:

  SpaMTP Spatial Metabolomic Seurat Object containing annotated m/z
  values.

- mz:

  Character string specifying the m/z value to return the metadata for.

- assay:

  Character string defining the Seurat assay that contains the annotated
  metadata corresponding to the m/z values (default = "Spatial").

- metadata.column:

  Character string corresponding to the `@meta.data` column to extract
  the data from (default = "all_IsomerNames").

- separate:

  Boolean indicating whether to separate the metadata string and return
  a vector. Note, if `TRUE` the metadata column should contain values
  separate by "; " (default = TRUE).

## Value

Vector of character strings containing the metadata for a specific mz
value.

## Examples

``` r
# GetMZMetadata(SpaMTP, mz = "mz-100", metadata.column = "all_IsomerNames")

##### Example for getting metabolite annotation IDs for a m/z value
# GetMZMetadata(SpaMTP, mz = "mz-100", metadata.column = "all_Isomers")
```
