# Create a singular multiomics assay by merging data from multiple assays.

Combines the `scale.data` slots from multiple assays in a SpaMTP Seurat
object into a single new assay. Useful for integrating multiple
modalities (e.g. transcriptomics, proteomics, metabolomics) that have
already been scaled.

## Usage

``` r
CreateMergedModalityAssay(
  SpaMTP,
  assays.to.merge,
  new.assay = "merged",
  return.original = TRUE,
  verbose = FALSE
)
```

## Arguments

- SpaMTP:

  A SpaMTP Seurat object that contains atleast two assays to be merged.

- assays.to.merge:

  A character vector specifying the names of assays whose `scale.data`
  slots should be merged. At least two assay names must be provided and
  both must contain the `scale.data` slot, for example: assays.to.merge
  = c("SPM", "SPT").

- new.assay:

  A character string specifying the name of the new assay to be created
  (default = "merged").

- return.original:

  Boolean value defining if the returned SpaMTP Seurat object will
  contain the original individual assays. If the data size is large it
  is recommended to set to `False` (default = TRUE).

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

## Value

A SpaMTP Seurat object containing a new assay with the merged scaled
data values.

## Details

This function assumes that each specified assay has been processed with
[`Seurat::ScaleData()`](https://satijalab.org/seurat/reference/ScaleData.html)
(or an alternative scaling method), and that their `scale.data` slots
contain numeric matrices. The merged assay will use the row-bound
`scale.data` matrices as the `counts`, `data`, and `scale.data` slots.

## Examples

``` r
# merged_obj <- CreateMergedModalityAssay(SpaMTP = spamtp_obj, assays.to.merge = c("SPM", "SPT"),new.assay = "merged")
```
