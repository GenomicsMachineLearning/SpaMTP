# Saves SpaMTP Object

This function saves a SpaMTP Seurat Object into a standard
single-cell/spatial file format. This includes a
filtered_feature_bc_matrix folder containing files storing the features,
barcode/pixels and intensity matrix. Metadata and sapatial files (such
as scale factors and hires/lowres images) are also stored.

## Usage

``` r
SaveSpaMTPData(
  data,
  outdir,
  assay = "Spatial",
  slot = "counts",
  image = NULL,
  annotations = FALSE,
  generate.h5 = TRUE,
  verbose = TRUE
)
```

## Arguments

- data:

  A Spatial Metabolomic SpaMTP Seurat Object being saved.

- outdir:

  Character string of the directory to save the mtx.mtx, barcode.tsv,
  features.tsv, barcode_metadata.csv and feature_metadata.csv in.

- assay:

  Character string defining the Seurat assay that contains the m/z count
  data (default = "Spatial").

- slot:

  Character string defining the Seurat assay slot that contains the m/z
  values directly (default = "counts").

- image:

  Character string defining the image stored within the SpaMTP Seurat
  object to save. If `NULL` no image will be saved (default = NULL).

- annotations:

  Boolean values defining if the Seurat Object contains annotations to
  be saved (default = FALSE).

- generate.h5:

  Boolean value indicating whether to generate a
  filtered_feature_bc_matrix.h5 file. Often used by data loading
  functions (e.g. scanpy.load_visium). If `FALSE`, only a
  filtered_feature_bc_matrix folder will be generated (default = TRUE).

- verbose:

  Boolean indicating whether to show informative processing messages. If
  TRUE the message will be show, else the message will be suppressed
  (default = TRUE).

  ### Details

  - This can be used for saving data for transfer to python. Can be read
    in as Anndata using scanpy.read_10x_mtx().

  - For saving in R saveRDS() is recommended.

## Examples

``` r
# saveSpaMTPData(SeuratObject, "../output", annotations = TRUE)
```
