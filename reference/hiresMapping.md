# Maps SM pixels to high resolution ST data

Function used by MapSpatialOmics to align SM data to ST spots with
higher resolution (i.e. Xenium cells)

## Usage

``` r
hiresMapping(
  SM.data,
  ST.data,
  SM.assay = "Spatial",
  ST.assay = "Xenium",
  SM.fov = "fov",
  ST.image = "fov",
  SM.pixel.width = NULL,
  annotations = TRUE,
  add.metadata = TRUE,
  map.data = FALSE,
  new_SPT.assay = "SPT",
  new_SPM.assay = "SPM",
  verbose = FALSE
)
```

## Arguments

- SM.data:

  A SpaMTP Seurat object representing the Spatial Metabolomics data.

- ST.data:

  A Seurat object representing the Xenium Spatial Transcriptomics data.

- SM.assay:

  Character string defining the Seurat assay that contains the annotated
  counts and metadata corresponding to the m/z values (default =
  "Spatial").

- ST.assay:

  Character string specifying the current assay to use to extract
  transcriptional data from (default = "Spatial").

- SM.fov:

  Character string of the image fov associated with the spatial
  metabolomic data (default = "fov").

- ST.image:

  Character string matching the image name associated with the ST data
  stored in the Xenium object (default = "slice1").

- SM.pixel.width:

  Numeric value defining the width of each SM pixel. If set to `NULL`,
  the median pixel width will be calculated based on the distance
  between each pixel (default = NULL).

- annotations:

  Boolean value indicating if the Spatial Metabolomics (MALDI) Seurat
  object contains annotations assigned to m/z values (default = TRUE).

- add.metadata:

  Boolean defining whether to add the current metadata stored in the SM
  object to the new mapped multi-omic SpaMTP object (default = TRUE)

- map.data:

  Boolean indicating whether to map normalised/additional data stored in
  the `data` slot of the SpaMTP assay. Note: this process is
  computationally expensive with large datasets (default = FALSE).

- new_SPT.assay:

  Character string defining the assay name of the new overlaid SpaMTP
  Seurat object containing all updated transcriptomics data (default =
  "SPT").

- new_SPM.assay:

  Character string defining the assay name of the new overlaid SpaMTP
  Seurat object containing all updated metabolomic data (default =
  "SPM").

- verbose:

  Boolean value indicating whether to print informative progression
  update messages and progress bars (default = TRUE).

## Value

A SpaMTP Seurat object with the Spatial Metabolomic data mapped to
equivalent Spatial Transcripomics (Xenium) cells.

## Examples

``` r
# Helper function for MapSpatialOmics
```
