# Maps Spatial Metabolomic (MALDI) data to corresponding Spatial Transcriptomics data and coordinates.

Maps Spatial Metabolomic (MALDI) data to corresponding Spatial
Transcriptomics data and coordinates.

## Usage

``` r
MapSpatialOmics(
  SM.data,
  ST.data,
  ST.hires = FALSE,
  SM.assay = "Spatial",
  ST.assay = "Spatial",
  SM.fov = "fov",
  ST.image = "slice1",
  ST.scale.factor = "hires",
  SM.pixel.width = NULL,
  overlap.threshold = 0.2,
  annotations = TRUE,
  add.metadata = TRUE,
  merge.unique.metadata = TRUE,
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

  A Seurat object representing the Spatial Transcriptomics data.

- ST.hires:

  Boolean string defining if the ST data is at a higher resolution
  compared to the SM pixel data. For example, generally Visium data will
  be lower res whereas Xenium/single-cell resolution spatial data will
  be a higher resolution (default = FALSE).

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
  such as 'fov' or 'slice1' object (default = "slice1").

- ST.scale.factor:

  Character string defining the image resolution associated with the
  Visium image pixel data. If `NULL` the full-res coordinates will be
  used and no scaling will be performed. Note: This parameter is only
  required for aligning lowres ST data (default = "hires").

- SM.pixel.width:

  Numeric value defining the width of each SM pixel. If set to `NULL`,
  the median pixel width will be calculated based on the distance
  between each pixel (default = NULL).

- overlap.threshold:

  Numeric value defining the overlap proportion threshold for a SM pixel
  to be associated with a ST spot. For example, if res_increase = 0.2
  then pixels that have at least 20% area overlap with the respective
  visium spot will be assigned a match. Note: This parameter is only
  required for aligning lowres ST data (default = 0.2).

- annotations:

  Boolean value indicating if the Spatial Metabolomics (MALDI) Seurat
  object contains annotations assigned to m/z values (default = TRUE).

- add.metadata:

  Boolean defining whether to add the current metadata stored in the SM
  object to the new mapped multi-omic SpaMTP object (default = TRUE)

- merge.unique.metadata:

  Boolean indicating whether to summaries duplicated metadata terms to
  only store unique values in the metadata. Note: This parameter is only
  required for aligning lowres ST data, and `add.metadata` must be set
  to `TRUE` for this functionality to be implemented (default = TRUE).

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

A SpaMTP Seurat object with Spatial Metabolomic data mapped to
equivalent Spatial Transcripomics coordinates (Visium spots/Xenium
cells).

## Examples

``` r

## Mapping MALDI data to equivalent Visium spots
# MapSpatialOmics(VisiumObj, SeuratObj, ST.scale.factor = "hires", SM.assay = "Spatial", ST.assay = "Spatial")

#' ## Mapping MALDI data to equivalent Xenium cells
# MapSpatialOmics(VisiumObj, SeuratObj, SM.assay = "Spatial", ST.assay = "Xenium")
```
