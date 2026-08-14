# Launch an Interactive ROI Annotation App for Seurat Spatial Data

This function opens a Shiny app to allow users to manually annotate
Regions of Interest (ROIs) on spatial transcriptomics or metabolomics
data stored in a Seurat object. Users can interactively select regions
using lasso selection, assign custom names to each ROI, and save the
results as binary (0/1) metadata columns in the Seurat object.

## Usage

``` r
SelectROIs(seurat_obj, image = "fov")
```

## Arguments

- seurat_obj:

  A `Seurat` object containing spatial coordinates and metadata.

- image:

  A character string specifying the spatial image source used by
  `GetTissueCoordinates()`. Default is `"fov"`. Set to `"VisiumV2"` for
  Visium-style spatial data.

## Value

A modified `Seurat` object with added metadata columns for each saved
ROI (1 = inside ROI, 0 = outside).

## Details

The app includes options to:

- Choose a metadata column to display.

- Adjust the spot size for display.

- Use lasso to select spots and name each ROI.

- Save each ROI to the Seurat object metadata.

- Export the final Seurat object upon clicking "Finish".

Spatial selection uses `sf` geometry tools. The plot is rendered using
`plotly`.
