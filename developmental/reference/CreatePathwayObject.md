# Create a SpaMTP Seurat Object containing expression values for all present pathways

This function computes pathway-level scores from analyte-level
expression data and stores the results as a new assay in a Seurat
object. Each pathway score is calculated as the scaled mean expression
of the analytes associated with that pathway, adjusted by the square
root of the pathway size.

## Usage

``` r
CreatePathwayObject(
  object,
  assay = SeuratObject::DefaultAssay(object),
  slot = "scale.data",
  new.assay = "pathway",
  remove.nans = TRUE
)
```

## Arguments

- object:

  A SpaMTP Seurat object containing the expression data.

- assay:

  Character. Name of the assay to extract analyte expression from. If no
  value is assigned, the DefaultAssay of the SpaMTP Seurat Object will
  be used (default = `DefaultAssay(object)`).

- slot:

  Character. Which data slot to use (e.g., "scale.data") (default =
  "scale.data").

- new.assay:

  Character. Name of the new assay where pathway scores will be stored
  (defaults = "pathway").

- remove.nans:

  Logical. Whether to remove pathways with all NaN values (e.g., no
  matched analytes) (defaults = TRUE).

## Value

A SpaMTP Seurat object with a new assay containing pathway-level
expression scores. Feature names are adjusted to use underscores instead
of dashes.

## Examples

``` r
#object <- CreatePathwayObject(seurat_obj, assay = "RNA", slot = "scale.data")
```
