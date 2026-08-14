# Add GESECA Scores to SpaMTP Object

Adds Feature Set Enrichment Score for a specific pathway to a SpaMTP
Seurat object's metadata. This is a helper function for
`PlotPathwaysSpatially`.

## Usage

``` r
addGesecaScores(
  pathways,
  object,
  assay = SeuratObject::DefaultAssay(object),
  slot = "scale.data",
  prefix = "",
  scale = FALSE
)
```

## Arguments

- pathways:

  List of gene/metabolite sets where the list names (pathway name)
  become metadata column names.

- object:

  SpaMTP Seurat object to add the pathway scores to.

- assay:

  Character string defining the name of assay to use for data extraction
  (default = \`SeuratObject::DefaultAssay(object)).

- slot:

  Character string stating the name of slot to extract data from
  (default = "scale.data").

- prefix:

  Character string to append before each pathway name in metadata
  columns (default = "").

- scale:

  Boolean logical value indicating whether to scale the GESECA score
  (default = FALSE).

## Value

SpaMTP Seurat object with added GESECA scores in metadata.

## Examples

``` r
### HELPER FUNCTION
#glycolysis_list <- list("Glycolysis" = c("RAMP_C_000218730","RAMP_G_000012583","RAMP_G_000001171","RAMP_G_000007564","RAMP_C_000218226","RAMP_C_000040403","RAMP_C_000001115","RAMP_G_000008859"))
#spamtp_obj <- addGesecaScores(pathways = glycolysis_list, object = spamtp_obj ,assay = "RNA", scale = TRUE)
```
