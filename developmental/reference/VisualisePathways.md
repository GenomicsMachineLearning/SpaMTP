# Visualise Significant Pathways

Displays the pathway analysis results form running the
'FishersPathwayAnalysis()' function

## Usage

``` r
VisualisePathways(
  SpaMTP,
  pathway_df,
  assay = "SPM",
  slot = "counts",
  min_n = 3,
  p_val_threshold = 0.1,
  method = "ward.D2",
  verbose = TRUE,
  ...
)
```

## Arguments

- SpaMTP:

  SpaMTP Seurat object used to run FishersPathwayAnalysis function.

- pathway_df:

  Dataframe containing the pathway enrichment results (output from
  SpaMTP::FishersPathwayAnalysis function).

- assay:

  Character string defining the SpaMTP assay that contains m/z values
  (default = "SPM").

- slot:

  Character string defining the assay slot contain the intensity values
  (default = "counts").

- min_n:

  Integer value specifying the minimum number of analytes required to be
  present in a pathway (default = 3).

- p_val_threshold:

  The p-val cutoff to keep the pathways generated from fisher exact test
  (default = "0.1").

- method:

  Character string defining the statistical method used to calculate
  hclust (default = "ward.D2").

- verbose:

  Boolean indicating whether to show informative messages. If FALSE
  these messages will be suppressed (default = TRUE).

- ...:

  The arguments pass to stats::hclust

## Value

A combined gg, ggplot object with pathway and dendrogram

## Examples

``` r
#SpaMTP:::VisualisePathways(SpaMTP =seurat,pathway_df = pathway_df,p_val_threshold = 0.1,assay = "Spatial",slot = "counts")
```
