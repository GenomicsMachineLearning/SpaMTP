# Regional Pathway Enrichment

This the function used to compute the gene/metabolites set enrichment
for multi-omics spatial data

## Usage

``` r
FindRegionalPathways(
  SpaMTP,
  ident,
  DE.list,
  analyte_types = c("genes", "metabolites"),
  SM_assay = "SPM",
  ST_assay = "SPT",
  SM_slot = "counts",
  ST_slot = "counts",
  min_path_size = 5,
  max_path_size = 500,
  pval_cutoff_mets = 0.05,
  pval_cutoff_genes = 0.05,
  annotation_score_threshold = 0.05,
  annotation_source = c("current", "auto", "legacy"),
  verbose = TRUE
)
```

## Arguments

- SpaMTP:

  A SpaMTP Seurat object contains spatial
  metabolomics(SM)/transcriptomics(ST) data or both, if contains SM
  data, it should be annotated via SpaMTP::AnnotateSM function.

- ident:

  A name character to specific the cluster vector for regions in
  `SpaMTP@meta.data` slot.

- DE.list:

  A list consisting of differential expression data.frames for each
  input modality. Within each data.frame column names MUST include
  'cluster', 'gene', ('avg_log2FC' or 'logFC') and ('p_val_adj' or
  'FDR').

- analyte_types:

  Vector of character strings defining which analyte types to use.
  Options can be c("genes"), c("metabolites") or both (default =
  c("genes", "metabolites")).

- SM_assay:

  A Character string defining describing slot name for spatial
  metabolomics data in SpaMTP to extract intensity values from (default
  = "SPM").

- ST_assay:

  A Character string defining describing slot name for spatial
  transcriptomics data in SpaMTP to extract RNA count values from
  (default = "SPT").

- SM_slot:

  The slot name containing the SM assay matrix data (default =
  "counts").

- ST_slot:

  The slot name containing the ST assay matrix data (default =
  "counts").

- min_path_size:

  The min number of metabolites in a specific pathway (default = 5).

- max_path_size:

  The max number of metabolites in a specific pathway (default = 500).

- pval_cutoff_mets:

  Adjusted p-value cutoff used when constructing metabolite ranks. Set
  to `1` to include annotated metabolites regardless of DE significance
  (default = 0.05).

- pval_cutoff_genes:

  A numerical value defining the adjusted p value cutoff for significant
  differentially expressed genes. If `NULL` cutoff = `0.05` (default =
  0.05).

- annotation_score_threshold:

  Minimum indexed annotation score used for pathway mapping. It can be
  changed without re-running
  [`AnnotateSM()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/AnnotateSM.md)
  when the object was annotated with `min_score = 0` (default = 0.05).

- annotation_source:

  Metabolite annotation provenance. The default, `"current"`, requires
  scored RaMP IDs from the indexed annotation pipeline. `"auto"` permits
  a warned fallback to legacy `@tools$db_3`, while `"legacy"` explicitly
  requests that compatibility path.

- verbose:

  Boolean indicating whether to show informative messages. If FALSE
  these messages will be suppressed (default = TRUE).

## Value

A SpaMTP object with set enrichment on given analyte types.

## Examples

``` r
# SpaMTP = FindRegionalPathways(SpaMTP, polarity = "positive")
```
