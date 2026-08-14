# Create a Pathway Assay from Gene or Metabolite Data

This function creates a new assay within the provided SpaMTP Seurat
object which contains features (either genes or metabolites) labeled by
their respective RAMP ID. This assay can be used for running feature set
co-regulation analysis (based on
[GSCA](https://doi.org/10.1093/bioinformatics/btp502)).

## Usage

``` r
CreatePathwayAssay(
  SpaMTP,
  analyte_type = "metabolites",
  assay = "Spatial",
  slot = "counts",
  new_assay = "pathway",
  annotation_score_threshold = 0.05,
  annotation_source = c("current", "auto", "legacy"),
  verbose = TRUE
)
```

## Arguments

- SpaMTP:

  A SpaMTP Seurat object containing either spatial metabolic or
  transcriptomic data

- analyte_type:

  Character string specifying the type of analytes to process.Must be
  either "genes" or "metabolites" (default = "metabolites").

- assay:

  Character string specifying the name of the assay to use as source
  data (default = "Spatial").

- slot:

  Character string specifying which slot in the assay to use as source
  data (default = "counts").

- new_assay:

  Character string specifying the name of the new assay to create
  (default = "pathway").

- annotation_score_threshold:

  Minimum indexed annotation score used to map m/z features to RaMP
  compounds (default = 0.05).

- annotation_source:

  Metabolite annotation provenance. `"current"` requires the indexed,
  scored RaMP output; `"auto"` and `"legacy"` enable compatibility with
  older serialized SpaMTP objects.

- verbose:

  Boolean logical value indicating whether to print verbose messages
  during execution. (default = TRUE).

## Value

A SpaMTP object with a new assay added, containing respective
gene/metabolite data formatted based on RAMP_db IDs.

## Examples

``` r
## Create a pathway assay from metabolite data
#spamtp_obj <- CreatePathwayAssay(spamtp_obj, analyte_type = "metabolites", assay = "SPM", new_assay = "pathway")

## Create a pathway assay from gene data with verbose output
#spamtp_obj <- CreatePathwayAssay(spamtp_obj, analyte_type = "genes", assay = "SPT", new_assay = "gene_pathway", verbose = TRUE)
```
