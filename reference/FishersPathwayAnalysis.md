# Calculates Significant Metabolic Pathways using a Fisher Exact Test

Calculates Significant Metabolic Pathways using a Fisher Exact Test

## Usage

``` r
FishersPathwayAnalysis(
  Analyte,
  max_path_size = 500,
  min_path_size = 5,
  alternative = "greater",
  pathway_all_info = FALSE,
  pval_cutoff = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- Analyte:

  A list of analytes containing a combination of three possible
  elements, namely "mzs", "genes" and/or "metabolites". The list must be
  named with these titles, corresponding to the relative input datasets.
  Read below for supported input formats.

- max_path_size:

  The max number of in a specific pathway (default = 500).

- min_path_size:

  The min number of in a specific pathway (default = 5).

- alternative:

  The hypothesis of the fisher exact test (default = "greater").

- pathway_all_info:

  Whether to included all genes/ screened in the return (default =
  FALSE).

- pval_cutoff:

  A numerical value defining the adjusted p value cutoff for keeing
  significant pathways (default = NULL).

- verbose:

  Boolean indicating whether to show informative messages. If FALSE
  these messages will be suppressed (default = TRUE).

- ...:

  Additional parameters that can be passed through to
  [`annotateTable()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/annotateTable.md)
  when running `mz`-based analysis. Please see documentation for
  [`annotateTable()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/annotateTable.md)
  for more details.

  ### Details

  - Supported `metabolites` format: strings which contain the metabolite
    ID with database name. For example = "hmdb:HMDBX", "chebi:X",
    "pubchem:X","wikidata:X" ,"kegg:X"
    ,"CAS:X","lipidbank:X","chemspider:X"," LIPIDMAPS:X" (where X stands
    for upper case of the cooresponding ID in each database)

  - Supported `genes` data format: strings which contain the gene name
    and formatting. For example = "entrez:X", "gene_symbol:X",
    "uniprot:X", "ensembl:X", "hmdb:HMDBPX"

  - Supported `mzs` format: any string or numeric vector containing m/z.
    If `mzs` values are provided, the current indexed annotation
    pipeline and the bundled RaMP `chem_props` table are used by
    default.

## Value

a dataframe with the relevant pathway information

## Examples

``` r
## Running in 'mzs' mode:
# FishersPathwayAnalysis(Analyte = list("mzs" = mz_values), ppm_error = 3)

## Running in 'metabolites' mode
# FishersPathwayAnalysis(Analyte = list("metabolites" = metabolite_ids))

## Running 'metabolites' and 'genes' combined
# FishersPathwayAnalysis(Analyte = list("metabolites" = metabolite_ids, "genes" = gene_names))
```
