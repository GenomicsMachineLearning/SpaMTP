# SpaMTP: Indexed Metabolite Annotation with RaMP 3.0

## Overview

SpaMTP’s current annotation pipeline separates candidate generation from
downstream filtering. It generates chemically valid adduct hypotheses
once, sorts their expected mass-to-charge ratios, queries that index
within a strict mass tolerance, and ranks the remaining candidates using
optional MS1 context.

The workflow is:

1.  Select a metabolite database and ion-mode-specific adduct rules.
2.  Build a reusable sorted annotation index.
3.  Query one or many observed m/z values within a ppm tolerance.
4.  Score candidates using mass accuracy, chemical heuristics, isotope
    evidence, and adduct-family support.
5.  Store all ppm-valid candidates and choose a score threshold at
    analysis time.

This article describes upstream candidate annotation. The
[`Annotation Refinement`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/articles/Annotation_Refinement.md)
article covers downstream methods for resolving candidates using
biological and spatial information.

``` r

library(SpaMTP)
```

## RaMP snapshot and identifiers

SpaMTP contains a deliberately pruned RaMP snapshot. The annotation
table keeps chemical structures, formulae, monoisotopic masses, names,
source identifiers, and current `RAMP_C_*` identifiers; pathway tables
retain only fields required by SpaMTP’s pathway and network functions.

``` r

ramp_summary <- data.frame(
  property = c(
    "RaMP version",
    "Chemical-property records",
    "Graph nodes after harmonisation",
    "Stale graph IDs remaining"
  ),
  value = c(
    ramp_db_metadata$ramp_version,
    format(ramp_db_metadata$snapshot_counts[["chem_props"]], big.mark = ","),
    format(
      ramp_db_metadata$graph_harmonisation$graph_nodes_after,
      big.mark = ","
    ),
    ramp_db_metadata$graph_harmonisation$stale_ids_remaining
  )
)
knitr::kable(ramp_summary, align = c("l", "r"))
```

| property                        |   value |
|:--------------------------------|--------:|
| RaMP version                    |   3.0.7 |
| Chemical-property records       | 289,754 |
| Graph nodes after harmonisation |  32,899 |
| Stale graph IDs remaining       |       0 |

The bundled chemical table can be used directly:

``` r

data.frame(column = names(chem_props)) |>
  knitr::kable(align = "l")
```

| column           |
|:-----------------|
| ramp_id          |
| chem_data_source |
| chem_source_id   |
| iso_smiles       |
| inchi_key_prefix |
| inchi_key        |
| inchi            |
| mw               |
| monoisotop_mass  |
| common_name      |
| mol_formula      |

## A small 2-hydroxyglutarate example

The examples below use the RaMP entry for D-2-hydroxyglutaric acid
(2HG), with neutral monoisotopic mass approximately 148.0372 Da. Two
carboxyl groups give this demonstration molecule a maximum of two
exchangeable protons.

``` r

target_db <- chem_props[
  chem_props$ramp_id == "RAMP_C_000219574" &
    grepl("HMDB0000606", chem_props$chem_source_id),
  c(
    "ramp_id", "chem_data_source", "chem_source_id", "inchi_key",
    "monoisotop_mass", "common_name", "mol_formula"
  ),
  drop = FALSE
]
stopifnot(nrow(target_db) > 0L)
target_db <- target_db[1L, , drop = FALSE]
target_db$max_exchangeable_protons <- 2L

knitr::kable(
  target_db[, c(
    "ramp_id", "chem_source_id", "common_name", "mol_formula",
    "monoisotop_mass", "max_exchangeable_protons"
  )],
  digits = 7
)
```

|  | ramp_id | chem_source_id | common_name | mol_formula | monoisotop_mass | max_exchangeable_protons |
|:---|:---|:---|:---|:---|---:|---:|
| 1441 | RAMP_C_000219574 | hmdb:HMDB0000606 | D-2-Hydroxyglutaric acid | C5H8O5 | 148.0372 | 2 |

`max_exchangeable_protons` is optional database metadata. When it is
present, SpaMTP removes impossible proton-loss states before indexing.
When it is absent, SpaMTP does not infer this property from molecular
formula alone.

## Validated adduct rules

[`AdductRules()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/AdductRules.md)
returns the ion masses, molecular multiplier, charge, proton-loss count,
chemical class, and prior used by the engine. Bracketed and unbracketed
forms such as `[M+H]+` and `M+H` are both accepted.

``` r

selected_adducts <- c("M+H", "M+Na", "M+2Na-H", "2M+H")
positive_rules <- AdductRules("positive")
rule_summary <- positive_rules[
  positive_rules$name %in% selected_adducts,
  c(
    "name", "n_molecules", "mass_shift", "charge", "loss_h",
    "complexity", "base_adduct", "prior"
  )
]
rule_summary <- rule_summary[match(selected_adducts, rule_summary$name), ]
knitr::kable(rule_summary, digits = 7, row.names = FALSE)
```

| name    | n_molecules | mass_shift | charge | loss_h | complexity     | base_adduct | prior |
|:--------|------------:|-----------:|-------:|-------:|:---------------|:------------|------:|
| M+H     |           1 |   1.007276 |      1 |      0 | simple         | M+H         |  1.00 |
| M+Na    |           1 |  22.989221 |      1 |      0 | salt           | M+Na        |  0.90 |
| M+2Na-H |           1 |  44.971165 |      1 |      1 | metal_exchange | M+Na        |  0.65 |
| 2M+H    |           2 |   1.007276 |      1 |      0 | multimer       | M+H         |  0.55 |

For each rule, SpaMTP calculates

``` math
m/z_{expected} =
\frac{nM + \Delta m_{adduct} - \Delta m_{loss}}{|z|}.
```

The stored `mass_shift` already contains the signed sum of additions and
losses. Custom rule tables must satisfy charge balance before candidate
generation.

## Build the reusable index

[`BuildMZAnnotationIndex()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/BuildMZAnnotationIndex.md)
performs the expensive expansion once. Expected ion masses are sorted,
allowing
[`QueryMZAnnotationIndex()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/QueryMZAnnotationIndex.md)
to use binary interval search rather than scanning every database record
for every observed peak.

``` r

positive_index <- BuildMZAnnotationIndex(
  target_db,
  polarity = "positive",
  adducts = selected_adducts
)
positive_index
```

    ## SpaMTP m/z annotation index
    ##   polarity: positive
    ##   metabolites: 1
    ##   adduct rules: 4
    ##   expected ions: 4

``` r

ion_table <- data.frame(
  adduct = positive_index$rules$name[positive_index$rule_idx],
  expected_mz = positive_index$expected_mz
)
knitr::kable(ion_table, digits = 9, row.names = FALSE)
```

| adduct  | expected_mz |
|:--------|------------:|
| M+H     |    149.0444 |
| M+Na    |    171.0264 |
| M+2Na-H |    193.0083 |
| 2M+H    |    297.0816 |

For large datasets, build one index for each database, polarity, and
adduct set, then reuse it for every spectrum or batch of observed peaks.

### Chemical validity pruning

The 2HG record permits at most two proton losses. Therefore `[M-3H]3-`
is not included in the indexed search space.

``` r

negative_index <- BuildMZAnnotationIndex(
  target_db,
  polarity = "negative",
  adducts = c("M-H", "M-2H", "M-3H")
)
indexed_negative_rules <- unique(
  negative_index$rules$name[negative_index$rule_idx]
)
data.frame(indexed_rule = indexed_negative_rules) |>
  knitr::kable(align = "l")
```

| indexed_rule |
|:-------------|
| M-2H         |
| M-H          |

## Query and rank candidates

The following observed peaks have small controlled mass offsets. The
contextual MS1 spectrum also contains a carbon-13 peak for `[M+H]+` at
+1.00335 Da. In real data, `ms1_spectrum` should represent the same
retention-time window, spatial pixel, or biologically coherent region as
the queried peaks.

``` r

ppm_offsets <- c(1.2, -1.4, 0.6, -0.8)
observed_mz <- ion_table$expected_mz * (1 + ppm_offsets * 1e-6)
names(observed_mz) <- ion_table$adduct

carbon13_spacing <- 1.00335483507
ms1_spectrum <- data.frame(
  mz = c(
    observed_mz,
    ion_table$expected_mz[ion_table$adduct == "M+H"] + carbon13_spacing
  ),
  intensity = c(1000, 700, 300, 200, 53.5)
)

candidates <- QueryMZAnnotationIndex(
  observed_mz,
  positive_index,
  ppm = 5,
  ms1_spectrum = ms1_spectrum,
  min_score = 0
)

candidate_view <- candidates[, c(
  "observed_mz", "expected_mz", "adduct", "ppm_error", "mass_score",
  "rule_prior", "chemical_score", "isotope_score",
  "adduct_network_score", "score", "ramp_ids"
)]
knitr::kable(candidate_view, digits = 4, row.names = FALSE)
```

| observed_mz | expected_mz | adduct | ppm_error | mass_score | rule_prior | chemical_score | isotope_score | adduct_network_score | score | ramp_ids |
|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|:---|
| 149.0446 | 149.0444 | M+H | 1.2 | 0.7717 | 1.00 | 1 | 1.0 | NA | 0.7717 | RAMP_C_000219574 |
| 171.0262 | 171.0264 | M+Na | 1.4 | 0.7027 | 0.90 | 1 | 0.5 | NA | 0.3162 | RAMP_C_000219574 |
| 193.0085 | 193.0083 | M+2Na-H | 0.6 | 0.9373 | 0.65 | 1 | 0.5 | 1 | 0.3046 | RAMP_C_000219574 |
| 297.0814 | 297.0816 | 2M+H | 0.8 | 0.8912 | 0.55 | 1 | 0.5 | 1 | 0.2451 | RAMP_C_000219574 |

Candidates outside the requested tolerance are removed before scoring:

``` r

outside_five_ppm <-
  ion_table$expected_mz[ion_table$adduct == "M+H"] * (1 + 5.1e-6)

data.frame(
  query = "5.1 ppm outside expected m/z",
  retained_candidates = nrow(QueryMZAnnotationIndex(
    outside_five_ppm,
    positive_index,
    ppm = 5
  ))
) |>
  knitr::kable(row.names = FALSE)
```

| query                        | retained_candidates |
|:-----------------------------|--------------------:|
| 5.1 ppm outside expected m/z |                   0 |

The final score is multiplicative:

``` math
S = S_{mass} \times P_{rule} \times S_{chemical}
    \times S_{isotope} \times S_{adduct\ family}.
```

Missing optional MS1 context is neutral rather than evidence of absence.
When a contextual spectrum is supplied, missing expected isotope or
base-adduct peaks can reduce the score.

### Adduct-family support

A complex ion is more credible when a corresponding base monomer is
present in the same spectrum. Here isotope scoring is disabled so that
only adduct-network support changes.

``` r

dimer_mz <- ion_table$expected_mz[ion_table$adduct == "2M+H"]
monomer_mz <- ion_table$expected_mz[ion_table$adduct == "M+H"]

isolated_dimer <- QueryMZAnnotationIndex(
  dimer_mz,
  positive_index,
  ppm = 5,
  ms1_spectrum = data.frame(mz = dimer_mz, intensity = 100),
  check_isotopes = FALSE
)
supported_dimer <- QueryMZAnnotationIndex(
  dimer_mz,
  positive_index,
  ppm = 5,
  ms1_spectrum = data.frame(
    mz = c(dimer_mz, monomer_mz),
    intensity = c(100, 500)
  ),
  check_isotopes = FALSE
)

data.frame(
  context = c("isolated dimer", "base monomer present"),
  adduct_network_score = c(
    isolated_dimer$adduct_network_score,
    supported_dimer$adduct_network_score
  ),
  final_score = c(isolated_dimer$score, supported_dimer$score)
) |>
  knitr::kable(digits = 3, row.names = FALSE)
```

| context              | adduct_network_score | final_score |
|:---------------------|---------------------:|------------:|
| isolated dimer       |                  0.1 |       0.055 |
| base monomer present |                  1.0 |       0.550 |

### Mass-defect heuristic

For a nominally CHO-only candidate, an unusually negative signed mass
defect is treated as background-like evidence. It receives a penalty
rather than being deleted, because halogenated, metal-containing, and
unusual biological ions must not be rejected by a universal
fractional-mass rule.

``` r

background_like_db <- data.frame(
  formula = "C2H2O2",
  exactmass = 99.942723533379,
  id = "synthetic-background-like-candidate",
  name = "Synthetic background-like candidate"
)
background_index <- BuildMZAnnotationIndex(
  background_like_db,
  polarity = "positive",
  adducts = "M+H"
)
background_mz <- background_index$expected_mz[[1L]]

with_penalty <- QueryMZAnnotationIndex(
  background_mz,
  background_index,
  ppm = 0,
  use_mass_defect = TRUE,
  check_isotopes = FALSE,
  check_adduct_network = FALSE
)
without_penalty <- QueryMZAnnotationIndex(
  background_mz,
  background_index,
  ppm = 0,
  use_mass_defect = FALSE,
  check_isotopes = FALSE,
  check_adduct_network = FALSE
)

data.frame(
  setting = c("mass-defect heuristic on", "mass-defect heuristic off"),
  signed_mass_defect = background_mz - round(background_mz),
  chemical_score = c(
    with_penalty$chemical_score,
    without_penalty$chemical_score
  ),
  retained = TRUE
) |>
  knitr::kable(digits = 3, row.names = FALSE)
```

| setting                   | signed_mass_defect | chemical_score | retained |
|:--------------------------|-------------------:|---------------:|:---------|
| mass-defect heuristic on  |              -0.05 |            0.2 | TRUE     |
| mass-defect heuristic off |              -0.05 |            1.0 | TRUE     |

## Annotate individual values or a SpaMTP object

[`AnnotateMZ()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/AnnotateMZ.md)
is a convenience wrapper for individual or small batches of observed
values. Supply a pre-built index to avoid rebuilding it.

``` r

single_result <- AnnotateMZ(
  observed_mz[["M+H"]],
  index = positive_index,
  polarity = "positive",
  ppm = 5,
  ms1_spectrum = ms1_spectrum,
  min_score = 0
)
knitr::kable(
  single_result[, c(
    "observed_mz", "metabolite_names", "ramp_ids", "adduct",
    "ppm_error", "score"
  )],
  digits = 4,
  row.names = FALSE
)
```

| observed_mz | metabolite_names         | ramp_ids         | adduct | ppm_error |  score |
|------------:|:-------------------------|:-----------------|:-------|----------:|-------:|
|    149.0446 | D-2-Hydroxyglutaric acid | RAMP_C_000219574 | M+H    |       1.2 | 0.7717 |

[`AnnotateSM()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/AnnotateSM.md)
applies the same engine to m/z feature metadata and preserves the ranked
candidate table in `@tools$mz_annotation`. The small object below is
only for demonstrating storage and provenance.

``` r

demo_counts <- Matrix::Matrix(
  matrix(
    seq_len(length(observed_mz) * 4L),
    nrow = length(observed_mz),
    ncol = 4L
  ),
  sparse = TRUE
)
rownames(demo_counts) <- paste0("mz-", observed_mz)
colnames(demo_counts) <- paste0("pixel_", seq_len(ncol(demo_counts)))

demo_object <- SeuratObject::CreateSeuratObject(
  counts = demo_counts,
  assay = "SPM"
)
feature_metadata <- data.frame(
  raw_mz = unname(observed_mz),
  row.names = rownames(demo_object[["SPM"]])
)
demo_object[["SPM"]] <- SeuratObject::AddMetaData(
  demo_object[["SPM"]],
  metadata = feature_metadata
)

annotated_object <- AnnotateSM(
  demo_object,
  db = NULL,
  assay = "SPM",
  raw.mz.column = "raw_mz",
  ppm_error = 5,
  adducts = selected_adducts,
  polarity = "positive",
  return.only.annotated = FALSE,
  save.intermediate = TRUE,
  min_score = 0,
  verbose = FALSE,
  index = positive_index,
  ms1_spectrum = ms1_spectrum
)

annotation_columns <- annotated_object[["SPM"]][[]][, c(
  "raw_mz", "all_IsomerNames", "all_Adducts", "all_Scores",
  "all_Ramp_IDs"
)]
knitr::kable(annotation_columns, digits = 4)
```

|  | raw_mz | all_IsomerNames | all_Adducts | all_Scores | all_Ramp_IDs |
|:---|---:|:---|:---|:---|:---|
| mz-149.044628685961 | 149.0446 | D-2-Hydroxyglutaric acid | M+H | 0.7717 | RAMP_C_000219574 |
| mz-171.026154631139 | 171.0262 | D-2-Hydroxyglutaric acid | M+Na | 0.3162 | RAMP_C_000219574 |
| mz-193.008454108564 | 193.0085 | D-2-Hydroxyglutaric acid | M+2Na-H | 0.3046 | RAMP_C_000219574 |
| mz-297.081385533322 | 297.0814 | D-2-Hydroxyglutaric acid | 2M+H | 0.2451 | RAMP_C_000219574 |

[`AnnotationInfo()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/AnnotationInfo.md)
records which engine and database provenance downstream functions will
consume.

``` r

annotation_info <- AnnotationInfo(annotated_object)
annotation_info_table <- data.frame(
  field = names(annotation_info),
  value = vapply(annotation_info, paste, collapse = "; ", character(1)),
  row.names = NULL
)
knitr::kable(annotation_info_table, row.names = FALSE)
```

| field          | value                      |
|:---------------|:---------------------------|
| schema_version | 2                          |
| engine         | indexed-chemical-v2        |
| ramp_version   | 3.0.7                      |
| generated_at   | 2026-08-11 02:59:15 UTC    |
| candidates     | 4                          |
| assay          | SPM                        |
| raw_mz_column  | raw_mz                     |
| polarity       | positive                   |
| ppm            | 5                          |
| min_score      | 0                          |
| adducts        | M+H; M+Na; M+2Na-H; 2M+H   |
| database       | pre-built annotation index |

## Choose thresholds without re-annotation

The recommended production pattern is to retain all ppm-valid candidates
using `min_score = 0`, then select a threshold appropriate to the
analysis. This makes coverage and specificity tunable without repeating
index construction.

``` r

stored_candidates <- annotated_object@tools$mz_annotation$results
score_thresholds <- c(0, 0.2, 0.5)
threshold_summary <- data.frame(
  score_threshold = score_thresholds,
  retained_candidates = vapply(
    score_thresholds,
    function(threshold) sum(stored_candidates$Score >= threshold),
    integer(1)
  )
)
knitr::kable(threshold_summary, row.names = FALSE)
```

| score_threshold | retained_candidates |
|----------------:|--------------------:|
|             0.0 |                   4 |
|             0.2 |                   4 |
|             0.5 |                   1 |

For pathway analysis, the annotation threshold and
differential-significance threshold answer different questions. Setting
`pval_cutoff_mets = 1` allows non-DE metabolites to contribute, while
`annotation_score_threshold` still controls which chemical candidates
are eligible.

``` r

regional_pathways <- FindRegionalPathways(
  object,
  analyte_types = c("genes", "metabolites"),
  ident = "region",
  DE.list = differential_results,
  pval_cutoff_mets = 1,
  annotation_score_threshold = 0.05
)

PathwayNetworkPlots(
  object,
  ident = "region",
  regpathway = regional_pathways,
  DE.list = differential_results,
  annotation_score_threshold = 0.05,
  metabolite_detection = "annotated"
)
```

## Production use with the complete RaMP table

For one call, `AnnotateSM(db = NULL)` automatically uses the bundled
RaMP chemical-property table. For repeated calls, explicitly build and
reuse an index. Restricting the adduct set to ions supported by the
experiment reduces memory use and false hypotheses.

``` r

ramp_positive_index <- BuildMZAnnotationIndex(
  chem_props,
  polarity = "positive",
  adducts = c("M+H", "M+Na", "M+K")
)

object <- AnnotateSM(
  object,
  db = NULL,
  index = ramp_positive_index,
  assay = "SPM",
  polarity = "positive",
  ppm_error = 5,
  adducts = c("M+H", "M+Na", "M+K"),
  min_score = 0,
  return.only.annotated = FALSE
)
```

## Interpretation and limitations

- An annotation is a ranked chemical hypothesis, not a confirmed
  identity, probability, or false-discovery rate.
- Accurate mass alone generally cannot resolve structural isomers.
- Isotope and adduct-family checks require a spectrum from an
  appropriate retention-time, pixel, or regional context.
- The mass-defect rule is a score penalty, not a universal exclusion
  rule.
- Exchangeable-proton pruning requires the relevant database property;
  it is not inferred from elemental formula alone.
- Orthogonal evidence such as authentic standards, retention time, ion
  mobility, or MS/MS remains necessary for high-confidence
  identification.

The
[`Annotation Refinement`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/articles/Annotation_Refinement.md)
workflow can then use spatial, pathway, lipid-class, or pseudo-MS/MS
information to further prioritise the candidates stored here.

``` r

sessionInfo()
```

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8    
    ##  [5] LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8    LC_PAPER=C.UTF-8       LC_NAME=C             
    ##  [9] LC_ADDRESS=C           LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] SpaMTP_1.1.0.9000
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3     rstudioapi_0.19.0      jsonlite_2.0.0         magrittr_2.0.5        
    ##   [5] spatstat.utils_3.2-4   farver_2.1.2           rmarkdown_2.31         fs_2.1.0              
    ##   [9] ragg_1.5.2             vctrs_0.7.3            ROCR_1.0-12            memoise_2.0.1         
    ##  [13] spatstat.explore_3.8-2 htmltools_0.5.9        sass_0.4.10            sctransform_0.4.3     
    ##  [17] parallelly_1.48.0      KernSmooth_2.23-26     bslib_0.12.0           htmlwidgets_1.6.4     
    ##  [21] desc_1.4.3             ica_1.0-3              plyr_1.8.9             plotly_4.12.1         
    ##  [25] zoo_1.9-0              cachem_1.1.0           whisker_0.4.1          igraph_2.3.3          
    ##  [29] mime_0.13              lifecycle_1.0.5        pkgconfig_2.0.3        Matrix_1.7-5          
    ##  [33] R6_2.6.1               fastmap_1.2.0          fitdistrplus_1.2-6     future_1.75.0         
    ##  [37] shiny_1.14.0           digest_0.6.39          S4Vectors_0.50.1       patchwork_1.3.2       
    ##  [41] Seurat_5.5.1           tensor_1.5.1           RSpectra_0.16-2        irlba_2.3.7           
    ##  [45] textshaping_1.0.5      progressr_1.0.0        spatstat.sparse_3.2-0  httr_1.4.8            
    ##  [49] polyclip_1.10-7        abind_1.4-8            compiler_4.6.1         proxy_0.4-29          
    ##  [53] withr_3.0.3            S7_0.2.2               BiocParallel_1.46.0    DBI_1.3.0             
    ##  [57] fastDummies_1.7.6      MASS_7.3-65            classInt_0.4-11        units_1.0-1           
    ##  [61] tools_4.6.1            lmtest_0.9-40          otel_0.2.0             httpuv_1.6.17         
    ##  [65] future.apply_1.20.2    goftest_1.2-3          glue_1.8.1             nlme_3.1-169          
    ##  [69] promises_1.5.0         sf_1.1-2               grid_4.6.1             Rtsne_0.17            
    ##  [73] cluster_2.1.8.2        reshape2_1.4.5         generics_0.1.4         gtable_0.3.6          
    ##  [77] spatstat.data_3.1-9    class_7.3-23           tidyr_1.3.2            data.table_1.18.4     
    ##  [81] sp_2.2-3               xml2_1.6.0             BiocGenerics_0.58.1    spatstat.geom_3.8-2   
    ##  [85] RcppAnnoy_0.0.23       ggrepel_0.9.8          RANN_2.6.2             pillar_1.11.1         
    ##  [89] stringr_1.6.0          spam_2.11-4            RcppHNSW_0.7.0         limma_3.68.4          
    ##  [93] later_1.4.8            splines_4.6.1          dplyr_1.2.1            lattice_0.22-9        
    ##  [97] survival_3.8-6         deldir_2.0-4           tidyselect_1.2.1       CardinalIO_1.10.0     
    ## [101] miniUI_0.1.2           pbapply_1.7-4          downlit_0.4.5          knitr_1.51            
    ## [105] gridExtra_2.3.1        ProtGenerics_1.44.0    matter_2.14.0          scattermore_1.2       
    ## [109] stats4_4.6.1           xfun_0.60              Biobase_2.72.0         statmod_1.5.2         
    ## [113] matrixStats_1.5.0      stringi_1.8.9          yaml_2.3.12            evaluate_1.0.5        
    ## [117] codetools_0.2-20       tibble_3.3.1           cli_3.6.6              ontologyIndex_2.12    
    ## [121] uwot_0.2.4             xtable_1.8-8           reticulate_1.46.0      systemfonts_1.3.2     
    ## [125] jquerylib_0.1.4        Rcpp_1.1.2             globals_0.19.1         spatstat.random_3.5-1 
    ## [129] zeallot_0.2.0          png_0.1-9              spatstat.univar_3.2-0  parallel_4.6.1        
    ## [133] pkgdown_2.2.1          ggplot2_4.0.3          dotCall64_1.2          listenv_1.0.0         
    ## [137] viridisLite_0.4.3      e1071_1.7-17           scales_1.4.0           ggridges_0.5.7        
    ## [141] SeuratObject_5.4.0     Cardinal_3.14.0        purrr_1.2.2            rlang_1.3.0           
    ## [145] cowplot_1.2.0          shinyjs_2.1.1
