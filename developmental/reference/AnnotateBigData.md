# Annotates vector of m/z values

This function assigns each valid m/z peak with one/multiple metabolite
names based on the mass difference between the observed value and the
theoretical value documented in the reference database. This function is
to be used when dealing with large datasets as a preprocessing step.
Users can annotate m/z values first and then subset their data
accordinly before loading it into a SpaMTP Seurat Object.

## Usage

``` r
AnnotateBigData(
  mzs,
  db = NULL,
  ppm_error = NULL,
  adducts = NULL,
  polarity = NULL,
  tof_resolution = 30000,
  verbose = TRUE,
  maldi_matrix = NULL,
  ...
)
```

## Arguments

- mzs:

  Vector containing m/z values for annotation.

- db:

  Reference metabolite dataset in the form of a data.frame. SpaMTP
  provides four pre-cleaned databases (`HMDB_db`, `Lipidmaps_db`,
  `Chebi_db`, `GNPS_db`). May be `NULL` when a pre-built `index` is
  supplied through `...`.

- ppm_error:

  Mass tolerance in ppm. If `NULL`, a strict 5 ppm maximum is used (or a
  smaller value inferred from `tof_resolution`). Set to zero for exact
  numerical matches.

- adducts:

  Optional adduct names/notations; see
  [`AdductRules()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/AdductRules.md).
  If `NULL`, use the complete rule space selected from `maldi_matrix`,
  or all validated general rules for the selected polarity when no
  matrix is given.

- polarity:

  Character string defining the ion mode. When `NULL`, infer it from a
  supplied index or MALDI matrix profile, otherwise use positive.

- tof_resolution:

  Instrument resolving power retained for compatibility; it can only
  tighten, not widen, the default 5 ppm mass-accuracy threshold.

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

- maldi_matrix:

  Optional MALDI matrix/reagent profile used for automatic rule
  selection. The `adducts` argument remains an optional restriction.

- ...:

  Additional indexed annotation/scoring arguments passed to
  [`annotateTable()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/annotateTable.md),
  such as `index`, `rules`, or `ms1_spectrum`.

## Value

A data.frame containing all successfully annotated m/z values, with
their corresponding annotation.

## Examples

``` r
#cardinal <- readImzML("./Test_Data/Spotted/test_data1")
#mzs <- data.frame(Cardinal::featureData(cardinal))$mz
#results <- AnnotateBigData(mzs, db = HMDB_db, ppm_error = 3, adducts = c("M-H", "M+Cl"), polarity = "negative")
#cardinal_subset <- Cardinal::subset(cardinal, mz %in% results$observed_mz)
#SpaMTP_data <- CardinalToSeurat(cardinal_subset)
```
