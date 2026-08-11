# Annotates m/z values stored in a SpaMTP Object

This function assigns each valid m/z peak with one/multiple metabolite
names based on the mass difference between the observed value and the
theoretical value documented in the reference database. SpaMTP contains
4 cleaned reference databases to choose from these include HMDB, Lipid
Maps, ChEBI and GNPS. These databases can also be combined for increased
coverage.

## Usage

``` r
AnnotateSM(
  data,
  db = NULL,
  assay = "Spatial",
  raw.mz.column = "raw_mz",
  ppm_error = NULL,
  adducts = NULL,
  polarity = NULL,
  tof_resolution = 30000,
  filepath = NULL,
  return.only.annotated = TRUE,
  save.intermediate = TRUE,
  min_score = 0,
  verbose = TRUE,
  maldi_matrix = NULL,
  ...
)
```

## Arguments

- data:

  Seurat Spatial Metabolomic Object containing m/z values for
  annotation.

- db:

  Reference metabolite dataset in the form of a data.frame. When `NULL`,
  the bundled current RaMP `chem_props` table is used, unless a
  pre-built `index` is supplied through `...`.

- assay:

  Character string defining the Seurat assay which contains the mz
  counts being annotated (default = "Spatial").

- raw.mz.column:

  Character string defining the Seurat assay slot which contains the raw
  mz values, this is without the 'mz-' and are a vector of integers.
  This is setup by default when running the cardinal_to_seurat()
  function (default = "raw_mz").

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

  Character string defining the ion mode. When `NULL`, use the MALDI
  matrix profile default, a supplied index's mode, or positive mode when
  neither is available.

- tof_resolution:

  Instrument resolving power retained for compatibility; it can only
  tighten, not widen, the default 5 ppm mass-accuracy threshold.

- filepath:

  Character string of the directory to store the
  \_annotated_mz_peaks.csv. If set to NULL no dataframe will be saved
  (default = NULL).

- return.only.annotated:

  Boolean value indicating if the annotated Seurat Object should only
  include m/z values that were successfully annotated (default = TRUE).

- save.intermediate:

  Boolean indicating whether to store the scored annotation result and
  its pipeline/RaMP provenance in `@tools$mz_annotation`. A
  compatibility copy is retained in `@tools$db_3` (default = TRUE).

- min_score:

  Minimum annotation score retained in the stored candidate table. The
  default `0` keeps all ppm-valid candidates so downstream pathway
  functions can apply a user-defined threshold without re-running
  annotation.

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

- maldi_matrix:

  Optional MALDI matrix or derivatization reagent name. When supplied,
  SpaMTP selects validated matrix-specific rules automatically.
  `adducts` is optional and only restricts that automatic search space
  when explicitly supplied.

- ...:

  Additional indexed annotation/scoring arguments passed to
  [`annotateTable()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/annotateTable.md),
  such as `index`, `rules`, or `ms1_spectrum`.

## Value

A Seurat Object with m/z values annotated. These annotations are stored
in the relative assay's meta.data (e.g. SeuratObj`[["Spatial"]][[]]`)

## Examples

``` r
# HMDB_db <- load("data/HMDB_1_names.rds")
# Annotated_SeuratObj <- AnnotateSM(SeuratObj, HMDB_db)
```
