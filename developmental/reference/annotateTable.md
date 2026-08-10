# Annotates m/z values stored in a data.frame based on reference metabolite dataset

Helper function for
[`AnnotateSM()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/AnnotateSM.md)
and
[`FishersPathwayAnalysis()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/FishersPathwayAnalysis.md).

## Usage

``` r
annotateTable(
  mz_df,
  db = NULL,
  ppm_error = NULL,
  adducts = NULL,
  polarity = "positive",
  tof_resolution = 30000,
  verbose = TRUE,
  index = NULL,
  rules = NULL,
  ms1_spectrum = NULL,
  use_mass_defect = TRUE,
  check_isotopes = TRUE,
  check_adduct_network = TRUE,
  min_score = 0
)
```

## Arguments

- mz_df:

  dataframe containing m/z values for annotation.

- db:

  Reference metabolite dataset in the form of a data.frame. May be
  `NULL` when `index` is supplied.

- ppm_error:

  Mass tolerance in ppm. If `NULL`, a strict 5 ppm maximum is used (or a
  smaller value inferred from `tof_resolution`). Set to zero for exact
  numerical matches.

- adducts:

  Optional adduct names/notations; see
  [`AdductRules()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/AdductRules.md).
  If `NULL`, all validated rules for the selected polarity are used.

- polarity:

  Character string defining the polarity of adducts to use, either
  "positive", "negative" or "neutral" (default = "positive").

- tof_resolution:

  Instrument resolving power retained for compatibility; it can only
  tighten, not widen, the default 5 ppm mass-accuracy threshold.

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

- index:

  Optional reusable index created by
  [`BuildMZAnnotationIndex()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/BuildMZAnnotationIndex.md).

- rules:

  Optional custom adduct rule table; see
  [`AdductRules()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/AdductRules.md).

- ms1_spectrum:

  Optional contextual MS1 spectrum with `mz` and `intensity` columns for
  isotope and adduct-family scoring.

- use_mass_defect:

  Apply the CHO negative mass-defect penalty.

- check_isotopes:

  Score carbon-13 and halogen isotope evidence when a contextual
  spectrum is supplied.

- check_adduct_network:

  Downweight complex ions whose base monomer family is absent from the
  contextual spectrum.

- min_score:

  Minimum final annotation score to retain.

## Value

Generates an intermediate annotated m/z dataframe

## Examples

``` r

### HelperFunction
```
