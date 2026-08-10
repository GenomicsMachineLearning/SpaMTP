# Query an indexed metabolite annotation search space

Candidates are pruned by ppm tolerance and proton/charge validity during
index construction, then ranked using mass error, rule prior,
mass-defect, isotope, and contextual adduct-family evidence.

## Usage

``` r
QueryMZAnnotationIndex(
  observed_mz,
  index,
  ppm = 5,
  ms1_spectrum = NULL,
  use_mass_defect = TRUE,
  check_isotopes = TRUE,
  check_adduct_network = TRUE,
  min_score = 0
)
```

## Arguments

- observed_mz:

  Numeric vector of observed m/z values.

- index:

  A `spamtp_mz_index` from
  [`BuildMZAnnotationIndex()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/BuildMZAnnotationIndex.md).

- ppm:

  Mass tolerance in parts per million.

- ms1_spectrum:

  Optional contextual spectrum with `mz` and `intensity` columns. It
  should represent the same retention-time window or spatial
  pixel/region as the queried peaks.

- use_mass_defect:

  Apply the CHO negative mass-defect penalty.

- check_isotopes:

  Score carbon-13 and, when relevant, chlorine/bromine M+2 patterns.

- check_adduct_network:

  Downweight complex ions whose base monomer family is absent from
  `ms1_spectrum`.

- min_score:

  Minimum final score to retain.

## Value

A ranked candidate data frame.
