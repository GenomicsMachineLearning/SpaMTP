# Find all library spectra10 above a cosine threshold for each imaging peak, restricting to a global precursor-mz range.

Find all library spectra10 above a cosine threshold for each imaging
peak, restricting to a global precursor-mz range.

## Usage

``` r
Pseudo_msms(
  raw_mz,
  counts,
  spectra10_list,
  ppm_tol = 10,
  frag_tol_da = 0.01,
  cos_threshold = 0.7,
  min_precursor = 0,
  max_precursor = Inf
)
```

## Arguments

- raw_mz:

  Numeric vector of length n_peaks

- counts:

  Numeric matrix n_peaks × n_pixels

- spectra10_list:

  List of Spectrum2 objects (length N_lib)

- ppm_tol:

  Precursor‐matching tolerance in ppm (default 10)

- frag_tol_da:

  Fragment‐matching tolerance in Da (default 0.01)

- cos_threshold:

  Minimum cosine similarity to report (default 0.7)

- min_precursor:

  Minimum allowed precursor m/z (default = 0, i.e. no lower bound)

- max_precursor:

  Maximum allowed precursor m/z (default = Inf, i.e. no upper bound)

## Value

A list of length n_peaks. Each element is an integer vector giving the
indices in spectra10_list of all library spectra10 with cosine ≥
cos_threshold (or integer(0) if none).
