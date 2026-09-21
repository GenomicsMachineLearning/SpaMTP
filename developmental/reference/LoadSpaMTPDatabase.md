# Load versioned SpaMTP annotation resources

Loads annotation resources bundled with SpaMTP without another data
package, a Hub lookup, or a download. Resources are cached for the
current R session. A named custom bundle or local RDS directory can be
supplied for user-curated workflows.

## Usage

``` r
LoadSpaMTPDatabase(
  resources = c("chem_props", "source_df", "analyte", "analytehaspathway", "pathway"),
  version = "latest",
  source = c("auto", "bundled", "local"),
  database = NULL,
  local_dir = NULL,
  hub = NULL,
  offline = FALSE,
  refresh = FALSE
)
```

## Arguments

- resources:

  Character vector of resource names. Use
  [`SpaMTPDatabaseInfo()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/SpaMTPDatabaseInfo.md)
  to list valid names.

- version:

  Database snapshot version, or `"latest"` for the bundled snapshot. For
  local files, declared metadata must match an explicit version.

- source:

  Database source. `"auto"` uses bundled data unless `local_dir` is
  supplied. `"bundled"` uses installed data; `"local"` requires local
  RDS files and never falls back to another source.

- database:

  Optional named list containing the requested resources. When supplied,
  no Hub lookup is performed.

- local_dir:

  Optional directory containing `<resource>.rds` files, either directly
  or in a version-named subdirectory. Without local version metadata,
  `"latest"` files are labelled `"unversioned"`.

- hub:

  Retained for call compatibility; must be `NULL` in this release.

- offline:

  Retained for call compatibility. All resource loading is offline,
  regardless of this argument.

- refresh:

  If `TRUE`, bypass SpaMTP's in-session resource cache.

## Value

A named list containing the requested resources.

## Details

Bundled RaMP 3.0.7 graphs already include corrected interaction labels
and directions. Exact affected local or custom topology resources
receive a checksum-guarded correction for historical interaction-code
and direction recycling. Source labels are retained in
`source_reaction_type`; the resource attribute
`spamtp_interaction_repair` records the correction separately from the
original download metadata. Published files are unchanged. See
[`vignette("Pathway_Database_Integration", package = "SpaMTP")`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/articles/Pathway_Database_Integration.md)
for provenance and reproducible build instructions.

## Examples

``` r
utils::str(formals(LoadSpaMTPDatabase))
#> Dotted pair list of 8
#>  $ resources: language c("chem_props", "source_df", "analyte", "analytehaspathway", "pathway")
#>  $ version  : chr "latest"
#>  $ source   : language c("auto", "bundled", "local")
#>  $ database : NULL
#>  $ local_dir: NULL
#>  $ hub      : NULL
#>  $ offline  : logi FALSE
#>  $ refresh  : logi FALSE
example_database <- list(
  ramp_db_metadata = list(ramp_version = "example")
)
database <- LoadSpaMTPDatabase(
  "ramp_db_metadata",
  database = example_database
)
names(database)
#> [1] "ramp_db_metadata"
```
