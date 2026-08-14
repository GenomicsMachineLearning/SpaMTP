# Metadata for the bundled pruned RaMP snapshot

Version, source database releases, row counts, and retained columns for
the SpaMTP RaMP-DB 3.0.7 snapshot. The full upstream SQLite database is
not bundled with SpaMTP.

## Usage

``` r
ramp_db_metadata
```

## Format

A named list containing `ramp_version`, `load_timestamp`,
`version_notes`, `source_versions`, `snapshot_counts`,
`upstream_repository`, `upstream_database_file`, `pruning`, and graph-ID
harmonisation metadata.

## Source

<https://github.com/ncats/RaMP-DB>
