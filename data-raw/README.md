# Rebuilding the pruned RaMP snapshot

SpaMTP does not distribute the full RaMP SQLite database. It stores only the
columns used by its annotation, pathway, and visualisation code as compressed R
data objects.

Download or obtain the official decompressed RaMP 3.x SQLite file, then run from
the package root:

```sh
Rscript data-raw/update_ramp_data.R /path/to/RaMP_SQLite_v3.0.7.sqlite
```

The script validates the upstream schema, preserves the existing SpaMTP object
and column contracts, updates the five pruned objects, and writes
`ramp_db_metadata.rda` with database/source versions and row counts. The large
SQLite file is a build input and must not be committed to SpaMTP.
