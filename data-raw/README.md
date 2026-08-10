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
`ramp_db_metadata.rda` with database/source versions and row counts. It retains
the previous `source_df` in memory and uses stable source identifiers (HMDB,
ChEBI, KEGG, LIPIDMAPS, gene symbols, UniProt, and related IDs) to migrate the
stored pathway-graph node IDs to the new RaMP release.

If the stored graph objects came from an older release than the currently
bundled `source_df`, supply the graph-era source table explicitly:

```sh
Rscript data-raw/update_ramp_data.R \
  /path/to/RaMP_SQLite_v3.0.7.sqlite \
  /path/to/graph-era/source_df.rda
```

Graph harmonisation can also be run or previewed independently:

```sh
Rscript data-raw/update_ramp_graph_data.R /path/to/graph-era/source_df.rda --dry-run
Rscript data-raw/update_ramp_graph_data.R /path/to/graph-era/source_df.rda
Rscript data-raw/audit_ramp_graph_data.R
```

The audit exits with a non-zero status when a graph node is absent from the
bundled `analyte` table or lacks both a current name and a stable source-ID
label. The large SQLite file is a build input and must not be committed to
SpaMTP.
