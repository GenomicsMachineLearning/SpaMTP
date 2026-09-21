# Pathway databases, interaction provenance, and rebuilding graphs

## Membership and topology are separate resources

SpaMTP combines identifiers and pathway memberships from [NCATS
RaMP-DB](https://github.com/ncats/RaMP-DB) with network topology from
[graphite](https://bioconductor.org/packages/graphite/). The `RAMP_*`
identifiers on graph nodes describe the common identifier system; the
interaction labels are derived from graphite’s source graphs. They are
not RaMP interaction IDs.

| Resource | Contents and origin |
|----|----|
| `analyte`, `source_df` | RaMP entities and their source identifier cross-references |
| `pathway`, `analytehaspathway` | Source-specific pathway descriptions and RaMP entity memberships |
| `chem_props` | Pruned RaMP chemical properties used for metabolite annotation |
| `ramp_kegg` | Human KEGG topology from graphite |
| `ramp_reactome` | Human Reactome topology from graphite |
| `ramp_wikipathway` | Human WikiPathways topology from graphite |
| `ramp_hmdb` | Human **SMPDB** topology from graphite, retained under its historical SpaMTP resource name |

The RaMP 3.0.7 membership snapshot and the graph snapshot have separate
dates. The currently supported topology correction uses graphite archive
19, dated 2022-11-01. Updating RaMP identifiers does not update the
biological topology. SpaMTP does not construct a consensus of the four
networks or treat equally named pathways as equivalent biological
definitions. Membership tables retain their source identifiers. The
network viewer selects a source graph; its legacy name-based selection
can collapse identically named result rows, so filter results to a
single database when comparing such pathways.

## Identifier consolidation

`data-raw/update_ramp_data.R` extracts the required columns from RaMP’s
SQLite `analyte`, `source`, `analytehaspathway`, `pathway`, and
`chem_props` tables. It records RaMP and upstream database versions in
`ramp_db_metadata`.

The original topology builder translated graphite namespaces, for
example `ENTREZID` to `entrez:` and `KEGGCOMP` to `kegg:`, then matched
those source IDs to RaMP identifiers. The historical builder selected
the first RaMP match; unmapped endpoints remain missing and the network
viewer excludes these edges. That historical mapping is reproduced by
the correction recipe to preserve the existing nodes, rather than
introducing new biological identifier choices.

`data-raw/ramp_graph_utils.R` migrates retired graph IDs to the current
RaMP snapshot using the maximum number of shared stable source
identifiers. Tied matches expand to multiple current entities. Two
retired heme identifiers use explicit, validated replacements. The
migration preserves interaction columns and deduplicates complete rows,
keeping different interactions between the same pair of nodes distinct.

## Interaction correction

The original graphite-to-SpaMTP conversion used code equivalent to:

``` r

ifelse(length(edge$type) == 0, NA, edge$type)
```

The test has length one, so the result contains only the first factor
code. Creating the edge table recycled that value across all rows. The
same issue affected `direction`. In addition, a factor code has meaning
only within its source’s factor levels; Reactome, WikiPathways, and
SMPDB use different levels from KEGG. Interpreting all four using KEGG’s
numeric dictionary was incorrect.

The corrected conversion in `data-raw/pathway_interaction_utils.R` reads
**every source label and direction**, maps labels to stable SpaMTP style
codes, and retains the exact label in `source_reaction_type`. Code 4
still means inhibition. Generic biochemical processes are not relabelled
as activation or inhibition. Undirected edges are displayed without an
arrowhead. These are source annotations, not causal relationships
inferred from expression data.

For KEGG Cell Cycle (`hsa:04110`), the repaired protein edge table has
1,009 rows and 11 interaction types, including 286 inhibition rows, 113
activation rows, 84 phosphorylation rows, and 224 binding rows. The
original table incorrectly assigned all 1,009 rows to inhibition.

The correction covers all four topology collections. Restoring labels
before ID migration also recovers parallel interactions that were
incorrectly collapsed when their corrupted labels became identical. The
network viewer draws parallel interactions on separate curves, including
opposite-direction edges and self-loops. Hover over each curve to
inspect its original source label; undirected interactions have no
arrowhead.

## Loading and inspecting corrected graphs

Use SpaMTP’s loader to obtain corrected graphs:

``` r

db <- LoadSpaMTPDatabase("ramp_kegg")
cell_cycle <- db$ramp_kegg[["Cell cycle"]]
table(cell_cycle$protEdges$source_reaction_type)
attr(db$ramp_kegg, "spamtp_interaction_repair")
```

The bundled graphs already include the correction. The database loader
and network workflows also repair local RDS files and custom bundles
containing the exact affected snapshot. Compact correction files in
`inst/extdata/*_interactions_v1.rds` are applied only when the
resource’s serialized contents match its pinned SHA-256 fingerprint.
Other custom graphs and future snapshots are not modified. Corrected
graphs carry `spamtp_interaction_repair` provenance and can be saved
with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html).

The immutable external 3.0.7 RDS files and their checksums remain
unchanged. Reading those original files directly bypasses SpaMTP’s
correction; supplying them through `LoadSpaMTPDatabase(local_dir = ...)`
applies it. No companion package is needed. The `spamtp_database`
attribute identifies the bundled or local resource;
`spamtp_interaction_repair` describes the additional in-memory
correction.

The membership-based enrichment functions do not use these edge labels.
This bug affects network interactions and any downstream analysis that
directly uses those interactions, rather than changing RaMP pathway
memberships.

## Reproducing the correction

The maintainer recipe is `data-raw/rebuild_pathway_interactions.R` in
the source checkout. Its input URLs and MD5 checksums are pinned in
`inst/extdata/pathway_interactions_provenance.csv`. It requires `digest`
plus base R; installing graphite is unnecessary because it reads
archived data.

1.  Download the four graphite archive-19 graphs and the pinned original
    RaMP source mapping using the URLs and filenames in the provenance
    CSV.
2.  Stage the listed immutable SpaMTPdb 3.0.7 resources.
3.  From the source checkout, run:

``` sh
Rscript data-raw/rebuild_pathway_interactions.R \
  /path/to/graphite-19 /path/to/initial-source_df.rda \
  /path/to/SpaMTPdb-resources/3.0.7
```

The recipe verifies input checksums, reproduces the old faulty
conversion, and requires every reconstructed endpoint, direction, and
code to match the stored input **before** generating corrections. It
then reconverts from source labels, repeats the same ID migration, and
stores corrected row references, directions, style codes, and source
labels. No new endpoints are invented.
`inst/extdata/pathway_interactions_audit.csv` records all pathways and
edge tables checked; edge counts there include historical all-missing
empty-table sentinels. Regression tests cover vector preservation,
factor order, source semantics, parallel interactions, input
fingerprints, and the staged Cell Cycle resource. Set
`SPAMTPDB_RESOURCE_DIR` to run the external-resource test.

For a future graph release, build from labelled source edges using the
conversion utilities before running ID migration. Publish rebuilt data
under a new resource version; do not overwrite the immutable 3.0.7 files
or reuse this correction’s fingerprints for a different snapshot.
