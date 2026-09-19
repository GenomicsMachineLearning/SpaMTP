`cell_cycle_interactions.rds` is an offline regression fixture for KEGG
`hsa:04110` (Cell cycle), entry 116 in the pinned collections.

* `legacy_topology` is the unmodified entry from the published SpaMTPdb 3.0.7
  `ramp_kegg.rds`.
* `source_reaction_type` and `source_direction` are character labels from the
  corresponding graphite archive-19 `protEdges` table (2022-11-01).

Both input URLs and checksums are recorded in
`inst/extdata/pathway_interactions_provenance.csv`. The fixture has 1,009
protein-edge rows. It tests the shipped correction against source labels,
including directions and labels on parallel edges, without network access.
