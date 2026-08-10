# SpaMTP (developmental)

* Pathway enrichment, pathway-assay construction, Fisher analysis of m/z
  inputs, and interactive pathway networks now consume the scored `Ramp_IDs`
  produced by the indexed annotation engine. `AnnotateSM()` defaults to the
  bundled RaMP 3.0.7 chemical-property table, records annotation provenance in
  `@tools$mz_annotation`, and retains `@tools$db_3` only for compatibility.
  Legacy annotations require an explicit `annotation_source` fallback.
* Annotation candidates can now be stored once with `AnnotateSM(min_score =
  0)` and filtered later with `annotation_score_threshold`. Interactive
  pathway networks add `metabolite_detection = "annotated"` to display
  score-filtered pathway metabolites without requiring DE significance, while
  retaining separate leading-edge evidence in node styling and tooltips.
* Updated the bundled pruned RaMP snapshot from 2.5.4 to 3.0.7, with explicit
  upstream/source version metadata and a reproducible `data-raw` build script.
* Added a chemically validated adduct rule table and reusable sorted m/z index.
* Added ppm pruning, proton/charge bounds, mass-defect scoring, isotope checks,
  and contextual adduct-family scoring while retaining the legacy annotation
  result columns.
* `chem_props` can now be used directly as a RaMP-backed annotation database.
* Rebuilt `PathwayNetworkPlots()` around cached topology lookup, precomputed
  edges, selective sparse-matrix extraction, and JSON serialization. The new
  responsive viewer adds zoom/pan, focused/full networks, label controls,
  shared reaction markers, spatial inspection, and SVG export.
* Made legacy pathway annotation robust to unequal `Isomers_IDs` and
  `IsomerNames` counts, and added a cached Mouse Brain/DHB network demo.
* Pathway network metabolite labels now remain visible in the default label
  mode; dragged nodes stay pinned until released, and edge legends persist
  when the focused view changes.
* Added interactive repulsion, balanced-force, radial, and gene/metabolite
  pathway layouts, with a more widely spaced repulsion layout as the default.
* Harmonised all bundled pathway-graph node IDs with RaMP-DB 3.0.7, added a
  reproducible cross-version graph updater/auditor, and added stable source-ID
  labels when RaMP itself has no common name.

# SpaMTP 1.1.0

***Updated SpaMTP Release (Oct 2025)***

Additional Features:

* Updated `LoadSM()` function for compatability with **Cardinal V3.8** 
* Implementation of GraphPCA in R.
* Additional functions for handling large datasets - `AnnotateBigData` and `SelectROI`.
* Function for refining m/z annotations based on correlated pathway activity or MS/MS profiles.
* Package wide update for compatability with **Cardinal V3.8** and **Seurat V5.3**.


# SpaMTP 1.0.0

* Initial SpaMTP Release (March 2025)
