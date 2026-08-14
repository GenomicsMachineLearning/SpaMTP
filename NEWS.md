# SpaMTP (developmental)

* Split versioned resources from the analysis code: `SpaMTPdb` now supplies
  the pruned RaMP annotation and pathway snapshot, and `SpaMTPData` provides
  named access to large experiment/vignette objects. `LoadSpaMTPDatabase()`
  and `SpaMTPDatabaseInfo()` form the public database interface. The software
  package retains only the small adduct and reaction-style constant tables.
* Updated the package citation to the peer-reviewed *Nature Methods* paper,
  designated Tianyao Lu as the current package maintainer, and added direct
  links to the developmental documentation and source branch. Andrew Causer
  remains credited as an original author and former maintainer.
* Added `ApplySpatialAlignment()` to apply existing SMINT coordinate outputs or
  homogeneous transforms, estimate landmark-based affine alignment in R, and
  run the SMINT-compatible STalign LDDMM workflow through an optional Python
  backend. Alignment provenance and nearest-target diagnostics are retained in
  the returned SpaMTP object without storing large Python velocity tensors.
* Pathway enrichment, pathway-assay construction, Fisher analysis of m/z
  inputs, and interactive pathway networks now consume the scored `Ramp_IDs`
  produced by the indexed annotation engine. `AnnotateSM()` resolves the
  versioned RaMP 3.0.7 chemical-property table through `SpaMTPdb`, records
  annotation provenance in `@tools$mz_annotation`, and retains `@tools$db_3`
  only for compatibility.
  Legacy annotations require an explicit `annotation_source` fallback.
* Annotation candidates can now be stored once with `AnnotateSM(min_score =
  0)` and filtered later with `annotation_score_threshold`. Interactive
  pathway networks add `metabolite_detection = "annotated"` to display
  score-filtered pathway metabolites without requiring DE significance, while
  retaining separate leading-edge evidence in node styling and tooltips.
* Moved the pruned RaMP snapshot from 2.5.4 to versioned SpaMTPdb 3.0.7
  resources, with explicit upstream/source metadata and reproducible resource
  staging scripts.
* Added a chemically validated adduct rule table and reusable sorted m/z index.
* Added ppm pruning, proton/charge bounds, mass-defect scoring, isotope checks,
  and contextual adduct-family scoring while retaining the legacy annotation
  result columns.
* Added dependency-free SMILES graph decomposition, functional-group and atom
  site reporting, mode-specific protonation/deprotonation and alkali-binding
  priors, and `PredictAdductsFromSMILES()`. Full RaMP runs join an independent
  precomputed `SpaMTPdb::smiles_features` resource; small custom databases are
  inferred automatically at runtime.
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
* Harmonised all versioned pathway-graph node IDs with RaMP-DB 3.0.7, added a
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
