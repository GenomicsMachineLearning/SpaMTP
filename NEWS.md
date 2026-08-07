# SpaMTP (developmental)

* Updated the bundled pruned RaMP snapshot from 2.5.4 to 3.0.7, with explicit
  upstream/source version metadata and a reproducible `data-raw` build script.
* Added a chemically validated adduct rule table and reusable sorted m/z index.
* Added ppm pruning, proton/charge bounds, mass-defect scoring, isotope checks,
  and contextual adduct-family scoring while retaining the legacy annotation
  result columns.
* `chem_props` can now be used directly as a RaMP-backed annotation database.

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
