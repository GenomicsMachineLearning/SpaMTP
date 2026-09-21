#' Bundled annotation and pathway databases
#'
#' The standalone SpaMTP release includes the pruned RaMP-DB 3.0.7 snapshot,
#' corrected pathway graphs, legacy chemical reference tables, and precomputed
#' SMILES features. No companion package, Hub, or download is required.
#'
#' @format Data frames (`chem_props`, `source_df`, `analyte`,
#'   `analytehaspathway`, `pathway`, `smiles_features`, `HMDB_db`, `Chebi_db`,
#'   `Lipidmaps_db`, `GNPS_db`, and `filtered_fmp10`) or named lists
#'   (`ramp_db_metadata` and the four `RAMP_*` topology collections).
#'   [SpaMTPDatabaseInfo()] reports each dataset's dimensions and resource name.
#' @details
#' `chem_props` contains masses, formulae, structures and RaMP/source IDs.
#' `source_df` maps source identifiers to RaMP IDs; `analyte` records analyte
#' type. `analytehaspathway` and `pathway` describe pathway membership and names.
#' The `RAMP_*` lists contain KEGG, Reactome, WikiPathways and SMPDB graphs,
#' including restored `source_reaction_type` edge labels and repair provenance.
#' `ramp_db_metadata` records upstream releases and pruning information.
#' `smiles_features` stores structure-derived fields keyed by SMILES.
#'
#' Prefer [LoadSpaMTPDatabase()] for resource loading and provenance attributes.
#' Direct `utils::data()` access is retained for existing standalone workflows.
#' @source RaMP-DB: <https://github.com/ncats/RaMP-DB>.
#'   The bundled snapshot is shared with SpaMTP's main release. Precomputed
#'   structure features come from the immutable resource collection
#'   \doi{10.5281/zenodo.22045311}. Topology sources and checksums are recorded
#'   in `inst/extdata/pathway_interactions_provenance.csv`.
#' @seealso [LoadSpaMTPDatabase()], [SpaMTPDatabaseInfo()],
#'   [BuildMZAnnotationIndex()]
#' @name BundledDatabases
#' @aliases chem_props source_df analyte analytehaspathway pathway
#'   ramp_db_metadata smiles_features RAMP_hmdb RAMP_kegg RAMP_Reactome
#'   RAMP_wikipathway HMDB_db Chebi_db Lipidmaps_db GNPS_db filtered_fmp10
#' @docType data
#' @keywords datasets
NULL
