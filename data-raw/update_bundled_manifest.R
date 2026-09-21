# Run from the SpaMTP source root after updating the bundled datasets.
# Only base R is needed; no companion package or network is consulted.
source("R/DatabaseResources.R", local = TRUE)
metadata <- new.env(parent = emptyenv())
load("data/ramp_db_metadata.rda", envir = metadata)
version <- metadata$ramp_db_metadata$ramp_version
manifest <- lapply(names(.spamtp_db_legacy_names), function(resource) {
  object <- unname(.spamtp_db_legacy_names[[resource]])
  path <- file.path("data", paste0(object, ".rda"))
  environment <- new.env(parent = emptyenv())
  load(path, envir = environment)
  value <- get(object, envir = environment, inherits = FALSE)
  data.frame(
    resource = resource, version = version, source = "bundled",
    category = if (resource == "smiles_features") "structure" else if (
      startsWith(resource, "ramp_") && resource != "ramp_db_metadata"
    ) "topology" else if (resource %in% c(
      "hmdb_db", "chebi_db", "lipidmaps_db", "gnps_db", "filtered_fmp10"
    )) "legacy" else "core",
    source_object = object, rdata_path = path,
    rows = if (is.data.frame(value)) nrow(value) else length(value),
    columns = if (is.data.frame(value)) ncol(value) else NA_integer_,
    serialized_bytes = unname(file.info(path)$size),
    md5 = unname(tools::md5sum(path)), stringsAsFactors = FALSE
  )
})
utils::write.csv(do.call(rbind, manifest),
                 "inst/extdata/bundled_database_manifest.csv", row.names = FALSE)
