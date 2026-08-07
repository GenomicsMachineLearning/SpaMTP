#!/usr/bin/env Rscript

# Rebuild the pruned RaMP snapshot distributed with SpaMTP.
#
# Usage:
#   Rscript data-raw/update_ramp_data.R /path/to/RaMP_SQLite_v3.0.7.sqlite

required_packages <- c("DBI", "RSQLite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install build dependency/dependencies: ",
       paste(missing_packages, collapse = ", "))
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop("Supply the path to a decompressed RaMP SQLite database.")
}
sqlite_path <- normalizePath(args[[1]], mustWork = TRUE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else {
  file.path(getwd(), "data-raw", "update_ramp_data.R")
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  stop("Could not locate the SpaMTP package root.")
}
data_dir <- file.path(repo_root, "data")

con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
on.exit(DBI::dbDisconnect(con), add = TRUE)

required_tables <- c(
  "db_version", "version_info", "analyte", "source", "analytehaspathway",
  "pathway", "chem_props"
)
missing_tables <- setdiff(required_tables, DBI::dbListTables(con))
if (length(missing_tables)) {
  stop("RaMP database is missing table(s): ",
       paste(missing_tables, collapse = ", "))
}

db_version <- DBI::dbGetQuery(con, "SELECT * FROM db_version")
if (nrow(db_version) != 1 || !grepl("^3\\.", db_version$ramp_version[[1]])) {
  stop("This updater requires a RaMP 3.x database.")
}
source_versions <- DBI::dbGetQuery(
  con,
  paste(
    "SELECT data_source_id, data_source_name, data_source_url,",
    "data_source_version FROM version_info WHERE status = 'current'",
    "ORDER BY data_source_id"
  )
)

queries <- list(
  analyte = "SELECT rampId, type FROM analyte",
  source_df = paste(
    "SELECT sourceId, rampId, IDtype, geneOrCompound, commonName,",
    "priorityHMDBStatus, dataSource, pathwayCount FROM source"
  ),
  analytehaspathway = paste(
    "SELECT rampId, pathwayRampId, pathwaySource FROM analytehaspathway"
  ),
  pathway = paste(
    "SELECT pathwayRampId, sourceId, type, pathwayCategory, pathwayName",
    "FROM pathway"
  ),
  chem_props = paste(
    "SELECT ramp_id, chem_data_source, chem_source_id, iso_smiles,",
    "inchi_key_prefix, inchi_key, inchi, mw, monoisotop_mass,",
    "common_name, mol_formula FROM chem_props"
  )
)

save_snapshot_object <- function(object_name, sql) {
  message("Extracting ", object_name, " ...")
  value <- DBI::dbGetQuery(con, sql)
  if (!nrow(value)) stop("Extracted table is empty: ", object_name)
  environment <- new.env(parent = emptyenv())
  environment[[object_name]] <- value
  save(
    list = object_name,
    envir = environment,
    file = file.path(data_dir, paste0(object_name, ".rda")),
    compress = "xz",
    version = 3
  )
  rows <- nrow(value)
  rm(value, environment)
  gc(verbose = FALSE)
  rows
}

snapshot_counts <- vapply(
  names(queries),
  function(object_name) save_snapshot_object(object_name, queries[[object_name]]),
  numeric(1)
)

ramp_db_metadata <- list(
  ramp_version = db_version$ramp_version[[1]],
  load_timestamp = as.character(db_version$load_timestamp[[1]]),
  version_notes = db_version$version_notes[[1]],
  source_versions = source_versions,
  snapshot_counts = snapshot_counts,
  upstream_repository = "https://github.com/ncats/RaMP-DB",
  upstream_database_file = basename(sqlite_path),
  pruning = list(
    analyte = c("rampId", "type"),
    source_df = c(
      "sourceId", "rampId", "IDtype", "geneOrCompound", "commonName",
      "priorityHMDBStatus", "dataSource", "pathwayCount"
    ),
    analytehaspathway = c("rampId", "pathwayRampId", "pathwaySource"),
    pathway = c(
      "pathwayRampId", "sourceId", "type", "pathwayCategory", "pathwayName"
    ),
    chem_props = c(
      "ramp_id", "chem_data_source", "chem_source_id", "iso_smiles",
      "inchi_key_prefix", "inchi_key", "inchi", "mw", "monoisotop_mass",
      "common_name", "mol_formula"
    )
  )
)
save(
  ramp_db_metadata,
  file = file.path(data_dir, "ramp_db_metadata.rda"),
  compress = "xz",
  version = 3
)

message(
  "SpaMTP RaMP snapshot updated to ", ramp_db_metadata$ramp_version,
  ": ", paste(names(snapshot_counts), snapshot_counts, sep = "=", collapse = ", ")
)
