#!/usr/bin/env Rscript

# Harmonise stored pathway topology IDs with the bundled RaMP snapshot.
#
# Usage:
#   Rscript data-raw/update_ramp_graph_data.R /path/to/graph-era/source_df.rda
#   Rscript data-raw/update_ramp_graph_data.R /path/to/source_df.rda --dry-run

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
args <- setdiff(args, "--dry-run")
if (length(args) != 1L) {
  stop("Supply the source_df.rda used when the graph objects were built.")
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else {
  file.path(getwd(), "data-raw", "update_ramp_graph_data.R")
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
data_dir <- file.path(repo_root, "data")
source(file.path(repo_root, "data-raw", "ramp_graph_utils.R"))

previous_source_path <- normalizePath(args[[1]], mustWork = TRUE)
previous_source_df <- .rgu_load_object(previous_source_path, "source_df")
result <- .rgu_migrate_graph_data(
  data_dir, previous_source_df, write = !dry_run
)

if (!dry_run) {
  metadata_path <- file.path(data_dir, "ramp_db_metadata.rda")
  ramp_db_metadata <- .rgu_load_object(metadata_path, "ramp_db_metadata")
  ramp_db_metadata$graph_harmonisation <- c(
    list(
      ramp_version = ramp_db_metadata$ramp_version,
      harmonised_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      previous_source_file = basename(previous_source_path),
      strategy = paste(
        "maximum shared stable source identifiers; tied entity splits are",
        "expanded; two retired heme identifiers use curated replacements"
      )
    ),
    result$report
  )
  .rgu_save_object(ramp_db_metadata, "ramp_db_metadata", metadata_path)
}

cat(if (dry_run) "RaMP graph dry-run\n" else "RaMP graph update complete\n")
for (name in names(result$report)) {
  cat("  ", name, ": ", result$report[[name]], "\n", sep = "")
}
