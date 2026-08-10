#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else {
  file.path(getwd(), "data-raw", "audit_ramp_graph_data.R")
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
data_dir <- file.path(repo_root, "data")
source(file.path(repo_root, "data-raw", "ramp_graph_utils.R"))

graphs <- stats::setNames(lapply(.rgu_graph_objects, function(object_name) {
  .rgu_load_object(file.path(data_dir, paste0(object_name, ".rda")), object_name)
}), .rgu_graph_objects)
analyte <- .rgu_load_object(file.path(data_dir, "analyte.rda"), "analyte")
source_df <- .rgu_load_object(file.path(data_dir, "source_df.rda"), "source_df")
chem_props <- .rgu_load_object(file.path(data_dir, "chem_props.rda"), "chem_props")
metadata <- .rgu_load_object(
  file.path(data_dir, "ramp_db_metadata.rda"), "ramp_db_metadata"
)

ids <- .rgu_graph_ids(graphs)
stale <- setdiff(ids, as.character(analyte$rampId))
names_audit <- .rgu_name_audit(ids, source_df, chem_props)
unnamed <- names_audit$graph_id[!names_audit$current_name_available]
timestamps <- unique(unlist(lapply(graphs, function(graph) {
  vapply(graph, function(topology) as.character(topology$timestamp), character(1))
}), use.names = FALSE))

cat("RaMP graph integrity audit\n")
cat("  bundled RaMP version: ", metadata$ramp_version, "\n", sep = "")
cat("  topology source timestamp(s): ", paste(timestamps, collapse = ", "), "\n", sep = "")
cat("  graph nodes: ", length(ids), "\n", sep = "")
cat("  compound nodes: ", sum(startsWith(ids, "RAMP_C_")), "\n", sep = "")
cat("  gene/protein nodes: ", sum(startsWith(ids, "RAMP_G_")), "\n", sep = "")
cat("  IDs absent from bundled analyte: ", length(stale), "\n", sep = "")
cat("  nodes without current names: ", length(unnamed), "\n", sep = "")
if (length(stale)) cat("  stale IDs: ", paste(stale, collapse = ", "), "\n", sep = "")
if (length(unnamed)) cat("  unnamed IDs: ", paste(unnamed, collapse = ", "), "\n", sep = "")

if (length(stale) || length(unnamed)) quit(status = 1L)
