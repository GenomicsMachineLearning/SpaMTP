#!/usr/bin/env Rscript
# Rscript data-raw/rebuild_pathway_interactions.R GRAPHITE_DIR OLD_SOURCE_RDA RESOURCE_DIR
# Inputs are pinned in inst/extdata/pathway_interactions_provenance.csv.
# The output is an in-memory correction for the immutable RaMP 3.0.7 graphs.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Supply GRAPHITE_DIR OLD_SOURCE_RDA RESOURCE_DIR.")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."))
source(file.path(root, "data-raw", "ramp_graph_utils.R"))
source(file.path(root, "data-raw", "pathway_interaction_utils.R"))
source(file.path(root, "R", "PathwayInteractionRepair.R"))
manifest <- read.csv(file.path(root, "inst", "extdata", "pathway_interactions_provenance.csv"))
check_input <- function(path, name = basename(path)) {
  expected <- manifest$md5[match(name, manifest$file)]
  if (is.na(expected) || !identical(unname(tools::md5sum(path)), expected)) {
    stop("Input checksum does not match the pinned snapshot: ", name)
  }
}
check_input(args[2L], "initial-source_df.rda")
for (name in c("source_df", "analyte", "ramp_kegg", "ramp_reactome", "ramp_wikipathway", "ramp_hmdb")) {
  check_input(file.path(args[3L], paste0(name, ".rds")))
}

old_source <- .rgu_load_object(args[2L], "source_df")
first <- !duplicated(old_source$sourceId) & !is.na(old_source$sourceId)
identifiers <- list2env(as.list(stats::setNames(
  as.character(old_source$rampId[first]), old_source$sourceId[first]
)), parent = emptyenv(), hash = TRUE)
current_source <- readRDS(file.path(args[3L], "source_df.rds"))
analyte <- readRDS(file.path(args[3L], "analyte.rds"))
databases <- c(ramp_kegg = "kegg", ramp_reactome = "reactome",
               ramp_wikipathway = "wikipathways", ramp_hmdb = "smpdb")
patches <- list()
audit <- list()
for (resource in names(databases)) {
  message("Restoring ", resource)
  raw_path <- file.path(args[1L], paste0("hsapiens-", databases[[resource]], "-19.rds"))
  check_input(raw_path)
  raw <- readRDS(raw_path)
  stopifnot(as.character(attr(raw, "timestamp")) == "2022-11-01")
  entries <- attr(raw, "entries")
  original <- readRDS(file.path(args[3L], paste0(resource, ".rds")))
  stopifnot(length(entries) == length(original))
  converted <- new.env(parent = emptyenv(), hash = TRUE)
  migrated <- new.env(parent = emptyenv(), hash = TRUE)
  edge_keys <- lapply(entries, function(p) {
    vapply(.rgu_edge_fields, function(field) {
      digest::digest(attr(p, field), algo = "sha256", serializeVersion = 2L)
    }, character(1))
  })
  convert_edge <- function(edge, key, legacy) {
    key <- paste0(legacy, key)
    if (!exists(key, converted, inherits = FALSE)) {
      converted[[key]] <- .pi_convert_edges(edge, identifiers, legacy)
    }
    converted[[key]]
  }
  # Recreate BOTH the buggy and corrected conversions, using identical ID
  # mappings and the same subsequent RaMP 3.0.7 identifier migration.
  convert <- function(legacy) lapply(seq_along(entries), function(i) {
    p <- entries[[i]]
    result <- attributes(p)[setdiff(names(attributes(p)), "class")]
    for (field in .rgu_edge_fields) {
      result[[field]] <- convert_edge(result[[field]], edge_keys[[i]][[field]], legacy)
    }
    result
  })
  broken <- convert(TRUE)
  migration <- .rgu_build_id_migration(
    .rgu_graph_ids(list(broken)), analyte, old_source, current_source
  )
  targets <- split(migration$new_id, migration$old_id)
  pathway_patches <- vector("list", length(entries))
  counts <- c(pathways = length(entries), tables = 0L, edges_before = 0L, edges_after = 0L)
  for (i in seq_along(entries)) {
    stopifnot(identical(broken[[i]]$id, original[[i]]$id))
    fixes <- list()
    for (field in .rgu_edge_fields) {
      key <- edge_keys[[i]][[field]]
      if (!exists(key, migrated, inherits = FALSE)) {
        restored <- convert_edge(attr(entries[[i]], field), key, FALSE)
        migrated[[key]] <- list(
          expected = .rgu_update_edge(broken[[i]][[field]], targets),
          corrected = .rgu_update_edge(restored, targets)
        )
      }
      expected <- migrated[[key]]$expected
      actual <- original[[i]][[field]]
      # Ignore empty-sentinel storage modes and row names, but require every
      # endpoint, direction and (corrupted) reaction code to match in order.
      equal <- vapply(names(actual), function(column) {
        identical(as.character(expected[[column]]), as.character(actual[[column]]))
      }, logical(1))
      if (!all(equal)) stop("Input reconstruction mismatch: ", resource, " ", i, " ", field,
                            " (", paste(names(equal)[!equal], collapse = ", "), ")")
      corrected <- migrated[[key]]$corrected
      fixes[[field]] <- .pi_make_patch(actual, corrected)
      counts["tables"] <- counts["tables"] + 1L
      counts["edges_before"] <- counts["edges_before"] + nrow(actual)
      counts["edges_after"] <- counts["edges_after"] + nrow(corrected)
    }
    pathway_patches[[i]] <- fixes
    if (i %% 5000L == 0L) message("  verified ", i, " / ", length(entries), " pathways")
  }
  patches[[resource]] <- list(
    input_sha256 = .spamtp_graph_fingerprint(original),
    provenance = list(repair = "pathway-interactions-v1", graphite_archive = 19L,
                       graphite_date = "2022-11-01", ramp_version = "3.0.7",
                       source_md5 = unname(tools::md5sum(raw_path))),
    pathways = pathway_patches
  )
  audit[[resource]] <- data.frame(resource = resource, as.list(counts))
  print(audit[[resource]])
  rm(raw, entries, original, broken, pathway_patches, converted, migrated, edge_keys)
  gc(FALSE)
}
output <- file.path(root, "inst", "extdata")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
for (resource in names(patches)) {
  saveRDS(patches[[resource]], file.path(output, paste0(resource, "_interactions_v1.rds")),
          compress = "xz", version = 2)
}
write.csv(do.call(rbind, audit), file.path(output, "pathway_interactions_audit.csv"), row.names = FALSE)
# Code 16 is a generic process, not an inferred indirect effect.
reaction_type <- .rgu_load_object(file.path(root, "data", "reaction_type.rda"), "reaction_type")
reaction_type$reaction_name[reaction_type$reaction_type == 16L] <- "Process"
reaction_type$linetype[reaction_type$reaction_type == 16L] <- "solid"
.rgu_save_object(reaction_type, "reaction_type", file.path(root, "data", "reaction_type.rda"))
