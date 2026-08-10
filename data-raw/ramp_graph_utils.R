.rgu_graph_objects <- c(
  "RAMP_hmdb", "RAMP_kegg", "RAMP_wikipathway", "RAMP_Reactome"
)

.rgu_edge_fields <- c(
  "protEdges", "protPropEdges", "metabolEdges", "metabolPropEdges",
  "mixedEdges"
)

.rgu_assert_columns <- function(value, columns, label) {
  missing <- setdiff(columns, names(value))
  if (length(missing)) {
    stop(label, " is missing column(s): ", paste(missing, collapse = ", "))
  }
}

.rgu_load_object <- function(path, object_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  if (!object_name %in% loaded) {
    stop(path, " does not contain object '", object_name, "'.")
  }
  environment[[object_name]]
}

.rgu_save_object <- function(value, object_name, path) {
  environment <- new.env(parent = emptyenv())
  environment[[object_name]] <- value
  temporary <- tempfile(pattern = paste0(object_name, "_"), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  save(
    list = object_name, envir = environment, file = temporary,
    compress = "xz", version = 3
  )
  if (!file.copy(temporary, path, overwrite = TRUE)) {
    stop("Failed to replace ", path)
  }
  invisible(path)
}

.rgu_topology_ids <- function(topology) {
  unique(unlist(lapply(intersect(.rgu_edge_fields, names(topology)), function(field) {
    edge <- topology[[field]]
    if (!is.data.frame(edge) || !nrow(edge)) return(character())
    c(as.character(edge$src), as.character(edge$dest))
  }), use.names = FALSE))
}

.rgu_graph_ids <- function(graphs) {
  value <- unique(unlist(lapply(graphs, function(graph) {
    unlist(lapply(graph, .rgu_topology_ids), use.names = FALSE)
  }), use.names = FALSE))
  value[!is.na(value) & nzchar(value)]
}

.rgu_default_overrides <- function() {
  data.frame(
    old_id = c("RAMP_C_000218485", "RAMP_C_000218486"),
    new_id = c("RAMP_C_000218291", "RAMP_C_000218617"),
    evidence = c(
      "ferroheme replaced by current RaMP ferroheme b(2-) entity",
      "ferriheme replaced by current RaMP ferriheme b entity"
    ),
    stringsAsFactors = FALSE
  )
}

.rgu_build_id_migration <- function(graph_ids, analyte, previous_source,
                                    current_source,
                                    manual_overrides = .rgu_default_overrides()) {
  .rgu_assert_columns(analyte, "rampId", "analyte")
  required <- c("sourceId", "rampId", "geneOrCompound", "commonName")
  .rgu_assert_columns(previous_source, required, "previous source_df")
  .rgu_assert_columns(current_source, required, "current source_df")

  current_ids <- unique(as.character(analyte$rampId))
  stale_ids <- setdiff(graph_ids, current_ids)
  if (!length(stale_ids)) {
    return(data.frame(
      old_id = character(), new_id = character(), shared_ids = integer(),
      method = character(), evidence = character(), split = logical(),
      stringsAsFactors = FALSE
    ))
  }

  old <- previous_source[
    previous_source$rampId %in% stale_ids &
      !is.na(previous_source$sourceId) & nzchar(previous_source$sourceId),
    c("rampId", "sourceId"), drop = FALSE
  ]
  current <- current_source[
    !is.na(current_source$sourceId) & nzchar(current_source$sourceId),
    c("rampId", "sourceId"), drop = FALSE
  ]
  old$source_key <- tolower(trimws(as.character(old$sourceId)))
  current$source_key <- tolower(trimws(as.character(current$sourceId)))
  old <- unique(old[c("rampId", "source_key")])
  current <- unique(current[c("rampId", "source_key")])

  candidates <- merge(
    old, current, by = "source_key", suffixes = c("_old", "_new"),
    sort = FALSE
  )
  candidates <- candidates[
    substr(candidates$rampId_old, 1L, 7L) ==
      substr(candidates$rampId_new, 1L, 7L),
    , drop = FALSE
  ]
  candidates <- unique(candidates[c("rampId_old", "rampId_new", "source_key")])

  migration <- candidates[0, c("rampId_old", "rampId_new"), drop = FALSE]
  migration$shared_ids <- integer()
  if (nrow(candidates)) {
    votes <- stats::aggregate(
      source_key ~ rampId_old + rampId_new, candidates, length
    )
    names(votes)[names(votes) == "source_key"] <- "shared_ids"
    maximum <- stats::aggregate(shared_ids ~ rampId_old, votes, max)
    migration <- merge(
      votes, maximum, by = c("rampId_old", "shared_ids"), sort = FALSE
    )
  }
  names(migration)[names(migration) == "rampId_old"] <- "old_id"
  names(migration)[names(migration) == "rampId_new"] <- "new_id"
  migration$method <- "shared_source_id"
  migration$evidence <- paste0(migration$shared_ids, " shared source ID(s)")

  missing <- setdiff(stale_ids, migration$old_id)
  override <- manual_overrides[
    manual_overrides$old_id %in% missing &
      manual_overrides$new_id %in% current_ids,
    , drop = FALSE
  ]
  if (nrow(override)) {
    override$shared_ids <- 0L
    override$method <- "curated_override"
    override <- override[
      c("old_id", "new_id", "shared_ids", "method", "evidence")
    ]
    migration <- rbind(
      migration[c("old_id", "new_id", "shared_ids", "method", "evidence")],
      override
    )
  }

  unresolved <- setdiff(stale_ids, migration$old_id)
  if (length(unresolved)) {
    stop(
      "No current RaMP ID could be resolved for graph node(s): ",
      paste(unresolved, collapse = ", ")
    )
  }
  target_count <- table(migration$old_id)
  migration$split <- unname(target_count[migration$old_id]) > 1L
  migration <- migration[order(migration$old_id, migration$new_id), , drop = FALSE]
  rownames(migration) <- NULL
  migration
}

.rgu_update_edge <- function(edge, target_map) {
  if (!is.data.frame(edge) || !nrow(edge)) return(edge)
  .rgu_assert_columns(edge, c("src", "dest"), "graph edge table")
  mapped_ids <- names(target_map)
  if (!any(edge$src %in% mapped_ids) && !any(edge$dest %in% mapped_ids)) {
    return(edge)
  }
  original_names <- names(edge)
  edge$src <- as.character(edge$src)
  edge$dest <- as.character(edge$dest)
  split_ids <- names(target_map)[lengths(target_map) > 1L]
  expanded_row <- edge$src %in% split_ids | edge$dest %in% split_ids

  simple <- edge[!expanded_row, , drop = FALSE]
  if (nrow(simple)) {
    src_index <- match(simple$src, names(target_map))
    dest_index <- match(simple$dest, names(target_map))
    replace_src <- !is.na(src_index)
    replace_dest <- !is.na(dest_index)
    simple$src[replace_src] <- vapply(
      target_map[src_index[replace_src]], `[[`, character(1), 1L
    )
    simple$dest[replace_dest] <- vapply(
      target_map[dest_index[replace_dest]], `[[`, character(1), 1L
    )
  }

  expanded <- lapply(which(expanded_row), function(index) {
    source <- target_map[[edge$src[[index]]]]
    destination <- target_map[[edge$dest[[index]]]]
    if (is.null(source)) source <- edge$src[[index]]
    if (is.null(destination)) destination <- edge$dest[[index]]
    combinations <- expand.grid(
      src = source, dest = destination, stringsAsFactors = FALSE
    )
    value <- edge[rep(index, nrow(combinations)), , drop = FALSE]
    value$src <- combinations$src
    value$dest <- combinations$dest
    value
  })
  result <- if (length(expanded)) {
    rbind(simple, do.call(rbind, expanded))
  } else {
    simple
  }
  result <- unique(result[original_names])
  rownames(result) <- NULL
  result
}

.rgu_update_graphs <- function(graphs, migration) {
  if (!nrow(migration)) return(graphs)
  target_map <- split(migration$new_id, migration$old_id)
  lapply(graphs, function(graph) {
    lapply(graph, function(topology) {
      for (field in intersect(.rgu_edge_fields, names(topology))) {
        topology[[field]] <- .rgu_update_edge(topology[[field]], target_map)
      }
      topology
    })
  })
}

.rgu_name_audit <- function(graph_ids, source_df, chem_props) {
  source_names <- unique(as.character(source_df$rampId[
    !is.na(source_df$commonName) & nzchar(trimws(source_df$commonName))
  ]))
  chemical_names <- unique(as.character(chem_props$ramp_id[
    !is.na(chem_props$common_name) & nzchar(trimws(chem_props$common_name))
  ]))
  source_identifiers <- unique(as.character(source_df$rampId[
    !is.na(source_df$sourceId) & nzchar(trimws(source_df$sourceId))
  ]))
  data.frame(
    graph_id = graph_ids,
    current_name_available = graph_ids %in%
      union(source_identifiers, union(source_names, chemical_names)),
    stringsAsFactors = FALSE
  )
}

.rgu_migrate_graph_data <- function(data_dir, previous_source_df,
                                    write = FALSE,
                                    manual_overrides = .rgu_default_overrides()) {
  graph_paths <- file.path(data_dir, paste0(.rgu_graph_objects, ".rda"))
  if (any(!file.exists(graph_paths))) {
    stop("Missing graph data object(s): ",
         paste(basename(graph_paths[!file.exists(graph_paths)]), collapse = ", "))
  }
  graphs <- stats::setNames(lapply(seq_along(graph_paths), function(index) {
    .rgu_load_object(graph_paths[[index]], .rgu_graph_objects[[index]])
  }), .rgu_graph_objects)
  analyte <- .rgu_load_object(file.path(data_dir, "analyte.rda"), "analyte")
  source_df <- .rgu_load_object(file.path(data_dir, "source_df.rda"), "source_df")
  chem_props <- .rgu_load_object(file.path(data_dir, "chem_props.rda"), "chem_props")

  ids_before <- .rgu_graph_ids(graphs)
  migration <- .rgu_build_id_migration(
    ids_before, analyte, previous_source_df, source_df, manual_overrides
  )
  updated <- .rgu_update_graphs(graphs, migration)
  ids_after <- .rgu_graph_ids(updated)
  stale_after <- setdiff(ids_after, as.character(analyte$rampId))
  if (length(stale_after)) {
    stop("Graph migration left stale RaMP IDs: ", paste(stale_after, collapse = ", "))
  }
  name_audit <- .rgu_name_audit(ids_after, source_df, chem_props)
  unnamed <- name_audit$graph_id[!name_audit$current_name_available]
  if (length(unnamed)) {
    stop("Current graph node(s) still have no name: ", paste(unnamed, collapse = ", "))
  }

  edges_before <- sum(vapply(graphs, function(graph) {
    sum(vapply(graph, function(topology) {
      sum(vapply(intersect(.rgu_edge_fields, names(topology)), function(field) {
        nrow(topology[[field]])
      }, integer(1)))
    }, integer(1)))
  }, integer(1)))
  edges_after <- sum(vapply(updated, function(graph) {
    sum(vapply(graph, function(topology) {
      sum(vapply(intersect(.rgu_edge_fields, names(topology)), function(field) {
        nrow(topology[[field]])
      }, integer(1)))
    }, integer(1)))
  }, integer(1)))

  report <- list(
    graph_nodes_before = length(ids_before),
    graph_nodes_after = length(ids_after),
    stale_ids_migrated = length(unique(migration$old_id)),
    current_targets = length(unique(migration$new_id)),
    split_ids_expanded = length(unique(migration$old_id[migration$split])),
    curated_overrides = sum(migration$method == "curated_override"),
    edges_before = edges_before,
    edges_after = edges_after,
    stale_ids_remaining = length(stale_after),
    unnamed_nodes_remaining = length(unnamed)
  )

  if (isTRUE(write)) {
    for (object_name in .rgu_graph_objects) {
      .rgu_save_object(
        updated[[object_name]], object_name,
        file.path(data_dir, paste0(object_name, ".rda"))
      )
    }
  }
  list(graphs = updated, migration = migration, name_audit = name_audit,
       report = report)
}
