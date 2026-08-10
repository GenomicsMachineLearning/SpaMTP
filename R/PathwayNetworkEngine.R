.spamtp_pathway_cache <- new.env(parent = emptyenv())

.pn_is_current_annotation_metadata <- function(metadata) {
  current_engines <- c(.spamtp_annotation_engine, "curated-fmp10-v2")
  is.list(metadata) && length(metadata$engine) == 1L &&
    metadata$engine %in% current_engines
}

.pn_assert_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop(label, " is missing required column(s): ", paste(missing, collapse = ", "))
  }
  invisible(x)
}

.pn_type_key <- function(x) {
  value <- tolower(trimws(as.character(x %||% "")))
  if (grepl("wiki", value, fixed = TRUE)) return("wiki")
  if (grepl("reactome", value, fixed = TRUE)) return("reactome")
  if (grepl("kegg", value, fixed = TRUE)) return("kegg")
  if (grepl("hmdb", value, fixed = TRUE)) return("hmdb")
  NA_character_
}

.pn_topology_signature <- function(topology_sources) {
  paste(
    names(topology_sources),
    lengths(topology_sources),
    vapply(topology_sources, function(x) {
      if (!length(x)) return("")
      paste(x[[1]]$id %||% "", x[[length(x)]]$id %||% "", sep = ":")
    }, character(1)),
    collapse = "|"
  )
}

.pn_topology_catalog <- function(topology_sources, use_cache = TRUE) {
  if (is.null(names(topology_sources)) || any(!nzchar(names(topology_sources)))) {
    stop("topology_sources must be a named list.")
  }
  signature <- .pn_topology_signature(topology_sources)
  if (isTRUE(use_cache) &&
      identical(.spamtp_pathway_cache$signature, signature) &&
      is.data.frame(.spamtp_pathway_cache$catalog)) {
    return(.spamtp_pathway_cache$catalog)
  }

  parts <- Map(function(collection, source) {
    if (!length(collection)) return(NULL)
    titles <- names(collection)
    if (is.null(titles)) {
      titles <- vapply(collection, function(x) as.character(x$title %||% ""), character(1))
    }
    ids <- vapply(collection, function(x) as.character(x$id %||% ""), character(1))
    data.frame(
      source = source,
      index = seq_along(collection),
      name_key = tolower(trimws(titles)),
      id_key = tolower(trimws(ids)),
      stringsAsFactors = FALSE
    )
  }, topology_sources, names(topology_sources))
  catalog <- do.call(rbind, parts)
  rownames(catalog) <- NULL
  if (isTRUE(use_cache)) {
    .spamtp_pathway_cache$signature <- signature
    .spamtp_pathway_cache$catalog <- catalog
  }
  catalog
}

.pn_select_pathway_rows <- function(regpathway, selected_pathways = NULL,
                                    top_n_pathways = 10L) {
  .pn_assert_columns(
    regpathway,
    c("pathwayName", "sourceId", "type", "Cluster_id", "leadingEdge"),
    "regpathway"
  )
  if (!nrow(regpathway)) stop("regpathway contains no rows.")

  name_key <- tolower(trimws(as.character(regpathway$pathwayName)))
  id_key <- tolower(trimws(as.character(regpathway$sourceId)))
  if (is.null(selected_pathways)) {
    .pn_assert_columns(regpathway, "NES", "regpathway")
    score <- suppressWarnings(abs(as.numeric(regpathway$NES)))
    score[!is.finite(score)] <- 0
    totals <- rowsum(score, group = name_key, reorder = FALSE, na.rm = TRUE)
    ranked <- rownames(totals)[order(totals[, 1], decreasing = TRUE)]
    selected_key <- utils::head(ranked, max(1L, as.integer(top_n_pathways)))
    keep <- name_key %in% selected_key
  } else {
    selected_key <- tolower(trimws(as.character(selected_pathways)))
    selected_key <- selected_key[nzchar(selected_key)]
    if (!length(selected_key)) stop("selected_pathways contains no usable values.")
    keep <- name_key %in% selected_key | id_key %in% selected_key
  }
  result <- regpathway[keep, , drop = FALSE]
  if (!nrow(result)) {
    stop("None of selected_pathways matched pathwayName or sourceId in regpathway.")
  }
  result
}

.pn_resolve_topologies <- function(pathway_rows, catalog, topology_sources) {
  unique_rows <- pathway_rows[!duplicated(tolower(pathway_rows$pathwayName)), , drop = FALSE]
  resolved <- vector("list", nrow(unique_rows))
  found <- logical(nrow(unique_rows))
  for (i in seq_len(nrow(unique_rows))) {
    source <- .pn_type_key(unique_rows$type[[i]])
    candidates <- catalog
    if (!is.na(source)) {
      candidates <- candidates[candidates$source == source, , drop = FALSE]
    }
    name_key <- tolower(trimws(as.character(unique_rows$pathwayName[[i]])))
    id_key <- tolower(trimws(as.character(unique_rows$sourceId[[i]])))
    hit <- which(candidates$name_key == name_key | candidates$id_key == id_key)
    if (length(hit)) {
      row <- candidates[hit[[1]], , drop = FALSE]
      resolved[[i]] <- topology_sources[[row$source]][[row$index]]
      found[[i]] <- TRUE
    }
  }
  names(resolved) <- as.character(unique_rows$pathwayName)
  list(topologies = resolved, found = found, pathways = unique_rows)
}

.pn_prepare_edges <- function(topology, reaction_types) {
  edge_sets <- c(
    "mixedEdges", "metabolEdges", "metabolPropEdges",
    "protEdges", "protPropEdges"
  )
  parts <- lapply(intersect(edge_sets, names(topology)), function(name) {
    edge <- topology[[name]]
    if (!is.data.frame(edge) || !nrow(edge) ||
        !all(c("src", "dest", "reaction_type") %in% names(edge))) {
      return(NULL)
    }
    edge[, intersect(c("src", "dest", "directed", "reaction_type"), names(edge)), drop = FALSE]
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(.pn_empty_links())

  edges <- unique(do.call(rbind, parts))
  edges$src <- as.character(edges$src)
  edges$dest <- as.character(edges$dest)
  keep <- !is.na(edges$src) & nzchar(edges$src) &
    !is.na(edges$dest) & nzchar(edges$dest)
  edges <- edges[keep, , drop = FALSE]
  if (!nrow(edges)) return(.pn_empty_links())

  reaction_index <- match(
    as.character(edges$reaction_type),
    as.character(reaction_types$reaction_type)
  )
  lookup <- function(column, default) {
    value <- as.character(reaction_types[[column]][reaction_index])
    value[is.na(value) | !nzchar(value)] <- default
    value
  }
  edges$reaction_name <- lookup("reaction_name", "Interaction")
  edges$linetype <- lookup("linetype", "solid")
  edges$arrowhead <- lookup("arrowhead", "arrow")
  edges$colour <- lookup("colour", "#64748b")

  degree <- table(c(edges$src, edges$dest))
  edges$weight <- as.integer(degree[edges$src]) + as.integer(degree[edges$dest])
  rownames(edges) <- NULL
  edges
}

.pn_empty_nodes <- function() {
  data.frame(
    id = character(), group = character(), expr = numeric(),
    name = character(), display = character(), detected = logical(),
    leading_edge = logical(), annotated = logical(),
    annotation_score = numeric(),
    stringsAsFactors = FALSE
  )
}

.pn_empty_links <- function() {
  data.frame(
    source = character(), target = character(), colour = character(),
    weight = integer(), style = character(), head = character(),
    reaction = character(), stringsAsFactors = FALSE
  )
}

.pn_prune_edges <- function(edges, detected, max_nodes) {
  if (!nrow(edges)) return(edges)
  nodes <- unique(c(edges$src, edges$dest))
  if (!is.finite(max_nodes)) return(edges)
  max_nodes <- as.integer(max_nodes)
  if (max_nodes <= 0L || length(nodes) <= max_nodes) {
    return(edges)
  }
  detected <- intersect(unique(as.character(detected)), nodes)
  incident <- edges$src %in% detected | edges$dest %in% detected
  neighbours <- unique(c(edges$src[incident], edges$dest[incident]))
  degree <- sort(table(c(edges$src, edges$dest)), decreasing = TRUE)
  cap <- max(max_nodes, length(detected))
  keep_nodes <- unique(c(detected, neighbours, names(degree)))
  keep_nodes <- utils::head(keep_nodes, cap)
  edges[edges$src %in% keep_nodes & edges$dest %in% keep_nodes, , drop = FALSE]
}

.pn_expand_values <- function(df, column, separator = ";") {
  values <- strsplit(as.character(df[[column]]), separator, fixed = TRUE)
  lengths <- lengths(values)
  keep <- lengths > 0L
  if (!any(keep)) return(df[0, , drop = FALSE])
  expanded <- df[rep(which(keep), lengths[keep]), , drop = FALSE]
  expanded[[column]] <- trimws(unlist(values[keep], use.names = FALSE))
  expanded[nzchar(expanded[[column]]) & !is.na(expanded[[column]]), , drop = FALSE]
}

.pn_normalise_source_id <- function(x) {
  x <- trimws(as.character(x))
  has_prefix <- grepl(":", x, fixed = TRUE)
  source <- ifelse(has_prefix, tolower(sub(":.*$", "", x)), "")
  value <- ifelse(has_prefix, sub("^[^:]+:", "", x), x)
  is_hmdb <- grepl("^HMDB", value, ignore.case = TRUE) | source == "hmdb"
  is_chebi <- grepl("^CHEBI", value, ignore.case = TRUE) | source == "chebi"
  is_lipid <- grepl("^(LM|LMPK)", value, ignore.case = TRUE) |
    source %in% c("lipidmaps", "lipidmap")
  value[is_chebi] <- sub("^CHEBI:", "", value[is_chebi], ignore.case = TRUE)
  result <- x
  result[is_hmdb] <- paste0("hmdb:", toupper(value[is_hmdb]))
  result[is_chebi] <- paste0("chebi:", value[is_chebi])
  result[is_lipid] <- paste0("lipidmaps:", value[is_lipid])
  result
}

.pn_normalise_de_list <- function(DE.list, analyte_types) {
  if (!is.list(DE.list) || length(DE.list) != length(analyte_types)) {
    stop("DE.list must contain exactly one data.frame for each analyte type.")
  }
  if (all(analyte_types %in% names(DE.list))) {
    result <- DE.list[analyte_types]
  } else {
    result <- stats::setNames(DE.list, analyte_types)
  }
  for (type in analyte_types) {
    de <- as.data.frame(result[[type]], stringsAsFactors = FALSE)
    if ("logFC" %in% names(de) && !"avg_log2FC" %in% names(de)) {
      names(de)[names(de) == "logFC"] <- "avg_log2FC"
    }
    if ("FDR" %in% names(de) && !"p_val_adj" %in% names(de)) {
      names(de)[names(de) == "FDR"] <- "p_val_adj"
    }
    .pn_assert_columns(de, c("cluster", "gene", "avg_log2FC", "p_val_adj"),
                       paste0("DE.list[['", type, "']]") )
    de$cluster <- as.character(de$cluster)
    de$gene <- as.character(de$gene)
    de$avg_log2FC <- suppressWarnings(as.numeric(de$avg_log2FC))
    de$p_val_adj <- suppressWarnings(as.numeric(de$p_val_adj))
    de$p_val_adj[!is.finite(de$p_val_adj)] <- Inf
    de$.cluster_key <- tolower(de$cluster)
    result[[type]] <- de
  }
  result
}

.pn_prepare_gene_de <- function(de, source_table) {
  direct <- startsWith(toupper(de$gene), "RAMP_G_")
  direct_de <- de[direct, , drop = FALSE]
  direct_de$rampId <- toupper(direct_de$gene)

  symbol_de <- de[!direct, , drop = FALSE]
  mapped <- symbol_de[0, , drop = FALSE]
  mapped$rampId <- character()
  if (nrow(symbol_de)) {
    symbols <- unique(toupper(symbol_de$gene))
    gene_source <- source_table[
      startsWith(source_table$rampId, "RAMP_G_") &
        toupper(source_table$commonName) %in% symbols,
      c("rampId", "commonName"), drop = FALSE
    ]
    gene_source$gene_key <- toupper(gene_source$commonName)
    gene_source <- unique(gene_source[c("rampId", "gene_key")])
    symbol_de$gene_key <- toupper(symbol_de$gene)
    mapped <- merge(symbol_de, gene_source, by = "gene_key", sort = FALSE)
    mapped$gene_key <- NULL
  }
  result <- rbind(direct_de, mapped[, names(direct_de), drop = FALSE])
  if (!nrow(result)) return(result)
  result$commonName <- toupper(result$gene)
  result
}

.pn_prepare_metabolite_map <- function(db_3, chemical_properties = NULL) {
  db_3 <- as.data.frame(db_3, stringsAsFactors = FALSE)
  if (!"mz_name" %in% names(db_3)) {
    .pn_assert_columns(db_3, "observed_mz", "SpaMTP@tools$db_3")
    db_3$mz_name <- paste0("mz-", db_3$observed_mz)
  }
  adduct <- if ("Adduct" %in% names(db_3)) as.character(db_3$Adduct) else ""
  db_3$Adduct <- adduct

  if ("ramp_id" %in% names(db_3) &&
      any(grepl("^RAMP_C_", as.character(db_3$ramp_id)))) {
    mapping <- db_3
    mapping$ramp_id <- toupper(as.character(mapping$ramp_id))
    mapping$common_name <- if ("IsomerNames" %in% names(mapping)) {
      as.character(mapping$IsomerNames)
    } else if ("common_name" %in% names(mapping)) {
      as.character(mapping$common_name)
    } else {
      NA_character_
    }
  } else if ("Ramp_IDs" %in% names(db_3) && any(nzchar(stats::na.omit(db_3$Ramp_IDs)))) {
    mapping <- .pn_expand_values(db_3, "Ramp_IDs")
    mapping$ramp_id <- mapping$Ramp_IDs
    mapping$common_name <- if ("IsomerNames" %in% names(mapping)) {
      as.character(mapping$IsomerNames)
    } else {
      NA_character_
    }
  } else {
    if (is.null(chemical_properties)) {
      stop("Legacy metabolite annotations require the chem_props lookup table.")
    }
    id_column <- if ("Isomers_IDs" %in% names(db_3)) "Isomers_IDs" else "Isomers"
    .pn_assert_columns(db_3, id_column, "SpaMTP@tools$db_3")
    mapping <- .pn_expand_values(db_3, id_column)
    mapping$chem_source_id <- .pn_normalise_source_id(mapping[[id_column]])
    chem_map <- unique(chemical_properties[c("chem_source_id", "ramp_id")])
    mapping <- merge(mapping, chem_map, by = "chem_source_id", sort = FALSE)
  }
  mapping <- mapping[startsWith(mapping$ramp_id, "RAMP_C_"), , drop = FALSE]
  if (!nrow(mapping)) {
    return(data.frame(
      ramp_id = character(), mz_name = character(), Adduct = character(),
      common_name = character(), stringsAsFactors = FALSE
    ))
  }
  if (!"common_name" %in% names(mapping)) mapping$common_name <- NA_character_
  if (!is.null(chemical_properties)) {
    name_map <- chemical_properties[
      chemical_properties$ramp_id %in% unique(mapping$ramp_id),
      c("ramp_id", "common_name"), drop = FALSE
    ]
    name_map <- name_map[!duplicated(name_map$ramp_id), , drop = FALSE]
    names(name_map)[names(name_map) == "common_name"] <- ".ramp_common_name"
    mapping <- merge(mapping, name_map, by = "ramp_id", all.x = TRUE, sort = FALSE)
    missing_name <- is.na(mapping$common_name) | !nzchar(mapping$common_name)
    mapping$common_name[missing_name] <- mapping$.ramp_common_name[missing_name]
    mapping$.ramp_common_name <- NULL
  }
  unique(mapping[c("ramp_id", "mz_name", "Adduct", "common_name")])
}

.pn_prepare_metabolite_de <- function(de, db_3, chemical_properties = NULL) {
  names(de)[names(de) == "gene"] <- "mz_name"
  mapping <- .pn_prepare_metabolite_map(db_3, chemical_properties)
  merge(mapping, de, by = "mz_name", sort = FALSE)
}

.pn_best_de_by_cluster <- function(de, id_column) {
  if (is.null(de) || !nrow(de)) return(list())
  id <- as.character(de[[id_column]])
  order_index <- order(de$.cluster_key, id, de$p_val_adj, na.last = TRUE)
  de <- de[order_index, , drop = FALSE]
  key <- paste(de$.cluster_key, de[[id_column]], sep = "\r")
  de <- de[!duplicated(key), , drop = FALSE]
  split(de, de$.cluster_key, drop = TRUE)
}

.pn_leading_edge <- function(x) {
  value <- if (is.list(x)) unlist(x, use.names = FALSE) else as.character(x)
  value <- trimws(value)
  unique(value[!is.na(value) & nzchar(value)])
}

.pn_name_map <- function(source_table, node_ids, chemical_properties = NULL) {
  value <- source_table[
    source_table$rampId %in% node_ids & !is.na(source_table$commonName) &
      nzchar(trimws(source_table$commonName)),
    c("rampId", "commonName"), drop = FALSE
  ]
  value <- value[!duplicated(value$rampId), , drop = FALSE]
  result <- stats::setNames(
    as.character(value$commonName), as.character(value$rampId)
  )
  if (!is.null(chemical_properties)) {
    .pn_assert_columns(
      chemical_properties, c("ramp_id", "common_name"), "chem_props"
    )
    missing_ids <- setdiff(node_ids, names(result))
    chemical <- chemical_properties[
      chemical_properties$ramp_id %in% missing_ids &
        !is.na(chemical_properties$common_name) &
        nzchar(trimws(chemical_properties$common_name)),
      c("ramp_id", "common_name"), drop = FALSE
    ]
    chemical <- chemical[!duplicated(chemical$ramp_id), , drop = FALSE]
    result <- c(
      result,
      stats::setNames(
        as.character(chemical$common_name), as.character(chemical$ramp_id)
      )
    )
  }
  missing_ids <- setdiff(node_ids, names(result))
  if (length(missing_ids) && "sourceId" %in% names(source_table)) {
    fallback <- source_table[
      source_table$rampId %in% missing_ids &
        !is.na(source_table$sourceId) & nzchar(trimws(source_table$sourceId)),
      intersect(c("rampId", "sourceId", "IDtype"), names(source_table)),
      drop = FALSE
    ]
    if (nrow(fallback)) {
      id_type <- if ("IDtype" %in% names(fallback)) {
        tolower(as.character(fallback$IDtype))
      } else {
        rep("", nrow(fallback))
      }
      preferred <- c(
        "gene_symbol", "entrez", "uniprot", "ensembl", "hmdb", "chebi",
        "kegg", "lipidmaps", "refmet"
      )
      priority <- match(id_type, preferred, nomatch = length(preferred) + 1L)
      fallback <- fallback[order(priority), , drop = FALSE]
      fallback <- fallback[!duplicated(fallback$rampId), , drop = FALSE]
      label <- as.character(fallback$sourceId)
      symbol <- tolower(as.character(fallback$IDtype %||% "")) == "gene_symbol"
      label[symbol] <- sub("^[^:]+:", "", label[symbol])
      result <- c(
        result,
        stats::setNames(label, as.character(fallback$rampId))
      )
    }
  }
  result
}

.pn_build_one_network <- function(edges, detected, cluster_key,
                                  gene_lookup, metabolite_lookup,
                                  name_map, max_nodes,
                                  leading_edge = detected,
                                  annotated = character(),
                                  annotation_scores = numeric()) {
  detected <- unique(as.character(detected))
  detected <- detected[grepl("^RAMP_[CG]_", detected)]
  annotated <- unique(as.character(annotated))
  annotated <- annotated[grepl("^RAMP_C_", annotated)]
  retained <- unique(c(detected, annotated))
  edges <- .pn_prune_edges(edges, retained, max_nodes)
  if (!nrow(edges) && !length(retained)) {
    return(list(nodes = .pn_empty_nodes(), links = .pn_empty_links()))
  }
  # A current RaMP leading-edge entity can be absent from an older stored
  # topology. Keep it as an isolated detected node instead of silently
  # discarding the current annotation result.
  node_ids <- unique(c(edges$src, edges$dest, retained))
  group <- ifelse(startsWith(node_ids, "RAMP_G_"), "rna", "mets")
  expr <- rep(NA_real_, length(node_ids))
  node_name <- unname(name_map[node_ids])
  node_name[is.na(node_name) | !nzchar(node_name)] <- node_ids[is.na(node_name) | !nzchar(node_name)]
  display <- node_name

  gene_rows <- gene_lookup[[cluster_key]]
  gene_index <- which(group == "rna")
  if (length(gene_index) && !is.null(gene_rows) && nrow(gene_rows)) {
    match_index <- match(node_ids[gene_index], gene_rows$rampId)
    present <- !is.na(match_index)
    expr[gene_index[present]] <- gene_rows$avg_log2FC[match_index[present]]
    label <- gene_rows$gene[match_index[present]]
    node_name[gene_index[present]] <- label
    display[gene_index[present]] <- label
  }

  met_rows <- metabolite_lookup[[cluster_key]]
  met_index <- which(group == "mets")
  if (length(met_index) && !is.null(met_rows) && nrow(met_rows)) {
    match_index <- match(node_ids[met_index], met_rows$ramp_id)
    present <- !is.na(match_index)
    target <- met_index[present]
    source <- match_index[present]
    expr[target] <- met_rows$avg_log2FC[source]
    valid_name <- !is.na(met_rows$common_name[source]) & nzchar(met_rows$common_name[source])
    node_name[target[valid_name]] <- met_rows$common_name[source[valid_name]]
    display[target] <- paste0(
      met_rows$mz_name[source], " [", met_rows$Adduct[source], "] - ",
      node_name[target]
    )
  }

  nodes <- data.frame(
    id = node_ids,
    group = group,
    expr = expr,
    name = node_name,
    display = display,
    detected = node_ids %in% detected,
    leading_edge = node_ids %in% leading_edge,
    annotated = node_ids %in% annotated,
    annotation_score = unname(annotation_scores[node_ids]),
    stringsAsFactors = FALSE
  )
  links <- if (nrow(edges)) {
    data.frame(
      source = edges$src,
      target = edges$dest,
      colour = edges$colour,
      weight = edges$weight,
      style = edges$linetype,
      head = edges$arrowhead,
      reaction = edges$reaction_name,
      stringsAsFactors = FALSE
    )
  } else {
    .pn_empty_links()
  }
  list(nodes = nodes, links = links)
}

.pn_annotation_score_map <- function(annotation_map) {
  if (is.null(annotation_map) || !nrow(annotation_map)) return(numeric())
  .pn_assert_columns(annotation_map, "ramp_id", "annotation_map")
  score <- if ("Score" %in% names(annotation_map)) {
    suppressWarnings(as.numeric(annotation_map$Score))
  } else {
    rep(NA_real_, nrow(annotation_map))
  }
  ids <- as.character(annotation_map$ramp_id)
  keep <- startsWith(ids, "RAMP_C_") & is.finite(score)
  if (!any(keep)) return(numeric())
  score <- tapply(score[keep], ids[keep], max, na.rm = TRUE)
  score[order(-score, names(score))]
}

.pn_annotated_metabolites_by_pathway <- function(
    annotation_map, pathway_names, pathway_rows,
    analyte_pathway, pathway_table) {
  result <- stats::setNames(vector("list", length(pathway_names)), pathway_names)
  if (is.null(annotation_map) || !nrow(annotation_map)) return(result)
  .pn_assert_columns(annotation_map, "ramp_id", "annotation_map")
  .pn_assert_columns(
    analyte_pathway, c("rampId", "pathwayRampId"), "analytehaspathway"
  )
  .pn_assert_columns(
    pathway_table, c("pathwayRampId", "pathwayName", "sourceId"), "pathway"
  )
  score <- if ("Score" %in% names(annotation_map)) {
    suppressWarnings(as.numeric(annotation_map$Score))
  } else {
    rep(0, nrow(annotation_map))
  }
  score[!is.finite(score)] <- -Inf
  annotation_map$.pathway_score <- score
  annotation_map <- annotation_map[
    order(-annotation_map$.pathway_score, annotation_map$ramp_id),
    ,
    drop = FALSE
  ]

  for (pathway_name in pathway_names) {
    selected_rows <- pathway_rows[
      tolower(as.character(pathway_rows$pathwayName)) == tolower(pathway_name),
      ,
      drop = FALSE
    ]
    name_keys <- unique(tolower(c(
      pathway_name,
      as.character(selected_rows$pathwayName)
    )))
    source_keys <- unique(tolower(as.character(selected_rows$sourceId)))
    pathway_ids <- unique(pathway_table$pathwayRampId[
      tolower(as.character(pathway_table$pathwayName)) %in% name_keys |
        tolower(as.character(pathway_table$sourceId)) %in% source_keys
    ])
    members <- unique(as.character(analyte_pathway$rampId[
      analyte_pathway$pathwayRampId %in% pathway_ids &
        startsWith(as.character(analyte_pathway$rampId), "RAMP_C_")
    ]))
    mapped <- annotation_map[
      annotation_map$ramp_id %in% members,
      ,
      drop = FALSE
    ]
    result[[pathway_name]] <- unique(as.character(mapped$ramp_id))
  }
  result
}

.pn_build_network_collection <- function(pathway_rows, pathway_names,
                                         prepared_edges, clusters,
                                         gene_lookup, metabolite_lookup,
                                         name_map, max_nodes,
                                         metabolite_detection = "leading_edge",
                                         annotated_metabolites = NULL,
                                         annotation_scores = numeric(),
                                         annotation_score_threshold = 0.05) {
  networks <- vector("list", length(clusters))
  names(networks) <- clusters
  matrix_ids <- character()
  for (i in seq_along(clusters)) {
    cluster_key <- tolower(as.character(clusters[[i]]))
    cluster_rows <- pathway_rows[
      tolower(as.character(pathway_rows$Cluster_id)) == cluster_key,
      , drop = FALSE
    ]
    cluster_networks <- vector("list", length(pathway_names))
    names(cluster_networks) <- pathway_names
    for (j in seq_along(pathway_names)) {
      pathway_key <- tolower(pathway_names[[j]])
      is_present <- any(tolower(cluster_rows$pathwayName) == pathway_key)
      if (!is_present) {
        cluster_networks[[j]] <- list(nodes = .pn_empty_nodes(), links = .pn_empty_links())
        next
      }
      pathway_cluster_rows <- cluster_rows[
        tolower(cluster_rows$pathwayName) == pathway_key,
        , drop = FALSE
      ]
      leading_edge <- unique(unlist(
        lapply(pathway_cluster_rows$leadingEdge, .pn_leading_edge),
        use.names = FALSE
      ))
      detected <- leading_edge
      if (identical(metabolite_detection, "annotated")) {
        pathway_annotated <- annotated_metabolites[[pathway_names[[j]]]]
        pathway_scores <- unname(annotation_scores[pathway_annotated])
        score_keep <- is.finite(pathway_scores) &
          pathway_scores >= annotation_score_threshold
        detected <- unique(c(
          detected,
          pathway_annotated[score_keep]
        ))
      }
      network <- .pn_build_one_network(
        prepared_edges[[pathway_names[[j]]]], detected, cluster_key,
        gene_lookup, metabolite_lookup, name_map, max_nodes,
        leading_edge = leading_edge,
        annotated = annotated_metabolites[[pathway_names[[j]]]],
        annotation_scores = annotation_scores
      )
      cluster_networks[[j]] <- network
      matrix_ids <- c(
        matrix_ids,
        network$nodes$id[
          network$nodes$leading_edge | network$nodes$annotated
        ]
      )
    }
    networks[[i]] <- cluster_networks
  }
  list(networks = networks, matrix_ids = unique(matrix_ids))
}

.pn_layer_data <- function(object, assay, layer) {
  assay_object <- object[[assay]]
  value <- tryCatch(
    SeuratObject::LayerData(assay_object, layer = layer),
    error = function(e) NULL
  )
  if (is.null(value) && methods::is(assay_object, "Assay5") &&
      layer %in% names(assay_object@layers)) {
    value <- assay_object@layers[[layer]]
  }
  if (is.null(value)) {
    value <- tryCatch(
      SeuratObject::GetAssayData(object, assay = assay, slot = layer),
      error = function(e) NULL
    )
  }
  if (is.null(value)) stop("Unable to read layer '", layer, "' from assay '", assay, "'.")
  value
}

.pn_coordinate_columns <- function(coordinates) {
  lower <- tolower(names(coordinates))
  pick <- function(candidates) {
    index <- match(candidates, lower, nomatch = 0L)
    index <- index[index > 0L]
    if (!length(index)) return(NA_integer_)
    index[[1]]
  }
  x <- pick(c("x", "imagecol", "col", "pxl_col_in_fullres"))
  y <- pick(c("y", "imagerow", "row", "pxl_row_in_fullres"))
  if (is.na(x) || is.na(y)) {
    stop("Tissue coordinates must contain x/y or imagecol/imagerow columns.")
  }
  c(x = x, y = y)
}

.pn_matrix_expression <- function(matrix, feature_names, cell_names) {
  if (!length(feature_names)) return(list())
  feature_index <- match(tolower(feature_names), tolower(rownames(matrix)))
  cell_index <- if (!is.null(colnames(matrix)) && length(cell_names)) {
    match(cell_names, colnames(matrix))
  } else {
    seq_len(min(ncol(matrix), length(cell_names)))
  }
  result <- vector("list", length(feature_names))
  for (i in seq_along(feature_names)) {
    if (is.na(feature_index[[i]])) next
    values <- rep(NA_real_, length(cell_names))
    present <- !is.na(cell_index)
    values[present] <- as.numeric(matrix[feature_index[[i]], cell_index[present], drop = TRUE])
    result[[i]] <- values
  }
  result
}

.pn_spatial_payload <- function(object, ident, image, matrix_ids,
                                gene_de, metabolite_de,
                                gene_matrix = NULL, metabolite_matrix = NULL,
                                max_spatial_points = 50000L) {
  coordinates <- as.data.frame(Seurat::GetTissueCoordinates(object, image = image))
  columns <- .pn_coordinate_columns(coordinates)
  x <- suppressWarnings(as.numeric(coordinates[[columns[["x"]]]]))
  y <- suppressWarnings(as.numeric(coordinates[[columns[["y"]]]]))
  cell_column <- match("cell", tolower(names(coordinates)), nomatch = 0L)
  coordinate_cells <- if (cell_column > 0L) {
    as.character(coordinates[[cell_column]])
  } else {
    rownames(coordinates)
  }
  metadata <- object@meta.data
  assignment <- if (length(coordinate_cells) && all(coordinate_cells %in% rownames(metadata))) {
    as.character(metadata[coordinate_cells, ident])
  } else {
    as.character(metadata[[ident]])[seq_len(nrow(coordinates))]
  }
  keep <- is.finite(x) & is.finite(y) & !is.na(assignment)
  keep_index <- which(keep)
  if (is.finite(max_spatial_points)) max_spatial_points <- as.integer(max_spatial_points)
  if (is.finite(max_spatial_points) && max_spatial_points > 0L &&
      length(keep_index) > max_spatial_points) {
    keep_index <- keep_index[unique(round(seq(
      1, length(keep_index), length.out = max_spatial_points
    )))]
  }
  x <- x[keep_index]
  y <- y[keep_index]
  assignment <- assignment[keep_index]
  if (!length(coordinate_cells)) {
    available_cells <- colnames(gene_matrix %||% metabolite_matrix)
    coordinate_cells <- available_cells[seq_len(nrow(coordinates))]
  }
  coordinate_cells <- coordinate_cells[keep_index]
  scale01 <- function(value) {
    range <- range(value, finite = TRUE)
    if (!is.finite(diff(range)) || diff(range) == 0) return(rep(0.5, length(value)))
    (value - range[[1]]) / diff(range)
  }

  gene_ids <- matrix_ids[startsWith(matrix_ids, "RAMP_G_")]
  met_ids <- matrix_ids[startsWith(matrix_ids, "RAMP_C_")]
  gene_expression <- list()
  metabolite_expression <- list()
  if (length(gene_ids) && !is.null(gene_matrix) && !is.null(gene_de) && nrow(gene_de)) {
    best <- gene_de[order(gene_de$p_val_adj), , drop = FALSE]
    best <- best[!duplicated(best$rampId), , drop = FALSE]
    features <- best$gene[match(gene_ids, best$rampId)]
    values <- .pn_matrix_expression(gene_matrix, features, coordinate_cells)
    names(values) <- gene_ids
    gene_expression <- Filter(Negate(is.null), values)
  }
  if (length(met_ids) && !is.null(metabolite_matrix) &&
      !is.null(metabolite_de) && nrow(metabolite_de)) {
    best <- metabolite_de[order(metabolite_de$p_val_adj), , drop = FALSE]
    best <- best[!duplicated(best$ramp_id), , drop = FALSE]
    features <- best$mz_name[match(met_ids, best$ramp_id)]
    values <- .pn_matrix_expression(metabolite_matrix, features, coordinate_cells)
    names(values) <- met_ids
    metabolite_expression <- Filter(Negate(is.null), values)
  }
  list(
    coordinates = data.frame(x = scale01(x), y = scale01(y)),
    clusters = assignment,
    expression = list(rna = gene_expression, mets = metabolite_expression),
    displayed_points = length(keep_index),
    total_points = sum(keep)
  )
}

.pn_fc_limit <- function(DE.list) {
  value <- unlist(lapply(DE.list, function(x) x$avg_log2FC), use.names = FALSE)
  value <- abs(value[is.finite(value)])
  if (!length(value)) return(1)
  max(1, as.numeric(stats::quantile(value, 0.98, names = FALSE, type = 7)))
}

.pn_template_path <- function() {
  installed <- system.file("templates", "pathway_network.html", package = "SpaMTP")
  if (nzchar(installed) && file.exists(installed)) return(installed)
  development <- file.path("inst", "templates", "pathway_network.html")
  if (file.exists(development)) return(development)
  stop("Cannot locate the pathway network HTML template.")
}

.pn_render_html <- function(payload, template_path = .pn_template_path()) {
  template <- paste(readLines(template_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  payload$pathways <- I(as.character(payload$pathways))
  payload$clusters <- I(as.character(payload$clusters))
  payload$networks <- unname(lapply(payload$networks, function(cluster) {
    unname(cluster)
  }))
  payload$spatial$clusters <- I(as.character(payload$spatial$clusters))
  payload$spatial$expression <- lapply(payload$spatial$expression, function(group) {
    lapply(group, function(values) I(as.numeric(values)))
  })
  payload_json <- jsonlite::toJSON(
    payload, auto_unbox = TRUE, dataframe = "rows", null = "null",
    na = "null", digits = 6, POSIXt = "ISO8601"
  )
  payload_json <- gsub("</", "<\\/", payload_json, fixed = TRUE)
  placeholder <- "__SPAMTP_PAYLOAD__"
  if (!grepl(placeholder, template, fixed = TRUE)) {
    stop("Pathway network template is missing its payload placeholder.")
  }
  sub(placeholder, payload_json, template, fixed = TRUE)
}
