#' Construct an interactive pathway network visualization
#'
#' Builds a self-contained data payload and an interactive D3 visualization for
#' selected pathways across spatial clusters. Topology lookup, edge preparation,
#' differential-expression matching, and matrix extraction are indexed so they
#' are not repeated for every cluster-pathway pair.
#'
#' @param SpaMTP A `SpaMTP` Seurat object containing spatial metabolomics and/or
#'   spatial transcriptomics data. Metabolomics data must first be annotated by
#'   [AnnotateSM()].
#' @param ident Metadata column used to identify spatial clusters or regions.
#' @param regpathway Data frame returned by [FindRegionalPathways()].
#' @param DE.list One differential-expression data frame per requested analyte
#'   type. Data frames must contain `cluster`, `gene`, `avg_log2FC` (or
#'   `logFC`), and `p_val_adj` (or `FDR`). A named list is recommended.
#' @param selected_pathways Optional pathway names or source IDs. Matching is
#'   case-insensitive. When `NULL`, the most important pathways are selected by
#'   summed absolute NES.
#' @param path Output directory for the generated HTML file.
#' @param SM_assay Spatial-metabolomics assay name.
#' @param ST_assay Spatial-transcriptomics assay name.
#' @param SM_slot Layer containing spatial-metabolomics values.
#' @param ST_slot Layer containing spatial-transcriptomics values.
#' @param colour_palette Colours used for the spatial abundance raster.
#' @param analyte_types One or both of `"genes"` and `"metabolites"`.
#' @param annotation_source Metabolite annotation provenance. The default,
#'   `"current"`, requires the indexed, scored RaMP annotation output.
#'   `"auto"` permits a warned legacy fallback and `"legacy"` requests it
#'   explicitly.
#' @param annotation_score_threshold Minimum indexed annotation score used for
#'   metabolite matching. When `NULL`, the threshold recorded by
#'   `FindRegionalPathways()` is reused (default = `NULL`).
#' @param annotation_score_floor Lowest annotation score embedded in the HTML
#'   for interactive filtering. It is automatically lowered when the initial
#'   `annotation_score_threshold` is smaller (default = `0.01`).
#' @param metabolite_detection Which metabolite nodes receive detected styling.
#'   `"leading_edge"` uses only GSEA leading-edge metabolites;
#'   `"annotated"` also includes every score-filtered annotation belonging to
#'   the selected pathway, regardless of DE significance.
#' @param image Spatial image/FOV passed to `Seurat::GetTissueCoordinates()`.
#' @param verbose Display progress messages.
#' @param top_n_pathways Number of pathways selected when `selected_pathways` is
#'   `NULL`.
#' @param max_nodes Maximum number of nodes per cluster-pathway view. Detected
#'   leading-edge analytes are always retained, followed by their neighbours and
#'   high-degree nodes. Use `Inf` to disable pruning.
#' @param label_mode Initial node-label mode: detected nodes only, all nodes, or
#'   no labels. It can also be changed interactively.
#' @param layout_mode Initial network layout. `"repulsion"` maximises spacing,
#'   `"force"` uses a balanced force-directed layout, `"radial"` separates
#'   analyte classes into rings, and `"bipartite"` separates genes and
#'   metabolites vertically. It can also be changed interactively.
#' @param max_spatial_points Maximum number of spatial points embedded in the
#'   HTML. Larger datasets are deterministically downsampled for responsive
#'   browser rendering. Use `Inf` to retain all points.
#' @param database Optional named list of database resources, normally created
#'   by [LoadSpaMTPDatabase()].
#' @param database_version SpaMTPdb/RaMP version used for pathway lookup.
#' @param database_source Database source; see [LoadSpaMTPDatabase()].
#' @param database_local_dir Optional staged SpaMTPdb resource directory.
#'
#' @return Invisibly returns the generated HTML file path.
#' @export
#'
#' @examples
#' utils::str(formals(PathwayNetworkPlots))
#' # PathwayNetworkPlots(
#' #   SpaMTP, ident = "Custom_ident", regpathway = regpathway,
#' #   DE.list = DE.list, selected_pathways = "WP1902"
#' # )
PathwayNetworkPlots <- function(SpaMTP,
                                ident,
                                regpathway,
                                DE.list,
                                selected_pathways = NULL,
                                path = getwd(),
                                SM_slot = "counts",
                                ST_slot = "counts",
                                colour_palette = NULL,
                                SM_assay = "SPM",
                                ST_assay = "SPT",
                                analyte_types = c("genes", "metabolites"),
                                annotation_source = c("current", "auto", "legacy"),
                                annotation_score_threshold = NULL,
                                annotation_score_floor = 0.01,
                                metabolite_detection = c("leading_edge", "annotated"),
                                image = "slice1",
                                verbose = TRUE,
                                top_n_pathways = 10L,
                                max_nodes = 500L,
                                label_mode = c("detected", "all", "none"),
                                max_spatial_points = 50000L,
                                layout_mode = c("repulsion", "force", "radial", "bipartite"),
                                database = NULL,
                                database_version = "latest",
                                database_source = c("auto", "spamtpdb"),
                                database_local_dir = NULL) {
  analyte_types <- match.arg(
    analyte_types, c("genes", "metabolites"), several.ok = TRUE
  )
  analyte_types <- unique(analyte_types)
  annotation_source <- match.arg(annotation_source)
  metabolite_detection <- match.arg(metabolite_detection)
  label_mode <- match.arg(label_mode)
  layout_mode <- match.arg(layout_mode)
  database_source <- match.arg(database_source)
  database_resources <- .spamtp_db_bundle(
    c(
      "chem_props", "source_df", "analytehaspathway", "pathway",
      "ramp_db_metadata"
    ),
    database = database,
    version = database_version,
    source = database_source,
    local_dir = database_local_dir
  )
  chem_props <- database_resources$chem_props
  source_df <- database_resources$source_df
  analytehaspathway <- database_resources$analytehaspathway
  pathway <- database_resources$pathway
  if (!is.character(ident) || length(ident) != 1L || !nzchar(ident)) {
    stop("ident must be one metadata column name.")
  }
  if (!ident %in% names(SpaMTP@meta.data)) {
    stop("Metadata column '", ident, "' was not found in SpaMTP@meta.data.")
  }
  if (!dir.exists(path)) stop("Output directory does not exist: ", path)
  if (file.access(path, mode = 2) != 0) {
    stop("Cannot write to output directory: ", path)
  }
  if (!is.numeric(max_nodes) || length(max_nodes) != 1L ||
      is.na(max_nodes) || max_nodes <= 0) {
    stop("max_nodes must be positive or Inf.")
  }
  if (!is.numeric(max_spatial_points) || length(max_spatial_points) != 1L ||
      is.na(max_spatial_points) || max_spatial_points <= 0) {
    stop("max_spatial_points must be positive or Inf.")
  }

  assays <- names(SpaMTP@assays)
  if ("genes" %in% analyte_types && !ST_assay %in% assays) {
    stop("Gene visualization requires assay '", ST_assay, "'.")
  }
  if ("metabolites" %in% analyte_types && !SM_assay %in% assays) {
    stop("Metabolite visualization requires assay '", SM_assay, "'.")
  }
  if (all(c("genes", "metabolites") %in% analyte_types) &&
      ncol(SpaMTP[[SM_assay]]) != ncol(SpaMTP[[ST_assay]])) {
    stop("SM and ST assays must contain the same number of aligned spatial observations.")
  }
  regpathway_annotation <- attr(regpathway, "annotation_metadata", exact = TRUE)
  if (is.null(annotation_score_threshold)) {
    annotation_score_threshold <- regpathway_annotation$analysis_score_threshold %||% 0.05
  }
  if (!is.numeric(annotation_score_threshold) ||
      length(annotation_score_threshold) != 1L ||
      !is.finite(annotation_score_threshold) || annotation_score_threshold < 0) {
    stop("annotation_score_threshold must be one non-negative finite number.")
  }
  if (!is.numeric(annotation_score_floor) ||
      length(annotation_score_floor) != 1L ||
      !is.finite(annotation_score_floor) || annotation_score_floor < 0) {
    stop("annotation_score_floor must be one non-negative finite number.")
  }
  annotation_score_floor <- min(
    annotation_score_floor, annotation_score_threshold
  )
  if ("metabolites" %in% analyte_types && annotation_source == "current") {
    if (!.pn_is_current_annotation_metadata(regpathway_annotation)) {
      stop(
        "regpathway was not generated with a current scored annotation ",
        "pipeline. Re-run FindRegionalPathways() before plotting, or ",
        "explicitly select annotation_source = 'legacy'."
      )
    }
    current_ramp <- as.character(
      database_resources$ramp_db_metadata$ramp_version %||%
        .annotation_ramp_version()
    )
    if (!identical(regpathway_annotation$ramp_version, current_ramp)) {
      stop(
        "regpathway uses RaMP ", regpathway_annotation$ramp_version %||% "unknown",
        " but this SpaMTP build contains RaMP ", current_ramp,
        ". Re-run FindRegionalPathways()."
      )
    }
  }

  verbose_message("Selecting enriched pathways ... ", verbose)
  pathway_rows <- .pn_select_pathway_rows(
    as.data.frame(regpathway), selected_pathways, top_n_pathways
  )
  requested_sources <- unique(vapply(pathway_rows$type, .pn_type_key, character(1)))
  if (anyNA(requested_sources)) {
    requested_sources <- c("wiki", "reactome", "kegg", "hmdb")
  }
  topology_resource_names <- c(
    wiki = "ramp_wikipathway",
    reactome = "ramp_reactome",
    kegg = "ramp_kegg",
    hmdb = "ramp_hmdb"
  )[requested_sources]
  topology_resources <- .spamtp_db_bundle(
    unname(topology_resource_names),
    database = database,
    version = database_version,
    source = database_source,
    local_dir = database_local_dir
  )
  topology_sources <- list()
  for (source_name in names(topology_resource_names)) {
    topology_sources[[source_name]] <- topology_resources[[topology_resource_names[[source_name]]]]
  }
  catalog <- .pn_topology_catalog(topology_sources)
  resolved <- .pn_resolve_topologies(pathway_rows, catalog, topology_sources)
  if (any(!resolved$found)) {
    missing <- resolved$pathways$pathwayName[!resolved$found]
    verbose_message(
      paste0(
        "No stored topology was found for: ",
        paste(missing, collapse = "; ")
      ),
      verbose
    )
  }
  topologies <- resolved$topologies[resolved$found]
  if (!length(topologies)) {
    stop("None of the selected pathways has a stored topological structure.")
  }
  pathway_names <- names(topologies)
  pathway_rows <- pathway_rows[
    tolower(pathway_rows$pathwayName) %in% tolower(pathway_names),
    , drop = FALSE
  ]

  verbose_message("Preparing pathway edges once per topology ... ", verbose)
  prepared_edges <- lapply(topologies, .pn_prepare_edges, reaction_types = reaction_type)
  usable <- lengths(lapply(prepared_edges, function(x) unique(c(x$src, x$dest)))) > 0L
  if (any(!usable)) {
    verbose_message(
      paste0(
        "No valid edges were found for: ",
        paste(names(prepared_edges)[!usable], collapse = "; ")
      ),
      verbose
    )
  }
  prepared_edges <- prepared_edges[usable]
  pathway_names <- names(prepared_edges)
  if (!length(pathway_names)) stop("Selected topologies contain no valid network edges.")
  pathway_rows <- pathway_rows[
    tolower(pathway_rows$pathwayName) %in% tolower(pathway_names),
    , drop = FALSE
  ]

  verbose_message("Indexing differential results ... ", verbose)
  differential <- .pn_normalise_de_list(DE.list, analyte_types)
  gene_de <- NULL
  metabolite_de <- NULL
  metabolite_annotations <- NULL
  annotation_metadata <- NULL
  if ("genes" %in% analyte_types) {
    gene_de <- .pn_prepare_gene_de(differential[["genes"]], source_df)
  }
  if ("metabolites" %in% analyte_types) {
    db_3 <- .resolve_pathway_metabolite_annotations(
      SpaMTP,
      annotation_source = annotation_source,
      score_threshold = annotation_score_floor,
      chemical_properties = chem_props
    )
    metabolite_annotations <- db_3
    annotation_metadata <- attr(db_3, "annotation_metadata")
    metabolite_de <- .pn_prepare_metabolite_de(
      differential[["metabolites"]], db_3, chem_props
    )
  }
  gene_lookup <- .pn_best_de_by_cluster(gene_de, "rampId")
  metabolite_lookup <- .pn_best_de_by_cluster(metabolite_de, "ramp_id")

  all_node_ids <- unique(unlist(lapply(prepared_edges, function(x) {
    c(x$src, x$dest)
  }), use.names = FALSE))
  annotated_metabolites <- stats::setNames(
    vector("list", length(pathway_names)), pathway_names
  )
  annotation_scores <- numeric()
  if ("metabolites" %in% analyte_types) {
    annotated_metabolites <- .pn_annotated_metabolites_by_pathway(
      metabolite_annotations,
      pathway_names,
      pathway_rows,
      analytehaspathway,
      pathway
    )
    annotation_scores <- .pn_annotation_score_map(metabolite_annotations)
  }
  all_node_ids <- unique(c(
    all_node_ids,
    unlist(annotated_metabolites, use.names = FALSE)
  ))
  node_names <- .pn_name_map(source_df, all_node_ids, chem_props)
  clusters <- naturalsort::naturalsort(unique(as.character(pathway_rows$Cluster_id)))
  built <- .pn_build_network_collection(
    pathway_rows, pathway_names, prepared_edges, clusters,
    gene_lookup, metabolite_lookup, node_names, max_nodes,
    metabolite_detection = metabolite_detection,
    annotated_metabolites = annotated_metabolites,
    annotation_scores = annotation_scores,
    annotation_score_threshold = annotation_score_threshold
  )

  verbose_message("Extracting only plotted assay features ... ", verbose)
  gene_matrix <- if ("genes" %in% analyte_types) {
    .pn_layer_data(SpaMTP, ST_assay, ST_slot)
  } else NULL
  metabolite_matrix <- if ("metabolites" %in% analyte_types) {
    .pn_layer_data(SpaMTP, SM_assay, SM_slot)
  } else NULL
  spatial <- .pn_spatial_payload(
    SpaMTP, ident, image, built$matrix_ids,
    gene_de, metabolite_de, gene_matrix, metabolite_matrix,
    max_spatial_points
  )

  if (is.null(colour_palette)) {
    colour_palette <- grDevices::colorRampPalette(
      c("#0d0887", "#7e03a8", "#cc4778", "#f89540", "#f0f921")
    )(100)
  }
  if (!is.character(colour_palette) || length(colour_palette) < 2L) {
    stop("colour_palette must contain at least two colours.")
  }

  payload <- list(
    title = paste("SpaMTP pathway networks -", ident),
    pathways = pathway_names,
    clusters = clusters,
    networks = built$networks,
    spatial = spatial,
    palette = colour_palette,
    fc_limit = .pn_fc_limit(differential),
    label_mode = label_mode,
    layout_mode = layout_mode,
    metadata = list(
      max_nodes = max_nodes,
      annotation = annotation_metadata,
      annotation_score_threshold = annotation_score_threshold,
      annotation_score_floor = annotation_score_floor,
      annotation_score_ceiling = max(
        c(annotation_score_threshold, annotation_scores), na.rm = TRUE
      ),
      metabolite_detection = metabolite_detection,
      annotation_controls = "metabolites" %in% analyte_types,
      generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
  )
  html <- .pn_render_html(payload)

  file_stem <- gsub("[^A-Za-z0-9._-]+", "_", ident)
  return_name <- paste0(file_stem, "_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".html")
  full_path <- file.path(path, return_name)
  writeLines(html, full_path, useBytes = TRUE)
  message("Pathway network written to ", full_path)
  invisible(full_path)
}
