.spamtp_db_cache <- new.env(parent = emptyenv())

.spamtp_db_legacy_names <- c(
  chem_props = "chem_props",
  source_df = "source_df",
  analyte = "analyte",
  analytehaspathway = "analytehaspathway",
  pathway = "pathway",
  ramp_db_metadata = "ramp_db_metadata",
  ramp_hmdb = "RAMP_hmdb",
  ramp_kegg = "RAMP_kegg",
  ramp_reactome = "RAMP_Reactome",
  ramp_wikipathway = "RAMP_wikipathway",
  hmdb_db = "HMDB_db",
  chebi_db = "Chebi_db",
  lipidmaps_db = "Lipidmaps_db",
  gnps_db = "GNPS_db",
  filtered_fmp10 = "filtered_fmp10"
)

# Structure features are intentionally separate from the pruned chemical table.
.spamtp_db_legacy_names <- c(
  .spamtp_db_legacy_names,
  smiles_features = "smiles_features"
)

.spamtp_db_normalise_resources <- function(resources) {
  resources <- tolower(trimws(as.character(resources)))
  resources <- unique(resources[nzchar(resources)])
  unknown <- setdiff(resources, names(.spamtp_db_legacy_names))
  if (length(unknown)) {
    stop(
      "Unknown SpaMTP database resource(s): ",
      paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }
  resources
}

.spamtp_db_cache_key <- function(resource, version, source, local_dir) {
  local_key <- if (is.null(local_dir)) "<default>" else {
    normalizePath(local_dir, mustWork = FALSE)
  }
  paste(
    resource,
    version %||% "latest",
    source,
    local_key,
    sep = "\r"
  )
}

.spamtp_db_label <- function(value, fallback = "user-supplied database") {
  metadata <- attr(value, "spamtp_database", exact = TRUE)
  if (!is.list(metadata)) return(fallback)
  version <- metadata$version %||% "unknown"
  paste0(metadata$source %||% "SpaMTP", " RaMP ", version, " ", metadata$resource)
}

.spamtp_db_resource <- function(resource,
                                version = "latest",
                                source = c("auto", "bundled", "local"),
                                local_dir = NULL,
                                hub = NULL,
                                offline = FALSE,
                                refresh = FALSE) {
  resource <- .spamtp_db_normalise_resources(resource)
  if (length(resource) != 1L) {
    stop("resource must identify exactly one database resource.", call. = FALSE)
  }
  source <- match.arg(source)
  if (!is.null(hub)) {
    stop("This SpaMTP release uses bundled or local databases, not a Hub.",
         call. = FALSE)
  }
  if (source == "auto") source <- if (is.null(local_dir)) "bundled" else "local"
  if (source == "local" && is.null(local_dir)) {
    stop("source = 'local' requires local_dir.", call. = FALSE)
  }
  if (source == "bundled" && !is.null(local_dir)) {
    stop("Use source = 'local' or 'auto' with local_dir.", call. = FALSE)
  }
  if (is.null(version)) version <- "latest"
  if (!is.character(version) || length(version) != 1L ||
      is.na(version) || !nzchar(version)) {
    stop("version must be a single non-empty string.", call. = FALSE)
  }
  key <- .spamtp_db_cache_key(resource, version, source, local_dir)
  if (!isTRUE(refresh) && exists(key, envir = .spamtp_db_cache, inherits = FALSE)) {
    return(get(key, envir = .spamtp_db_cache, inherits = FALSE))
  }

  if (source == "bundled") {
    registry <- SpaMTPDatabaseInfo(version = version)
    row <- registry[registry$resource == resource, , drop = FALSE]
    if (nrow(row) != 1L) {
      stop("No bundled resource '", resource, "' for version '", version, "'.",
           call. = FALSE)
    }
    environment <- new.env(parent = emptyenv())
    utils::data(list = row$source_object, package = "SpaMTP", envir = environment)
    if (!exists(row$source_object, envir = environment, inherits = FALSE)) {
      stop("Bundled database is missing: ", resource, ".", call. = FALSE)
    }
    value <- get(row$source_object, envir = environment, inherits = FALSE)
    resolved_version <- row$version
  } else {
    # Local RDS files are an explicit user override, never a network fallback.
    directory <- normalizePath(local_dir, mustWork = TRUE)
    local_version <- if (version == "latest") {
      unique(SpaMTPDatabaseInfo()$version)
    } else version
    paths <- file.path(c(directory, file.path(directory, local_version)),
                       paste0(resource, ".rds"))
    paths <- paths[file.exists(paths)]
    if (!length(paths)) {
      stop("Local database resource is missing: ", resource, ".rds.",
           call. = FALSE)
    }
    value <- readRDS(paths[[1L]])
    metadata_path <- file.path(dirname(paths[[1L]]), "ramp_db_metadata.rds")
    local_metadata <- if (file.exists(metadata_path)) readRDS(metadata_path) else NULL
    resolved_version <- if (is.list(local_metadata)) local_metadata$ramp_version else NULL
    if (!is.null(resolved_version) && version != "latest" &&
        !identical(as.character(resolved_version), version)) {
      stop("Local database metadata does not match requested version '",
           version, "'.", call. = FALSE)
    }
    resolved_version <- resolved_version %||%
      if (version == "latest") "unversioned" else version
  }

  value <- .spamtp_repair_pathway_interactions(value, resource)
  attr(value, "spamtp_database") <- list(
    resource = resource,
    version = as.character(resolved_version),
    source = source,
    local_dir = local_dir,
    offline = TRUE
  )
  assign(key, value, envir = .spamtp_db_cache)
  value
}

.spamtp_db_bundle <- function(resources,
                              database = NULL,
                              version = "latest",
                              source = c("auto", "bundled", "local"),
                              local_dir = NULL,
                              hub = NULL,
                              offline = FALSE,
                              refresh = FALSE) {
  resources <- .spamtp_db_normalise_resources(resources)
  source <- match.arg(source)
  if (!is.null(database)) {
    if (!is.list(database) || is.data.frame(database) || is.null(names(database))) {
      stop("database must be a named list of SpaMTP database resources.", call. = FALSE)
    }
    missing <- setdiff(resources, tolower(names(database)))
    if (length(missing)) {
      stop(
        "database is missing resource(s): ", paste(missing, collapse = ", "),
        ".", call. = FALSE
      )
    }
    names(database) <- tolower(names(database))
    values <- lapply(resources, function(resource) {
      .spamtp_repair_pathway_interactions(database[[resource]], resource)
    })
    return(stats::setNames(values, resources))
  }

  values <- lapply(resources, function(resource) {
    .spamtp_db_resource(
      resource = resource,
      version = version,
      source = source,
      local_dir = local_dir,
      hub = hub,
      offline = offline,
      refresh = refresh
    )
  })
  stats::setNames(values, resources)
}

#' Load versioned SpaMTP annotation resources
#'
#' Loads annotation resources bundled with SpaMTP without another data package,
#' a Hub lookup, or a download. Resources are cached for the current R session.
#' A named custom bundle or local RDS directory can be supplied for
#' user-curated workflows.
#'
#' @details
#' Bundled RaMP 3.0.7 graphs already include corrected interaction labels and
#' directions. Exact affected local or custom topology resources receive a checksum-guarded
#' correction for historical interaction-code and direction recycling. Source
#' labels are retained in `source_reaction_type`; the resource attribute
#' `spamtp_interaction_repair` records the correction separately from the
#' original download metadata. Published files are unchanged. See
#' `vignette("Pathway_Database_Integration", package = "SpaMTP")` for provenance
#' and reproducible build instructions.
#'
#' @param resources Character vector of resource names. Use
#'   [SpaMTPDatabaseInfo()] to list valid names.
#' @param version Database snapshot version, or `"latest"` for the bundled
#'   snapshot. For local files, declared metadata must match an explicit version.
#' @param source Database source. `"auto"` uses bundled data unless `local_dir`
#'   is supplied. `"bundled"` uses installed data; `"local"` requires local RDS
#'   files and never falls back to another source.
#' @param database Optional named list containing the requested resources. When
#'   supplied, no Hub lookup is performed.
#' @param local_dir Optional directory containing `<resource>.rds` files,
#'   either directly or in a version-named subdirectory. Without local version
#'   metadata, `"latest"` files are labelled `"unversioned"`.
#' @param hub Retained for call compatibility; must be `NULL` in this release.
#' @param offline Retained for call compatibility. All resource loading is
#'   offline, regardless of this argument.
#' @param refresh If `TRUE`, bypass SpaMTP's in-session resource cache.
#'
#' @return A named list containing the requested resources.
#' @export
#'
#' @examples
#' utils::str(formals(LoadSpaMTPDatabase))
#' example_database <- list(
#'   ramp_db_metadata = list(ramp_version = "example")
#' )
#' database <- LoadSpaMTPDatabase(
#'   "ramp_db_metadata",
#'   database = example_database
#' )
#' names(database)
LoadSpaMTPDatabase <- function(
    resources = c(
      "chem_props", "source_df", "analyte", "analytehaspathway", "pathway"
    ),
    version = "latest",
    source = c("auto", "bundled", "local"),
    database = NULL,
    local_dir = NULL,
    hub = NULL,
    offline = FALSE,
    refresh = FALSE) {
  .spamtp_db_bundle(
    resources = resources,
    database = database,
    version = version,
    source = match.arg(source),
    local_dir = local_dir,
    hub = hub,
    offline = offline,
    refresh = refresh
  )
}

#' Inspect SpaMTP database resources
#'
#' @param version Optional bundled database snapshot version. `NULL` and
#'   `"latest"` list the snapshot shipped with this SpaMTP release.
#'
#' @return A data frame describing bundled resources, including their version,
#'   dimensions, dataset names, and source-package file sizes and MD5 checksums.
#'   An unavailable version returns a zero-row data frame.
#' @export
#'
#' @examples
#' utils::str(formals(SpaMTPDatabaseInfo))
#' SpaMTPDatabaseInfo()
SpaMTPDatabaseInfo <- function(version = NULL) {
  path <- system.file("extdata", "bundled_database_manifest.csv", package = "SpaMTP")
  if (!nzchar(path)) stop("Bundled database manifest is missing.", call. = FALSE)
  registry <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!is.null(version)) {
    if (!is.character(version) || length(version) != 1L ||
        is.na(version) || !nzchar(version)) {
      stop("version must be NULL or a single non-empty string.", call. = FALSE)
    }
    if (version != "latest") registry <- registry[registry$version == version, , drop = FALSE]
  }
  registry
}
