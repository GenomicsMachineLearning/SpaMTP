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
                                source = c("auto", "spamtpdb"),
                                local_dir = NULL,
                                hub = NULL,
                                offline = FALSE,
                                refresh = FALSE) {
  resource <- .spamtp_db_normalise_resources(resource)
  if (length(resource) != 1L) {
    stop("resource must identify exactly one database resource.", call. = FALSE)
  }
  source <- match.arg(source)
  key <- .spamtp_db_cache_key(resource, version, source, local_dir)
  if (!isTRUE(refresh) && exists(key, envir = .spamtp_db_cache, inherits = FALSE)) {
    return(get(key, envir = .spamtp_db_cache, inherits = FALSE))
  }

  if (!requireNamespace("SpaMTPdb", quietly = TRUE)) {
    stop(
      "Loading versioned annotation resources requires the SpaMTPdb package. ",
      "Install SpaMTPdb or supply a named custom database bundle.",
      call. = FALSE
    )
  }

  value <- SpaMTPdb::SpaMTPdbResource(
    resource = resource,
    version = version,
    local_dir = local_dir,
    hub = hub,
    offline = offline
  )
  resource_metadata <- SpaMTPdb::SpaMTPdbResource(
    resource = resource,
    version = version,
    metadata = TRUE
  )

  attr(value, "spamtp_database") <- list(
    resource = resource,
    version = as.character(resource_metadata$version[[1L]]),
    source = "spamtpdb",
    local_dir = local_dir,
    offline = offline
  )
  assign(key, value, envir = .spamtp_db_cache)
  value
}

.spamtp_db_bundle <- function(resources,
                              database = NULL,
                              version = "latest",
                              source = c("auto", "spamtpdb"),
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
    return(database[resources])
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
#' Loads a coherent group of annotation resources from [SpaMTPdb]. Retrieved
#' resources are cached for the current R session. A named custom bundle can be
#' supplied for offline, testing, or user-curated workflows.
#'
#' @param resources Character vector of resource names. Use
#'   [SpaMTPDatabaseInfo()] to list valid names.
#' @param version SpaMTPdb/RaMP resource version, or `"latest"`.
#' @param source Database source. `"auto"` and `"spamtpdb"` resolve versioned
#'   resources through SpaMTPdb.
#' @param database Optional named list containing the requested resources. When
#'   supplied, no Hub lookup is performed.
#' @param local_dir Optional directory containing staged SpaMTPdb `.rds` files.
#' @param hub Optional pre-created `AnnotationHub` object passed to SpaMTPdb.
#' @param offline If `TRUE`, do not query AnnotationHub.
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
    source = c("auto", "spamtpdb"),
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
#' @param version Optional SpaMTPdb resource version. `NULL` lists every
#'   available external version.
#'
#' @return A data frame describing resources in SpaMTPdb.
#' @export
#'
#' @examples
#' utils::str(formals(SpaMTPDatabaseInfo))
#' SpaMTPDatabaseInfo()
SpaMTPDatabaseInfo <- function(version = NULL) {
  if (!requireNamespace("SpaMTPdb", quietly = TRUE)) {
    stop(
      "SpaMTPDatabaseInfo() requires the SpaMTPdb package.",
      call. = FALSE
    )
  }
  SpaMTPdb::SpaMTPdbResources(version = version)
}
