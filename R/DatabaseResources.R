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

.spamtp_db_load_bundled <- function(resource) {
  legacy_name <- unname(.spamtp_db_legacy_names[[resource]])
  package_namespace <- tryCatch(asNamespace("SpaMTP"), error = function(e) NULL)
  if (!is.null(package_namespace)) {
    value <- get0(legacy_name, envir = package_namespace, inherits = FALSE)
    if (!is.null(value)) return(value)
  }

  data_environment <- new.env(parent = emptyenv())
  suppressWarnings(utils::data(
    list = legacy_name,
    package = "SpaMTP",
    envir = data_environment
  ))
  value <- get0(legacy_name, envir = data_environment, inherits = FALSE)
  if (is.null(value)) {
    stop(
      "Bundled compatibility resource '", legacy_name,
      "' is unavailable. Install SpaMTPdb or configure its local resource ",
      "directory.",
      call. = FALSE
    )
  }
  value
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
                                source = c("auto", "spamtpdb", "bundled"),
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

  value <- NULL
  resolved_source <- source
  external_error <- NULL
  use_external <- source != "bundled" &&
    requireNamespace("SpaMTPdb", quietly = TRUE)

  if (source == "spamtpdb" && !use_external) {
    stop(
      "source = 'spamtpdb' requires the SpaMTPdb package.",
      call. = FALSE
    )
  }
  if (use_external) {
    value <- tryCatch(
      SpaMTPdb::SpaMTPdbResource(
        resource = resource,
        version = version,
        local_dir = local_dir,
        hub = hub,
        offline = offline
      ),
      error = function(e) {
        external_error <<- conditionMessage(e)
        NULL
      }
    )
    resolved_source <- "spamtpdb"
  }

  if (is.null(value)) {
    if (source == "spamtpdb") {
      stop(external_error, call. = FALSE)
    }
    value <- .spamtp_db_load_bundled(resource)
    resolved_source <- "bundled"
    if (!is.null(external_error)) {
      warning(
        "SpaMTPdb could not provide '", resource, "' (", external_error,
        "); using the bundled compatibility copy.",
        call. = FALSE
      )
    }
  }

  attr(value, "spamtp_database") <- list(
    resource = resource,
    version = if (resolved_source == "spamtpdb") {
      as.character(version %||% "latest")
    } else {
      metadata <- if (identical(resource, "ramp_db_metadata")) value else {
        tryCatch(
          .spamtp_db_load_bundled("ramp_db_metadata"),
          error = function(e) NULL
        )
      }
      if (is.list(metadata)) as.character(metadata$ramp_version %||% NA_character_) else NA_character_
    },
    source = resolved_source,
    local_dir = local_dir,
    offline = offline
  )
  assign(key, value, envir = .spamtp_db_cache)
  value
}

.spamtp_db_bundle <- function(resources,
                              database = NULL,
                              version = "latest",
                              source = c("auto", "spamtpdb", "bundled"),
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
#' Loads a coherent group of annotation resources from [SpaMTPdb] when it is
#' installed. During the database-package transition, `source = "auto"` falls
#' back to compatibility copies bundled with SpaMTP. Retrieved resources are
#' cached for the current R session.
#'
#' @param resources Character vector of resource names. Use
#'   [SpaMTPDatabaseInfo()] to list valid names.
#' @param version SpaMTPdb/RaMP resource version, or `"latest"`.
#' @param source Database source. `"auto"` prefers SpaMTPdb and falls back to
#'   bundled compatibility data; `"spamtpdb"` requires SpaMTPdb; `"bundled"`
#'   explicitly uses the compatibility data.
#' @param local_dir Optional directory containing staged SpaMTPdb `.rds` files.
#' @param hub Optional pre-created `AnnotationHub` object passed to SpaMTPdb.
#' @param offline If `TRUE`, do not query AnnotationHub.
#' @param refresh If `TRUE`, bypass SpaMTP's in-session resource cache.
#'
#' @return A named list containing the requested resources.
#' @export
#'
#' @examples
#' database <- LoadSpaMTPDatabase("ramp_db_metadata")
#' names(database)
LoadSpaMTPDatabase <- function(
    resources = c(
      "chem_props", "source_df", "analyte", "analytehaspathway", "pathway"
    ),
    version = "latest",
    source = c("auto", "spamtpdb", "bundled"),
    local_dir = NULL,
    hub = NULL,
    offline = FALSE,
    refresh = FALSE) {
  .spamtp_db_bundle(
    resources = resources,
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
#' @return A data frame describing resources in SpaMTPdb, or the bundled
#'   compatibility mapping when SpaMTPdb is not installed.
#' @export
#'
#' @examples
#' SpaMTPDatabaseInfo()
SpaMTPDatabaseInfo <- function(version = NULL) {
  if (requireNamespace("SpaMTPdb", quietly = TRUE)) {
    return(SpaMTPdb::SpaMTPdbResources(version = version))
  }
  data.frame(
    resource = names(.spamtp_db_legacy_names),
    legacy_object = unname(.spamtp_db_legacy_names),
    source = "bundled compatibility",
    stringsAsFactors = FALSE
  )
}
