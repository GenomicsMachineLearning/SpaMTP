.spamtp_annotation_schema <- 2L
.spamtp_annotation_engine <- "indexed-chemical-v2"

.annotation_ramp_version <- function() {
  metadata <- tryCatch(
    .spamtp_db_resource("ramp_db_metadata"),
    error = function(e) get0("ramp_db_metadata", inherits = TRUE)
  )
  if (is.list(metadata) && length(metadata$ramp_version) == 1L) {
    return(as.character(metadata$ramp_version))
  }
  NA_character_
}

.annotation_has_current_schema <- function(x) {
  required <- c(
    "observed_mz", "Adduct", "Ramp_IDs", "Score", "MassScore",
    "ChemicalScore", "IsotopeScore", "AdductNetworkScore"
  )
  is.data.frame(x) && all(required %in% names(x))
}

.annotation_store <- function(results, metadata = list()) {
  metadata <- utils::modifyList(
    list(
      schema_version = .spamtp_annotation_schema,
      engine = .spamtp_annotation_engine,
      ramp_version = .annotation_ramp_version(),
      generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      candidates = nrow(results)
    ),
    metadata
  )
  list(metadata = metadata, results = results)
}

.select_mz_annotations <- function(tools,
                                   annotation_source = c("current", "auto", "legacy"),
                                   warn_legacy = TRUE) {
  annotation_source <- match.arg(annotation_source)
  current <- tools[["mz_annotation"]]
  if (annotation_source != "legacy" && is.list(current) &&
      .annotation_has_current_schema(current$results)) {
    result <- current$results
    attr(result, "annotation_metadata") <- current$metadata
    attr(result, "annotation_source") <- "mz_annotation"
    return(result)
  }

  compatibility <- tools[["db_3"]]
  if (annotation_source != "legacy" &&
      .annotation_has_current_schema(compatibility)) {
    attr(compatibility, "annotation_metadata") <- list(
      schema_version = .spamtp_annotation_schema,
      engine = .spamtp_annotation_engine,
      ramp_version = .annotation_ramp_version(),
      provenance = "current-compatible db_3"
    )
    attr(compatibility, "annotation_source") <- "db_3-current"
    return(compatibility)
  }

  if (annotation_source == "current") {
    stop(
      "No current indexed RaMP annotation result was found. Run AnnotateSM() ",
      "again (db = chem_props, save.intermediate = TRUE), or explicitly set ",
      "annotation_source = 'legacy' for an older @tools$db_3 result."
    )
  }
  if (!is.data.frame(compatibility)) {
    stop("No metabolite annotation result was found in the SpaMTP object.")
  }
  if (isTRUE(warn_legacy) && !.annotation_has_current_schema(compatibility)) {
    warning(
      "Using a legacy @tools$db_3 annotation without scored RaMP IDs. ",
      "Re-run AnnotateSM() to use the current indexed annotation pipeline.",
      call. = FALSE
    )
  }
  attr(compatibility, "annotation_metadata") <- list(
    schema_version = 1L,
    engine = "legacy-mass-match",
    ramp_version = NA_character_,
    provenance = "legacy db_3"
  )
  attr(compatibility, "annotation_source") <- "db_3-legacy"
  compatibility
}

.get_mz_annotations <- function(SpaMTP,
                                annotation_source = c("current", "auto", "legacy"),
                                warn_legacy = TRUE) {
  .select_mz_annotations(
    SpaMTP@tools,
    annotation_source = annotation_source,
    warn_legacy = warn_legacy
  )
}

.expand_annotation_column <- function(x, column) {
  if (!column %in% names(x)) {
    stop("Annotation result is missing required column: ", column)
  }
  x <- tidyr::separate_rows(
    as.data.frame(x, stringsAsFactors = FALSE),
    dplyr::all_of(column),
    sep = "\\s*;\\s*"
  )
  x[[column]] <- trimws(as.character(x[[column]]))
  x[!is.na(x[[column]]) & nzchar(x[[column]]), , drop = FALSE]
}

.preferred_chemical_metadata <- function(chemical_properties, ramp_ids) {
  if (is.null(chemical_properties) || !length(ramp_ids)) {
    return(data.frame(
      ramp_id = character(), common_name = character(),
      chem_source_id = character(), stringsAsFactors = FALSE
    ))
  }
  required <- c("ramp_id", "chem_source_id", "common_name")
  missing <- setdiff(required, names(chemical_properties))
  if (length(missing)) {
    stop("chem_props is missing required column(s): ", paste(missing, collapse = ", "))
  }
  value <- as.data.frame(
    chemical_properties[chemical_properties$ramp_id %in% ramp_ids, required],
    stringsAsFactors = FALSE
  )
  if (!nrow(value)) return(value)
  source_type <- tolower(sub(":.*$", "", value$chem_source_id))
  source_priority <- match(
    source_type,
    c("hmdb", "chebi", "kegg", "lipidmaps", "refmet", "pubchem", "wikidata")
  )
  source_priority[is.na(source_priority)] <- 99L
  has_name <- !is.na(value$common_name) & nzchar(trimws(value$common_name))
  value <- value[order(value$ramp_id, !has_name, source_priority), , drop = FALSE]
  value <- value[!duplicated(value$ramp_id), , drop = FALSE]
  rownames(value) <- NULL
  value
}

.resolve_pathway_metabolite_annotations <- function(
    SpaMTP,
    annotation_source = c("current", "auto", "legacy"),
    score_threshold = 0,
    chemical_properties = NULL,
    warn_legacy = TRUE) {
  annotation_source <- match.arg(annotation_source)
  if (!is.numeric(score_threshold) || length(score_threshold) != 1L ||
      !is.finite(score_threshold) || score_threshold < 0) {
    stop("score_threshold must be one non-negative finite number.")
  }
  annotations <- .get_mz_annotations(
    SpaMTP, annotation_source = annotation_source, warn_legacy = warn_legacy
  )
  source_used <- attr(annotations, "annotation_source")
  metadata <- attr(annotations, "annotation_metadata")

  if (.annotation_has_current_schema(annotations)) {
    score <- suppressWarnings(as.numeric(annotations$Score))
    annotations <- annotations[
      is.finite(score) & score >= score_threshold,
      ,
      drop = FALSE
    ]
    value <- .expand_annotation_column(annotations, "Ramp_IDs")
    value$ramp_id <- toupper(as.character(value$Ramp_IDs))
    value <- value[grepl("^RAMP_C_", value$ramp_id), , drop = FALSE]
    if (!"mz_name" %in% names(value)) {
      value$mz_name <- paste0("mz-", value$observed_mz)
    }
    if (!"Adduct" %in% names(value)) value$Adduct <- ""
    if (!"IsomerNames" %in% names(value)) value$IsomerNames <- NA_character_

    chemical <- .preferred_chemical_metadata(
      chemical_properties, unique(value$ramp_id)
    )
    if (nrow(chemical)) {
      names(chemical)[names(chemical) == "common_name"] <- ".current_name"
      names(chemical)[names(chemical) == "chem_source_id"] <- ".current_source_id"
      value <- merge(value, chemical, by = "ramp_id", all.x = TRUE, sort = FALSE)
      named <- !is.na(value$.current_name) & nzchar(trimws(value$.current_name))
      value$IsomerNames[named] <- value$.current_name[named]
      value$chem_source_id <- value$.current_source_id
      value$.current_name <- NULL
      value$.current_source_id <- NULL
    } else {
      value$chem_source_id <- NA_character_
    }
    value$inputid <- value$ramp_id
  } else {
    if (score_threshold > 0) {
      warning(
        "score_threshold cannot be applied to a legacy unscored annotation.",
        call. = FALSE
      )
    }
    if (is.null(chemical_properties)) {
      stop("Legacy metabolite annotations require the chem_props lookup table.")
    }
    value <- .expand_annotation_column(annotations, "Isomers_IDs")
    value$chem_source_id <- trimws(as.character(value$Isomers_IDs))
    value$inputid <- value$chem_source_id
    value <- merge(chemical_properties, value, by = "chem_source_id", sort = FALSE)
    if (!"mz_name" %in% names(value)) {
      value$mz_name <- paste0("mz-", value$observed_mz)
    }
    if (!"Adduct" %in% names(value)) value$Adduct <- ""
    if (!"IsomerNames" %in% names(value)) value$IsomerNames <- NA_character_
    canonical <- !is.na(value$common_name) & nzchar(trimws(value$common_name))
    value$IsomerNames[canonical] <- value$common_name[canonical]
  }

  if (!nrow(value)) {
    stop(
      "The selected annotation result contains no current RaMP compound IDs ",
      "at score_threshold = ", score_threshold, "."
    )
  }
  order_columns <- intersect(c("Score", "Error"), names(value))
  if (length(order_columns)) {
    score <- if ("Score" %in% names(value)) -suppressWarnings(as.numeric(value$Score)) else 0
    error <- if ("Error" %in% names(value)) suppressWarnings(as.numeric(value$Error)) else 0
    value <- value[order(value$mz_name, value$ramp_id, score, error, na.last = TRUE), , drop = FALSE]
  }
  # The same neutral metabolite can match one observed feature through several
  # adduct hypotheses. Keep the highest-scoring interpretation so pathway
  # statistics and labels consume the ranking instead of arbitrary row order.
  key <- paste(value$mz_name, value$ramp_id, sep = "\r")
  value <- value[!duplicated(key), , drop = FALSE]
  attr(value, "annotation_source") <- source_used
  if (is.null(metadata)) metadata <- list()
  metadata$analysis_score_threshold <- score_threshold
  attr(value, "annotation_metadata") <- metadata
  rownames(value) <- NULL
  value
}

#' Inspect the metabolite annotation provenance stored in a SpaMTP object
#'
#' Reports whether downstream pathway functions will use the current indexed,
#' scored RaMP annotation pipeline or a legacy mass-match result.
#'
#' @param SpaMTP A SpaMTP Seurat object.
#'
#' @return A named list containing annotation schema, engine, RaMP version,
#'   provenance, and candidate count where available.
#' @export
AnnotationInfo <- function(SpaMTP) {
  current <- SpaMTP@tools[["mz_annotation"]]
  if (is.list(current) && .annotation_has_current_schema(current$results)) {
    return(current$metadata)
  }
  compatibility <- SpaMTP@tools[["db_3"]]
  if (.annotation_has_current_schema(compatibility)) {
    return(list(
      schema_version = .spamtp_annotation_schema,
      engine = .spamtp_annotation_engine,
      ramp_version = .annotation_ramp_version(),
      provenance = "current-compatible db_3",
      candidates = nrow(compatibility)
    ))
  }
  list(
    schema_version = if (is.data.frame(compatibility)) 1L else NA_integer_,
    engine = if (is.data.frame(compatibility)) "legacy-mass-match" else NA_character_,
    ramp_version = NA_character_,
    provenance = if (is.data.frame(compatibility)) "legacy db_3" else "not annotated",
    candidates = if (is.data.frame(compatibility)) nrow(compatibility) else 0L
  )
}
