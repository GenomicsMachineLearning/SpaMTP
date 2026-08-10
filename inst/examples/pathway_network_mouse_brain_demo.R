# SpaMTP pathway-network demo using the real Mouse Brain/DHB vignette data.
#
# From the SpaMTP repository root:
#   devtools::install(dependencies = FALSE)
#   source("inst/examples/pathway_network_mouse_brain_demo.R")
#
# From an installed SpaMTP package:
#   source(system.file(
#     "examples/pathway_network_mouse_brain_demo.R", package = "SpaMTP"
#   ))

# Avoid socket-based parallel backends on restricted HPC compute nodes. There
# are only two regions in this demo, so a deterministic serial backend is
# sufficient. Remove this block if you have configured another backend.
suppressPackageStartupMessages(library(BiocParallel))
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
devtools::load_all("/vast/scratch/users/lu.t/SpaMTP")
suppressPackageStartupMessages({
  library(Seurat)
})

if (!"max_nodes" %in% names(formals(SpaMTP::PathwayNetworkPlots))) {
  stop("This demo requires the developmental SpaMTP PathwayNetworkPlots().")
}

repo_object <- file.path(
  "vignettes", "vignette_data_files", "Mouse_Brain", "DHB_data",
  "striatum.dhb.data.RDS"
)
if (file.exists(repo_object)) {
  object_file <- repo_object
  output_dir <- file.path(dirname(object_file), "pathway_network_demo_output")
} else {
  download_dir <- file.path(tempdir(), "SpaMTP_mouse_brain_demo")
  dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
  object_file <- file.path(download_dir, "striatum.dhb.data.RDS")
  output_dir <- file.path(getwd(), "pathway_network_demo_output")
  expected_min_bytes <- 60 * 1024^2
  object_is_complete <- file.exists(object_file) &&
    isTRUE(file.info(object_file)$size >= expected_min_bytes)
  if (!object_is_complete) {
    message("Downloading the Mouse Brain/DHB demo object from Zenodo ...")
    download_file <- paste0(object_file, ".part")
    options(timeout = max(600, getOption("timeout", 60)))
    download.file(
      paste0(
        "https://zenodo.org/records/17246900/files/",
        "striatum.dhb.data.RDS?download=1"
      ),
      download_file,
      mode = "wb",
      quiet = FALSE
    )
    if (!file.exists(download_file) ||
        file.info(download_file)$size < expected_min_bytes) {
      stop("The Zenodo download is incomplete: ", download_file)
    }
    if (!file.rename(download_file, object_file)) {
      stop("Could not move the completed download to: ", object_file)
    }
  }
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
marker_cache_file <- file.path(output_dir, "differential_results.rds")
cache_file <- file.path(output_dir, "regional_pathway_results.rds")
annotation_cache_file <- file.path(
  output_dir, "indexed_ramp_metabolite_annotations.rds"
)

message("Loading the 876-spot Mouse Brain/DHB demo object ...")
striatum <- readRDS(object_file)
striatum <- SeuratObject::`Idents<-`(
  striatum,
  value = "striatum"
)

# The published demo object contains a legacy mass-only @tools$db_3 result.
# Re-annotate its SPM features with the current indexed/scored RaMP pipeline.
# The cache stores candidates plus provenance, not a second Seurat object.
current_ramp_version <- as.character(SpaMTP::ramp_db_metadata$ramp_version)
analysis_score_threshold <- 0.05
frontend_score_floor <- 0.01
metabolite_de_pval_cutoff <- 0.05
cached_annotation <- if (file.exists(annotation_cache_file)) {
  readRDS(annotation_cache_file)
} else NULL
cache_is_current <- is.list(cached_annotation) &&
  is.list(cached_annotation$metadata) &&
  identical(cached_annotation$metadata$engine, "indexed-chemical-v2") &&
  identical(cached_annotation$metadata$ramp_version, current_ramp_version) &&
  identical(cached_annotation$metadata$min_score, 0) &&
  is.data.frame(cached_annotation$results)

if (cache_is_current) {
  message("Using cached indexed RaMP metabolite annotations ...")
  striatum@tools$mz_annotation <- cached_annotation
  striatum@tools$db_3 <- cached_annotation$results
} else {
  message("Annotating SPM features with the current indexed RaMP pipeline ...")
  spm_counts <- SeuratObject::LayerData(striatum, assay = "SPM", layer = "counts")
  aggregate_spectrum <- data.frame(
    mz = as.numeric(striatum[["SPM"]][[]]$raw_mz),
    intensity = Matrix::rowMeans(spm_counts)
  )
  striatum <- SpaMTP::AnnotateSM(
    striatum,
    db = NULL,
    assay = "SPM",
    raw.mz.column = "raw_mz",
    ppm_error = 5,
    polarity = "positive",
    ms1_spectrum = aggregate_spectrum,
    # Store every ppm-valid candidate. Pathway functions apply a user-tunable
    # score threshold below, so changing it does not require re-annotation.
    min_score = 0,
    return.only.annotated = FALSE,
    save.intermediate = TRUE,
    verbose = TRUE
  )
  saveRDS(striatum@tools$mz_annotation, annotation_cache_file, compress = "xz")
}
annotation_info <- SpaMTP::AnnotationInfo(striatum)
annotation_signature <- paste(
  annotation_info$engine,
  annotation_info$ramp_version,
  annotation_info$polarity,
  annotation_info$ppm,
  annotation_info$min_score,
  analysis_score_threshold,
  metabolite_de_pval_cutoff,
  annotation_info$candidates,
  sep = "|"
)
message(
  "Annotation source: ", annotation_info$engine,
  "; RaMP ", annotation_info$ramp_version,
  "; candidates ", annotation_info$candidates
)

# Differential analysis and pathway enrichment are cached because they are the
# expensive part. A cache made from an older annotation pipeline is ignored.
cached_results <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
if (is.list(cached_results) &&
    identical(cached_results$annotation_signature, annotation_signature)) {
  message("Using cached differential and regional-pathway results ...")
  results <- cached_results
} else {
  if (file.exists(marker_cache_file)) {
    message("Using cached differential results ...")
    de_list <- readRDS(marker_cache_file)
    deg <- de_list$genes
    dem <- de_list$metabolites
  } else {
    message("Finding transcript markers ...")
    deg <- Seurat::FindAllMarkers(
      striatum,
      test.use = "wilcox_limma",
      assay = "SPT",
      only.pos = FALSE,
      verbose = FALSE
    )

    message("Finding metabolite markers ...")
    dem <- Seurat::FindAllMarkers(
      striatum,
      assay = "SPM",
      only.pos = FALSE,
      verbose = FALSE
    )

    de_list <- list(genes = deg, metabolites = dem)
    saveRDS(de_list, marker_cache_file, compress = "xz")
  }

  message("Running multi-omic regional pathway enrichment ...")
  regpathway <- SpaMTP::FindRegionalPathways(
    striatum,
    analyte_types = c("genes", "metabolites"),
    SM_slot = "counts",
    ST_slot = "counts",
    ident = "striatum",
    DE.list = de_list,
    annotation_score_threshold = analysis_score_threshold,
    # Change to 1 to construct metabolite ranks without a DE-significance
    # cutoff. The coverage network below does not require this recomputation.
    pval_cutoff_mets = metabolite_de_pval_cutoff,
    verbose = TRUE
  )

  results <- list(
    deg = deg,
    dem = dem,
    de_list = de_list,
    regpathway = regpathway,
    annotation_signature = annotation_signature,
    annotation_metadata = annotation_info
  )
  saveRDS(results, cache_file, compress = "xz")
}

# Prefer pathways highlighted in the published vignette. If none survive the
# current RaMP/database analysis, PathwayNetworkPlots selects the top four by
# summed absolute NES instead.
preferred_pathways <- c(
  "Dopamine beta-hydroxylase deficiency",
  "Tryptophan metabolism",
  "Dopamine metabolism",
  "Sphingolipid metabolism: integrated pathway"
)
selected_pathways <- intersect(
  preferred_pathways,
  unique(as.character(results$regpathway$pathwayName))
)
if (!length(selected_pathways)) selected_pathways <- NULL

message("Rendering the strict GSEA leading-edge network ...")
strict_html_file <- SpaMTP::PathwayNetworkPlots(
  striatum,
  ident = "striatum",
  regpathway = results$regpathway,
  DE.list = results$de_list,
  selected_pathways = selected_pathways,
  path = output_dir,
  SM_slot = "counts",
  ST_slot = "counts",
  analyte_types = c("genes", "metabolites"),
  annotation_score_threshold = analysis_score_threshold,
  annotation_score_floor = frontend_score_floor,
  metabolite_detection = "leading_edge",
  image = "slice1",
  top_n_pathways = 4,
  max_nodes = 350,
  label_mode = "detected",
  layout_mode = "repulsion",
  max_spatial_points = 10000,
  verbose = TRUE
)

Sys.sleep(1)
message("Rendering the annotation-coverage network (DE significance not required) ...")
coverage_html_file <- SpaMTP::PathwayNetworkPlots(
  striatum,
  ident = "striatum",
  regpathway = results$regpathway,
  DE.list = results$de_list,
  selected_pathways = selected_pathways,
  path = output_dir,
  SM_slot = "counts",
  ST_slot = "counts",
  analyte_types = c("genes", "metabolites"),
  annotation_score_threshold = analysis_score_threshold,
  annotation_score_floor = frontend_score_floor,
  metabolite_detection = "annotated",
  image = "slice1",
  top_n_pathways = 4,
  max_nodes = 350,
  label_mode = "detected",
  layout_mode = "repulsion",
  max_spatial_points = 10000,
  verbose = TRUE
)

message(
  "\nStrict leading-edge view:\n", normalizePath(strict_html_file),
  "\n\nNon-DE annotation coverage view:\n", normalizePath(coverage_html_file),
  "\n\nBoth HTML files include front-end evidence and score controls (floor = ",
  frontend_score_floor, ")."
)
# In an interactive R session, uncomment the next line:
# browseURL(coverage_html_file)
