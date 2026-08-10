current_annotation_fixture <- function() {
  data.frame(
    observed_mz = 149.04445,
    Adduct = "M+H",
    Ramp_IDs = "RAMP_C_000219574",
    IsomerNames = "legacy display name",
    Score = 0.95,
    MassScore = 0.99,
    ChemicalScore = 1,
    IsotopeScore = NA_real_,
    AdductNetworkScore = NA_real_,
    Error = 0.2,
    stringsAsFactors = FALSE
  )
}

annotation_object_fixture <- function() {
  counts <- Matrix::Matrix(
    matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("mz-149.04445", "mz-150"), c("spot1", "spot2"))
    ),
    sparse = TRUE
  )
  SeuratObject::CreateSeuratObject(counts = counts)
}

test_that("current annotation store is preferred over legacy db_3", {
  current <- current_annotation_fixture()
  legacy <- data.frame(
    observed_mz = 149.04445,
    Isomers_IDs = "hmdb:HMDB0000606",
    stringsAsFactors = FALSE
  )
  tools <- list(
    mz_annotation = SpaMTP:::.annotation_store(current),
    db_3 = legacy
  )

  result <- SpaMTP:::.select_mz_annotations(tools, "current")

  expect_equal(result$Ramp_IDs, "RAMP_C_000219574")
  expect_equal(attr(result, "annotation_source"), "mz_annotation")
  expect_equal(
    attr(result, "annotation_metadata")$engine,
    "indexed-chemical-v2"
  )
})

test_that("current mode refuses an unversioned legacy annotation", {
  tools <- list(db_3 = data.frame(
    observed_mz = 149.04445,
    Isomers_IDs = "hmdb:HMDB0000606",
    stringsAsFactors = FALSE
  ))

  expect_error(
    SpaMTP:::.select_mz_annotations(tools, "current"),
    "Run AnnotateSM"
  )
  expect_warning(
    legacy <- SpaMTP:::.select_mz_annotations(tools, "auto"),
    "legacy"
  )
  expect_equal(attr(legacy, "annotation_source"), "db_3-legacy")
})

test_that("current Ramp IDs map directly to current chemical metadata", {
  object <- annotation_object_fixture()
  object@tools$mz_annotation <- SpaMTP:::.annotation_store(
    current_annotation_fixture()
  )
  chemical <- data.frame(
    ramp_id = "RAMP_C_000219574",
    chem_source_id = "hmdb:HMDB0000606",
    common_name = "D-2-Hydroxyglutaric acid",
    stringsAsFactors = FALSE
  )

  result <- SpaMTP:::.resolve_pathway_metabolite_annotations(
    object, annotation_source = "current", chemical_properties = chemical
  )

  expect_equal(result$ramp_id, "RAMP_C_000219574")
  expect_equal(result$chem_source_id, "hmdb:HMDB0000606")
  expect_equal(result$IsomerNames, "D-2-Hydroxyglutaric acid")
  expect_equal(result$Score, 0.95)
})

test_that("pathway resolver keeps the best scored adduct hypothesis", {
  object <- annotation_object_fixture()
  annotations <- current_annotation_fixture()[c(1, 1), ]
  annotations$Adduct <- c("M+Na", "M+H")
  annotations$Score <- c(0.2, 0.9)
  annotations$Error <- c(4.5, 0.4)
  object@tools$mz_annotation <- SpaMTP:::.annotation_store(annotations)
  chemical <- data.frame(
    ramp_id = "RAMP_C_000219574",
    chem_source_id = "hmdb:HMDB0000606",
    common_name = "D-2-Hydroxyglutaric acid",
    stringsAsFactors = FALSE
  )

  result <- SpaMTP:::.resolve_pathway_metabolite_annotations(
    object, annotation_source = "current", chemical_properties = chemical
  )

  expect_equal(nrow(result), 1L)
  expect_equal(result$Adduct, "M+H")
  expect_equal(result$Score, 0.9)
})

test_that("pathway annotation score threshold is user adjustable", {
  object <- annotation_object_fixture()
  annotations <- current_annotation_fixture()[c(1, 1), ]
  annotations$observed_mz <- c(149.04445, 150.00000)
  annotations$Ramp_IDs <- c("RAMP_C_HIGH", "RAMP_C_LOW")
  annotations$Score <- c(0.8, 0.02)
  object@tools$mz_annotation <- SpaMTP:::.annotation_store(annotations)
  chemical <- data.frame(
    ramp_id = c("RAMP_C_HIGH", "RAMP_C_LOW"),
    chem_source_id = c("chebi:1", "chebi:2"),
    common_name = c("high score", "low score"),
    stringsAsFactors = FALSE
  )

  strict <- SpaMTP:::.resolve_pathway_metabolite_annotations(
    object, score_threshold = 0.05, chemical_properties = chemical
  )
  permissive <- SpaMTP:::.resolve_pathway_metabolite_annotations(
    object, score_threshold = 0.01, chemical_properties = chemical
  )

  expect_equal(strict$ramp_id, "RAMP_C_HIGH")
  expect_setequal(permissive$ramp_id, c("RAMP_C_HIGH", "RAMP_C_LOW"))
  expect_equal(
    attr(strict, "annotation_metadata")$analysis_score_threshold,
    0.05
  )
})

test_that("AnnotationInfo reports legacy objects without treating them as current", {
  object <- annotation_object_fixture()
  object@tools$db_3 <- data.frame(
    observed_mz = 100,
    Isomers_IDs = "chebi:1",
    stringsAsFactors = FALSE
  )

  info <- AnnotationInfo(object)

  expect_equal(info$schema_version, 1L)
  expect_equal(info$engine, "legacy-mass-match")
})

test_that("curated FMP10 annotations are mapped into the current RaMP schema", {
  fmp <- data.frame(
    mass = c(179.0831, 200),
    observed_mz = c(179.0832, 200.001),
    Adduct = c("[M+K]", "[M+FMP10]+"),
    Formula = c("C9H16O", "C5H10O5"),
    Isomers = c("LMFA05000118", "custom"),
    Isomers_IDs = c("LIPIDMAPS:LMFA05000118", "HMDB:HMDBTEST"),
    IsomerNames = c("test lipid", "test metabolite"),
    stringsAsFactors = FALSE
  )
  chemical <- data.frame(
    ramp_id = c("RAMP_C_000000123", "RAMP_C_000000456"),
    chem_source_id = c("lipidmaps:LMFA05000118", "hmdb:HMDBTEST"),
    stringsAsFactors = FALSE
  )

  result <- SpaMTP:::.fmp10_current_annotations(
    fmp,
    mass_threshold = 0.05,
    chemical_properties = chemical
  )
  store <- SpaMTP:::.annotation_store(
    result,
    metadata = list(engine = "curated-fmp10-v2")
  )

  expect_true(SpaMTP:::.annotation_has_current_schema(result))
  expect_setequal(
    result$Ramp_IDs,
    c("RAMP_C_000000123", "RAMP_C_000000456")
  )
  expect_true(all(result$Score > 0 & result$Score <= 1))
  expect_equal(store$metadata$engine, "curated-fmp10-v2")
})

test_that("annotation search falls back to current stored labels", {
  object <- annotation_object_fixture()
  object@tools$mz_annotation <- SpaMTP:::.annotation_store(
    current_annotation_fixture()
  )

  result <- SearchAnnotations(object, "display name", assay = "RNA")

  expect_equal(nrow(result), 1L)
  expect_equal(result$mz_names, "mz-149.04445")
})
