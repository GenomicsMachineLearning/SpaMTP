test_that("custom database resources can be used without a Hub lookup", {
  example_database <- list(
    ramp_db_metadata = list(ramp_version = "example")
  )
  database <- LoadSpaMTPDatabase(
    "ramp_db_metadata",
    database = example_database,
    refresh = TRUE
  )

  expect_named(database, "ramp_db_metadata")
  expect_type(database$ramp_db_metadata, "list")
  expect_identical(database$ramp_db_metadata$ramp_version, "example")
})

test_that("local resources load without a companion package", {
  staging <- tempfile("spamtp-local-")
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  fixture <- list(ramp_version = "3.0.7")
  saveRDS(fixture, file.path(staging, "ramp_db_metadata.rds"))

  database <- LoadSpaMTPDatabase(
    "ramp_db_metadata",
    source = "local",
    local_dir = staging,
    offline = TRUE,
    refresh = TRUE
  )

  expect_identical(database$ramp_db_metadata$ramp_version, "3.0.7")
  expect_identical(
    attr(database$ramp_db_metadata, "spamtp_database")$source,
    "local"
  )
})

test_that("custom resource bundles are validated and subset", {
  custom <- list(
    chem_props = data.frame(exactmass = 100),
    pathway = data.frame(pathwayRampId = "RAMP_P_1")
  )

  selected <- SpaMTP:::.spamtp_db_bundle("chem_props", database = custom)
  expect_named(selected, "chem_props")
  expect_equal(selected$chem_props$exactmass, 100)

  expect_error(
    SpaMTP:::.spamtp_db_bundle("source_df", database = custom),
    "missing resource"
  )
})

test_that("database registry reports canonical resource names", {
  registry <- SpaMTPDatabaseInfo()
  expect_s3_class(registry, "data.frame")
  expect_true("resource" %in% names(registry))
  expect_true(all(
    c("chem_props", "source_df", "analytehaspathway", "pathway") %in%
      registry$resource
  ))
  expect_true(all(registry$source == "bundled"))
  expect_identical(unique(registry$version), "3.0.7")
  expect_identical(SpaMTPDatabaseInfo("latest"), registry)
  expect_identical(SpaMTPDatabaseInfo("3.0.7"), registry)
  expect_equal(nrow(SpaMTPDatabaseInfo("not-bundled")), 0L)
  expect_error(SpaMTPDatabaseInfo(NA_character_), "single non-empty")
  expect_error(SpaMTPDatabaseInfo(c("3.0.7", "latest")), "single non-empty")
})

test_that("default resources are bundled and need no download", {
  before <- loadedNamespaces()
  database <- LoadSpaMTPDatabase("ramp_db_metadata", refresh = TRUE)
  expect_identical(database$ramp_db_metadata$ramp_version, "3.0.7")
  metadata <- attr(database$ramp_db_metadata, "spamtp_database")
  expect_identical(metadata$source, "bundled")
  expect_identical(metadata$version, "3.0.7")
  expect_true(metadata$offline)
  expect_identical(
    LoadSpaMTPDatabase("ramp_db_metadata", source = "bundled", offline = TRUE),
    database
  )
  expect_false(any(c("SpaMTPdb", "SpaMTPData") %in%
                     setdiff(loadedNamespaces(), before)))
  expect_error(LoadSpaMTPDatabase("ramp_db_metadata", version = "9.9"),
               "No bundled resource")
  expect_error(LoadSpaMTPDatabase("ramp_db_metadata", hub = list()), "not a Hub")
})

test_that("local overrides are explicit and version declarations are checked", {
  staging <- tempfile("spamtp-local-")
  dir.create(file.path(staging, "3.0.7"), recursive = TRUE)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  saveRDS(list(ramp_version = "3.0.7"),
          file.path(staging, "3.0.7", "ramp_db_metadata.rds"))
  result <- LoadSpaMTPDatabase("ramp_db_metadata", local_dir = staging)
  expect_identical(attr(result$ramp_db_metadata, "spamtp_database")$source, "local")
  expect_identical(attr(result$ramp_db_metadata, "spamtp_database")$version, "3.0.7")
  expect_error(LoadSpaMTPDatabase("chem_props", local_dir = staging),
               "Local database resource is missing")
  expect_error(LoadSpaMTPDatabase("ramp_db_metadata", source = "local"),
               "requires local_dir")
  expect_error(LoadSpaMTPDatabase("ramp_db_metadata", source = "bundled",
                                  local_dir = staging), "Use source")
  expect_error(LoadSpaMTPDatabase("ramp_db_metadata", version = "9.9",
                                  local_dir = file.path(staging, "3.0.7")),
               "does not match requested version")

  saveRDS(data.frame(mass = 100), file.path(staging, "chem_props.rds"))
  local <- LoadSpaMTPDatabase("chem_props", local_dir = staging)
  expect_identical(attr(local$chem_props, "spamtp_database")$version, "unversioned")
  expect_equal(local$chem_props$mass, 100)
  saveRDS(data.frame(mass = 200), file.path(staging, "chem_props.rds"))
  expect_equal(LoadSpaMTPDatabase("chem_props", local_dir = staging)$chem_props$mass, 100)
  expect_equal(LoadSpaMTPDatabase("chem_props", local_dir = staging,
                                  refresh = TRUE)$chem_props$mass, 200)
})

test_that("bundled graphs retain the Cell Cycle interaction correction", {
  graphs <- LoadSpaMTPDatabase("ramp_kegg")$ramp_kegg
  cell <- graphs[[which(vapply(graphs, function(p) p$id == "hsa:04110", logical(1)))]]
  expect_equal(nrow(cell$protEdges), 1009L)
  expect_length(unique(cell$protEdges$reaction_type), 11L)
  expect_equal(sum(cell$protEdges$reaction_type == 4L), 286L)
  expect_equal(sum(cell$protEdges$reaction_type == 6L), 113L)
  expect_identical(attr(graphs, "spamtp_interaction_repair")$graphite_archive, 19L)
})

test_that("bundled chemistry joins precomputed structure features", {
  database <- LoadSpaMTPDatabase("chem_props")$chem_props
  small <- database[database$chem_source_id == "hmdb:HMDB0000606", , drop = FALSE]
  expect_gt(nrow(small), 0L)
  # A zero runtime limit ensures missing precomputed fields cannot be hidden
  # by the automatic SMILES inference fallback.
  withr::local_options(SpaMTP.max_runtime_smiles = 0L)
  expect_warning(
    normalized <- SpaMTP:::.normalise_metabolite_db(small, collapse_isomers = FALSE),
    NA
  )
  expect_true(all(normalized$structure_valid))
  expect_true(all(normalized$carboxyl_sites == 2L))
})

test_that("standalone installation and website have no companion dependency", {
  description <- utils::packageDescription("SpaMTP")
  dependencies <- unlist(description[c("Depends", "Imports", "Suggests")])
  expect_false(any(grepl("SpaMTPdb|SpaMTPData", dependencies)))
  workflow <- testthat::test_path("..", "..", ".github", "workflows", "pkgdown.yaml")
  skip_if_not(file.exists(workflow), "Website workflow is not installed")
  lines <- readLines(workflow)
  expect_false(any(grepl("SpaMTPdb|SpaMTPData|spamtpdb", lines)))
  expect_true(any(grepl("bundled_database_manifest.csv", lines, fixed = TRUE)))
})

test_that("bundled manifest describes the source-package files", {
  registry <- SpaMTPDatabaseInfo()
  paths <- testthat::test_path("..", "..", registry$rdata_path)
  skip_if_not(all(file.exists(paths)), "Source data files are not installed")
  expect_identical(unname(tools::md5sum(paths)), registry$md5)
  expect_equal(unname(file.info(paths)$size), registry$serialized_bytes)
  expect_setequal(registry$resource, names(SpaMTP:::.spamtp_db_legacy_names))
})
