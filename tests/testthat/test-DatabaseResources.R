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

test_that("SpaMTPdb resources load from an offline staging directory", {
  skip_if_not_installed("SpaMTPdb")
  staging <- tempfile("spamtpdb-")
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  fixture <- list(ramp_version = "3.0.7")
  saveRDS(fixture, file.path(staging, "ramp_db_metadata.rds"))

  database <- LoadSpaMTPDatabase(
    "ramp_db_metadata",
    source = "spamtpdb",
    local_dir = staging,
    offline = TRUE,
    refresh = TRUE
  )

  expect_identical(database$ramp_db_metadata$ramp_version, "3.0.7")
  expect_identical(
    attr(database$ramp_db_metadata, "spamtp_database")$source,
    "spamtpdb"
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
  skip_if_not_installed("SpaMTPdb")
  registry <- SpaMTPDatabaseInfo()
  expect_s3_class(registry, "data.frame")
  expect_true("resource" %in% names(registry))
  expect_true(all(
    c("chem_props", "source_df", "analytehaspathway", "pathway") %in%
      registry$resource
  ))
})
