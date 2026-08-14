test_that("bundled database resources remain available during migration", {
  database <- LoadSpaMTPDatabase(
    "ramp_db_metadata",
    source = "bundled",
    refresh = TRUE
  )

  expect_named(database, "ramp_db_metadata")
  expect_type(database$ramp_db_metadata, "list")
  expect_identical(database$ramp_db_metadata$ramp_version, "3.0.7")
  expect_identical(
    attr(database$ramp_db_metadata, "spamtp_database")$source,
    "bundled"
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
})
