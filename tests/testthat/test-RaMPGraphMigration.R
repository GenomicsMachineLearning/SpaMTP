test_that("RaMP graph migration expands entities split by a newer release", {
  utility <- testthat::test_path(
    "..", "..", "data-raw", "ramp_graph_utils.R"
  )
  skip_if_not(file.exists(utility), "data-raw utilities are excluded from built packages")
  environment <- new.env(parent = globalenv())
  sys.source(utility, envir = environment)

  analyte <- data.frame(
    rampId = c("RAMP_C_new1", "RAMP_C_new2"), type = "compound"
  )
  previous <- data.frame(
    sourceId = c("hmdb:A", "chebi:B"), rampId = "RAMP_C_old",
    geneOrCompound = "compound", commonName = "Old merged entity"
  )
  current <- data.frame(
    sourceId = c("hmdb:A", "chebi:B"),
    rampId = c("RAMP_C_new1", "RAMP_C_new2"),
    geneOrCompound = "compound", commonName = c("Entity 1", "Entity 2")
  )
  migration <- environment$.rgu_build_id_migration(
    "RAMP_C_old", analyte, previous, current,
    manual_overrides = data.frame(
      old_id = character(), new_id = character(), evidence = character()
    )
  )

  expect_setequal(migration$new_id, c("RAMP_C_new1", "RAMP_C_new2"))
  expect_true(all(migration$split))

  edge <- data.frame(
    src = "RAMP_C_old", dest = "RAMP_C_fixed", directed = 1L,
    reaction_type = 1L
  )
  updated <- environment$.rgu_update_edge(
    edge, split(migration$new_id, migration$old_id)
  )
  expect_setequal(updated$src, c("RAMP_C_new1", "RAMP_C_new2"))
  expect_true(all(updated$dest == "RAMP_C_fixed"))
})
