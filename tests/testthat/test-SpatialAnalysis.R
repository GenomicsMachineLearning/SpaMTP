test_that("FindSpatiallyVariableMetabolites validates sampling controls", {
  expect_error(
    FindSpatiallyVariableMetabolites(NULL, max_spots = 1),
    "max_spots must be NULL or one finite number >= 2",
    fixed = TRUE
  )
  expect_error(
    FindSpatiallyVariableMetabolites(NULL, seed = Inf),
    "seed must be one finite number",
    fixed = TRUE
  )
})
