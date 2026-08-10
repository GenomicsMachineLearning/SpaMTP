test_that("MultiOmicIntegration restores future connection diagnostics", {
  previous <- getOption("future.connections.onMisuse")

  suppressWarnings(
    try(
      MultiOmicIntegration(
        NULL,
        reduction.list = list("missing-a", "missing-b"),
        dims.list = list(1, 1)
      ),
      silent = TRUE
    )
  )

  expect_identical(getOption("future.connections.onMisuse"), previous)
})
