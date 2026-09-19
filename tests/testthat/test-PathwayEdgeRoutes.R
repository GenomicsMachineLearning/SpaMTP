test_that("the viewer separates parallel, reciprocal and self-loop interactions", {
  node <- Sys.which("node")
  skip_if_not(nzchar(node), "Node.js is needed for viewer geometry tests")
  script <- testthat::test_path("..", "javascript", "pathway_edge_routes.test.cjs")
  template <- SpaMTP:::.pn_template_path()
  output <- system2(node, c(shQuote(script), shQuote(template)), stdout = TRUE, stderr = TRUE)
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
  expect_equal(sum(startsWith(output, "PASS ")), 8L)
})
