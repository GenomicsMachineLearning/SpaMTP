test_that("legacy annotation IDs expand independently of mismatched names", {
  annotations <- data.frame(
    observed_mz = 130.0504,
    Isomers_IDs = "hmdb:HMDB0000267; hmdb:HMDB0001369",
    IsomerNames = "Pyroglutamic acid; extra alias; another alias",
    stringsAsFactors = FALSE
  )

  expanded <- .expand_pathway_annotation_ids(annotations)

  expect_equal(
    expanded$Isomers_IDs,
    c("hmdb:HMDB0000267", "hmdb:HMDB0001369")
  )
  expect_equal(nrow(expanded), 2L)
  expect_true(all(expanded$IsomerNames == annotations$IsomerNames))
})

test_that("legacy annotation expansion removes empty identifiers", {
  annotations <- data.frame(
    Isomers_IDs = c("chebi:123; ; chebi:456", NA_character_),
    IsomerNames = c("one; two", NA_character_),
    stringsAsFactors = FALSE
  )

  expanded <- .expand_pathway_annotation_ids(annotations)

  expect_equal(expanded$Isomers_IDs, c("chebi:123", "chebi:456"))
})
