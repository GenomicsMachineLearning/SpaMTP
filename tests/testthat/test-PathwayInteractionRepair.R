interaction_utils <- function() {
  path <- testthat::test_path("..", "..", "data-raw", "pathway_interaction_utils.R")
  skip_if_not(file.exists(path), "Maintainer source utilities are not installed")
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)
  env
}

test_that("conversion preserves every interaction and direction, independent of factor levels", {
  util <- interaction_utils()
  edge <- data.frame(
    src_type = "ENTREZID", src = c("1", "1", "2"),
    dest_type = "ENTREZID", dest = c("2", "3", "3"),
    direction = factor(c("directed", "undirected", "directed"),
                       levels = c("undirected", "directed")),
    type = factor(c("Process(inhibition)", "Process(activation)", "Process(phosphorylation)"),
                  levels = c("Process(phosphorylation)", "Process(inhibition)", "Process(activation)"))
  )
  identifiers <- list2env(list("entrez:1" = "RAMP_G_1", "entrez:2" = "RAMP_G_2",
                              "entrez:3" = "RAMP_G_3"), parent = emptyenv())
  result <- util$.pi_convert_edges(edge, identifiers)
  expect_identical(result$reaction_type, c(4L, 6L, 5L))
  expect_identical(result$directed, c(1L, 2L, 1L))
  expect_identical(result$source_reaction_type, as.character(edge$type))
  expect_identical(result$src, c("RAMP_G_1", "RAMP_G_1", "RAMP_G_2"))
  expect_identical(result$dest, c("RAMP_G_2", "RAMP_G_3", "RAMP_G_3"))
  empty <- util$.pi_convert_edges(edge[FALSE, ], identifiers)
  expect_true(all(is.na(empty)))
})

test_that("database-specific labels share style codes without confusing activation and inhibition", {
  util <- interaction_utils()
  expect_identical(util$.pi_reaction_codes(c(
    "Binding", "Control(In: ACTIVATION of BiochemicalReaction)",
    "Control(Out: INHIBITION of BiochemicalReaction)", "Process(BiochemicalReaction)",
    "Control(indirect)"
  )), c(2L, 6L, 4L, 16L, 1L))
  expect_error(util$.pi_reaction_codes("new unrecognised interaction"), "Unmapped source")
  expect_error(util$.pi_direction_codes("new direction"), "Unknown source")
})

test_that("network labels retain source semantics and undirected edges have no arrow", {
  topology <- list(protEdges = data.frame(
    src = "a", dest = c("b", "c"), directed = c(1L, 2L), reaction_type = c(4L, 2L),
    source_reaction_type = c("Control(Out: INHIBITION of BiochemicalReaction)", "Binding")
  ), mixedEdges = data.frame(src = "c", dest = "d", reaction_type = 6L))
  constants <- new.env(parent = emptyenv())
  utils::data("reaction_type", package = "SpaMTP", envir = constants)
  reaction <- constants$reaction_type
  edges <- SpaMTP:::.pn_prepare_edges(topology, reaction)
  expect_identical(edges$reaction_name[edges$dest == "b"], topology$protEdges$source_reaction_type[1])
  expect_identical(edges$arrowhead[edges$dest == "b"], "thead")
  expect_identical(edges$arrowhead[edges$dest == "c"], "none")
  expect_identical(edges$reaction_name[edges$dest == "d"], "Process(activation)")
})

test_that("bundled Cell Cycle preserves the archived interaction labels", {
  fixture <- readRDS(testthat::test_path("fixtures", "cell_cycle_interactions.rds"))
  constants <- new.env(parent = emptyenv())
  utils::data("RAMP_kegg", package = "SpaMTP", envir = constants)
  fixed <- constants$RAMP_kegg[["Cell cycle"]]
  expect_true(all(fixture$legacy_topology$protEdges$reaction_type == 4L))
  expect_identical(fixed$protEdges$source_reaction_type, fixture$source_reaction_type)
  expect_identical(fixed$protEdges$directed,
    match(fixture$source_direction, c("directed", "undirected")))
  expect_length(unique(fixed$protEdges$reaction_type), 11L)
  expect_equal(sum(fixed$protEdges$reaction_type == 4L), 286L)
  expect_equal(sum(fixed$protEdges$reaction_type == 6L), 113L)
  expect_equal(nrow(fixed$protEdges), 1009L)
})

test_that("all bundled graph collections carry the correction and audited edge counts", {
  objects <- c(RAMP_kegg = 215713L, RAMP_Reactome = 2568277L,
               RAMP_wikipathway = 88284L, RAMP_hmdb = 2868575L)
  fields <- c("protEdges", "protPropEdges", "metabolEdges", "metabolPropEdges", "mixedEdges")
  for (object in names(objects)) {
    constants <- new.env(parent = emptyenv())
    utils::data(list = object, package = "SpaMTP", envir = constants)
    graph <- constants[[object]]
    provenance <- attr(graph, "spamtp_interaction_repair", exact = TRUE)
    expect_identical(provenance$repair, "pathway-interactions-v1")
    expect_identical(provenance$graphite_archive, 19L)
    rows <- sum(vapply(graph, function(p) {
      sum(vapply(p[fields], nrow, integer(1)))
    }, integer(1)))
    expect_equal(rows, unname(objects[[object]]))
  }
})
