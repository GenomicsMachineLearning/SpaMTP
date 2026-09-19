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

test_that("corrections restore collapsed parallel interactions and leave other resources untouched", {
  value <- list(list(id = "hsa:04110", protEdges = data.frame(
    src = "RAMP_G_1", dest = "RAMP_G_2", directed = 1L, reaction_type = 4L
  )))
  patch <- list(
    input_sha256 = SpaMTP:::.spamtp_graph_fingerprint(value),
    provenance = list(repair = "fixture"),
    pathways = list(list(protEdges = list(
      rows = c(1L, 1L), directed = c(1L, 2L), reaction_type = c(4L, 6L),
      source_reaction_type = c("Process(inhibition)", "Process(activation)")
    )))
  )
  attr(value, "spamtp_database") <- list(version = "3.0.7")
  result <- SpaMTP:::.spamtp_apply_interaction_patch(value, patch)
  expect_identical(result[[1]]$protEdges$reaction_type, c(4L, 6L))
  expect_identical(result[[1]]$protEdges$src, rep("RAMP_G_1", 2))
  expect_identical(attr(result, "spamtp_database"), attr(value, "spamtp_database"))
  expect_identical(SpaMTP:::.spamtp_apply_interaction_patch(result, patch), result)
  modified <- value
  modified[[1]]$protEdges$dest <- "CUSTOM"
  expect_identical(SpaMTP:::.spamtp_apply_interaction_patch(modified, patch), modified)
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

test_that("the Cell Cycle correction recovers every archived source label offline", {
  fixture <- readRDS(testthat::test_path("fixtures", "cell_cycle_interactions.rds"))
  patch <- readRDS(system.file("extdata", "ramp_kegg_interactions_v1.rds", package = "SpaMTP"))
  # Exercise the real correction on a small, unmodified excerpt. The staged
  # resource test below separately verifies the full-resource fingerprint.
  value <- list(fixture$legacy_topology)
  patch$input_sha256 <- SpaMTP:::.spamtp_graph_fingerprint(value)
  patch$pathways <- patch$pathways[fixture$pathway_index]
  fixed <- SpaMTP:::.spamtp_apply_interaction_patch(value, patch)[[1L]]
  expect_true(all(fixture$legacy_topology$protEdges$reaction_type == 4L))
  expect_identical(fixed$protEdges$source_reaction_type, fixture$source_reaction_type)
  expect_identical(fixed$protEdges$directed,
                   match(fixture$source_direction, c("directed", "undirected")))
  expect_length(unique(fixed$protEdges$reaction_type), 11L)
  expect_equal(sum(fixed$protEdges$reaction_type == 4L), 286L)
  expect_equal(sum(fixed$protEdges$reaction_type == 6L), 113L)
})

test_that("Cell Cycle is repaired when the pinned resource is staged", {
  staging <- Sys.getenv("SPAMTPDB_RESOURCE_DIR", "")
  path <- file.path(staging, "ramp_kegg.rds")
  if (!file.exists(path)) path <- file.path(staging, "3.0.7", "ramp_kegg.rds")
  skip_if_not(nzchar(staging) && file.exists(path), "Pinned external graph is not staged")
  original <- readRDS(path)
  fixed <- SpaMTP:::.spamtp_db_bundle("ramp_kegg", database = list(ramp_kegg = original))$ramp_kegg
  cell <- fixed[[which(vapply(fixed, function(p) p$id == "hsa:04110", logical(1)))]]
  expect_length(unique(cell$protEdges$reaction_type), 11L)
  expect_equal(sum(cell$protEdges$reaction_type == 4L), 286L)
  expect_equal(sum(cell$protEdges$reaction_type == 6L), 113L)
  expect_equal(nrow(cell$protEdges), 1009L)
  expect_identical(attr(fixed, "spamtp_interaction_repair")$graphite_archive, 19L)
  expect_identical(SpaMTP:::.spamtp_repair_pathway_interactions(fixed, "ramp_kegg"), fixed)
})
