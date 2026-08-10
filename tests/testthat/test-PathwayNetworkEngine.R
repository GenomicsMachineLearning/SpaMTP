test_that("topology catalog resolves names and IDs without rescanning lists", {
  topology <- function(id, title) list(id = id, title = title)
  sources <- list(
    wiki = list("Pathway A" = topology("WP1", "Pathway A")),
    hmdb = list("Pathway B" = topology("SMP2", "Pathway B"))
  )
  catalog <- SpaMTP:::.pn_topology_catalog(sources, use_cache = FALSE)
  rows <- data.frame(
    pathwayName = c("pathway a", "unknown title"),
    sourceId = c("wp1", "smp2"),
    type = c("WikiPathways", "HMDB"),
    stringsAsFactors = FALSE
  )
  result <- SpaMTP:::.pn_resolve_topologies(rows, catalog, sources)

  expect_true(all(result$found))
  expect_equal(result$topologies[[1]]$id, "WP1")
  expect_equal(result$topologies[[2]]$id, "SMP2")
})

test_that("edge preparation includes direct protein edges and shared metadata", {
  topology <- list(
    mixedEdges = data.frame(
      src = "RAMP_G_1", dest = "RAMP_C_1", directed = 1L,
      reaction_type = 1L
    ),
    protEdges = data.frame(
      src = "RAMP_G_1", dest = "RAMP_G_2", directed = 1L,
      reaction_type = 2L
    )
  )
  reaction <- data.frame(
    reaction_type = 1:2,
    reaction_name = c("Indirect", "Binding"),
    linetype = c("dashed", "solid"),
    arrowhead = c("arrow", "revarrow"),
    colour = c("black", "blue")
  )
  edges <- SpaMTP:::.pn_prepare_edges(topology, reaction)

  expect_equal(nrow(edges), 2)
  expect_setequal(edges$reaction_name, c("Indirect", "Binding"))
  expect_true(all(edges$weight >= 2L))
})

test_that("node pruning preserves detected analytes and incident neighbours", {
  edges <- data.frame(
    src = paste0("n", 1:20),
    dest = paste0("n", 2:21),
    reaction_type = 1L,
    reaction_name = "Interaction",
    linetype = "solid",
    arrowhead = "arrow",
    colour = "black",
    weight = 2L
  )
  result <- SpaMTP:::.pn_prune_edges(edges, detected = "n20", max_nodes = 6)
  retained <- unique(c(result$src, result$dest))

  expect_true("n20" %in% retained)
  expect_true(any(c("n19", "n21") %in% retained))
  expect_lte(length(retained), 6)
})

test_that("current leading-edge IDs missing from old topology stay visible", {
  edges <- data.frame(
    src = "RAMP_G_1", dest = "RAMP_C_1",
    reaction_type = 1L, reaction_name = "Interaction",
    linetype = "solid", arrowhead = "arrow", colour = "black", weight = 2L,
    stringsAsFactors = FALSE
  )
  metabolite <- data.frame(
    ramp_id = "RAMP_C_NEW", mz_name = "mz-149.04445", Adduct = "M+H",
    common_name = "Current metabolite", avg_log2FC = 1.5,
    stringsAsFactors = FALSE
  )
  result <- SpaMTP:::.pn_build_one_network(
    edges = edges,
    detected = "RAMP_C_NEW",
    cluster_key = "region1",
    gene_lookup = list(),
    metabolite_lookup = list(region1 = metabolite),
    name_map = c(RAMP_G_1 = "GENE1", RAMP_C_1 = "Old metabolite"),
    max_nodes = 10
  )

  node <- result$nodes[result$nodes$id == "RAMP_C_NEW", , drop = FALSE]
  expect_equal(nrow(node), 1L)
  expect_true(node$detected)
  expect_equal(node$name, "Current metabolite")
  expect_match(node$display, "mz-149.04445", fixed = TRUE)
  expect_true(node$leading_edge)
  expect_false(node$annotated)
  expect_true(is.na(node$annotation_score))
})

test_that("annotated pathway mode includes non-DE RaMP compounds", {
  annotations <- data.frame(
    ramp_id = c("RAMP_C_LEADING", "RAMP_C_NONDE", "RAMP_C_OTHER"),
    Score = c(0.8, 0.7, 0.9),
    stringsAsFactors = FALSE
  )
  pathway_rows <- data.frame(
    pathwayName = "Pathway A", sourceId = "WP1",
    stringsAsFactors = FALSE
  )
  membership <- data.frame(
    rampId = c("RAMP_C_LEADING", "RAMP_C_NONDE", "RAMP_C_OTHER"),
    pathwayRampId = c("RAMP_P_1", "RAMP_P_1", "RAMP_P_2"),
    stringsAsFactors = FALSE
  )
  pathways <- data.frame(
    pathwayRampId = c("RAMP_P_1", "RAMP_P_2"),
    pathwayName = c("Pathway A", "Pathway B"),
    sourceId = c("WP1", "WP2"),
    stringsAsFactors = FALSE
  )

  result <- SpaMTP:::.pn_annotated_metabolites_by_pathway(
    annotations, "Pathway A", pathway_rows, membership, pathways
  )

  expect_setequal(
    result[["Pathway A"]],
    c("RAMP_C_LEADING", "RAMP_C_NONDE")
  )

  scores <- SpaMTP:::.pn_annotation_score_map(rbind(
    annotations,
    data.frame(ramp_id = "RAMP_C_NONDE", Score = 0.95)
  ))
  expect_equal(unname(scores["RAMP_C_NONDE"]), 0.95)
})

test_that("network JSON template is valid and contains new controls", {
  payload <- list(
    title = "Test network",
    pathways = "Pathway A",
    clusters = "Cluster 1",
    networks = list(list(list(
      nodes = data.frame(
        id = "RAMP_G_1", group = "rna", expr = 1.2,
        name = "GENE1", display = "GENE1", detected = TRUE
      ),
      links = SpaMTP:::.pn_empty_links()
    ))),
    spatial = list(
      coordinates = data.frame(x = 0.5, y = 0.5),
      clusters = "Cluster 1",
      expression = list(rna = list(RAMP_G_1 = 3), mets = list()),
      displayed_points = 1L,
      total_points = 1L
    ),
    palette = c("#000000", "#ffffff"),
    fc_limit = 2,
    label_mode = "detected",
    metadata = list(
      max_nodes = 500L,
      annotation_score_threshold = 0.05,
      annotation_score_floor = 0.01,
      annotation_score_ceiling = 1,
      metabolite_detection = "annotated",
      annotation_controls = TRUE
    )
  )
  html <- SpaMTP:::.pn_render_html(payload)

  expect_false(grepl("__SPAMTP_PAYLOAD__", html, fixed = TRUE))
  expect_match(html, "Focused network", fixed = TRUE)
  expect_match(html, "Metabolites + detected", fixed = TRUE)
  expect_match(html, "Balanced force", fixed = TRUE)
  expect_match(html, "Gene ↔ metabolite", fixed = TRUE)
  expect_match(html, "d3.forceSimulation", fixed = TRUE)
  expect_match(html, "layoutParameters", fixed = TRUE)
  expect_match(html, "fitNetworkView", fixed = TRUE)
  expect_match(html, "Min annotation score", fixed = TRUE)
  expect_match(html, "All annotations", fixed = TRUE)
  expect_match(html, "visibleWorldBounds", fixed = TRUE)
  expect_match(html, "invertX", fixed = TRUE)
  expect_match(html, "updateSimulationViewport", fixed = TRUE)
  expect_match(html, "Dragged nodes remain pinned", fixed = TRUE)
  expect_match(html, "item.pinned = true", fixed = TRUE)
  expect_match(html, "updateReactionLegend(currentRawNetwork().links", fixed = TRUE)
  expect_match(html, "node-gene, .node-metabolite", fixed = TRUE)
  expect_match(html, "stroke: #64748b", fixed = TRUE)
  expect_match(html, "id=\"spatial-legend\"", fixed = TRUE)
  expect_match(html, "id=\"spatial-gradient\"", fixed = TRUE)
  expect_match(html, "Colours are clipped at the 99th percentile", fixed = TRUE)
  expect_match(html, '"RAMP_G_1"', fixed = TRUE)
  expect_match(html, '"pathways":["Pathway A"]', fixed = TRUE)
  expect_match(html, '"networks":[[', fixed = TRUE)
})

test_that("pathway network API exposes interactive layout selection", {
  layouts <- eval(formals(PathwayNetworkPlots)$layout_mode)

  expect_equal(layouts, c("repulsion", "force", "radial", "bipartite"))
  expect_equal(eval(formals(PathwayNetworkPlots)$annotation_score_floor), 0.01)
})

test_that("standard and curated FMP10 annotation engines are current", {
  expect_true(SpaMTP:::.pn_is_current_annotation_metadata(
    list(engine = "indexed-chemical-v2")
  ))
  expect_true(SpaMTP:::.pn_is_current_annotation_metadata(
    list(engine = "curated-fmp10-v2")
  ))
  expect_false(SpaMTP:::.pn_is_current_annotation_metadata(
    list(engine = "legacy-mass-match")
  ))
  expect_false(SpaMTP:::.pn_is_current_annotation_metadata(list()))
})

test_that("network nodes retain low-score candidates for frontend filtering", {
  edges <- data.frame(
    src = "RAMP_G_1", dest = "RAMP_C_LOW",
    reaction_type = 1L, reaction_name = "Interaction",
    linetype = "solid", arrowhead = "arrow", colour = "black", weight = 2L,
    stringsAsFactors = FALSE
  )
  result <- SpaMTP:::.pn_build_one_network(
    edges = edges,
    detected = "RAMP_G_1",
    cluster_key = "region1",
    gene_lookup = list(),
    metabolite_lookup = list(),
    name_map = c(RAMP_G_1 = "GENE1", RAMP_C_LOW = "Low score"),
    max_nodes = 10,
    leading_edge = "RAMP_G_1",
    annotated = "RAMP_C_LOW",
    annotation_scores = c(RAMP_C_LOW = 0.02)
  )

  node <- result$nodes[result$nodes$id == "RAMP_C_LOW", , drop = FALSE]
  expect_true(node$annotated)
  expect_false(node$detected)
  expect_equal(node$annotation_score, 0.02)
})

test_that("graph node names fall back to chemical properties", {
  source <- data.frame(
    rampId = c("RAMP_C_1", "RAMP_G_1", "RAMP_C_2"),
    sourceId = c("chebi:1", "gene_symbol:GENE1", "chebi:2"),
    IDtype = c("chebi", "gene_symbol", "chebi"),
    commonName = c(NA_character_, "GENE1", NA_character_)
  )
  chemical <- data.frame(
    ramp_id = "RAMP_C_1", common_name = "Metabolite 1"
  )

  result <- SpaMTP:::.pn_name_map(
    source, c("RAMP_C_1", "RAMP_G_1", "RAMP_C_2"), chemical
  )

  expect_equal(unname(result[c("RAMP_C_1", "RAMP_G_1", "RAMP_C_2")]),
               c("Metabolite 1", "GENE1", "chebi:2"))
})

test_that("matrix extraction aligns requested cells without transposing", {
  matrix <- Matrix::Matrix(
    matrix(1:12, nrow = 3, dimnames = list(
      c("geneA", "geneB", "geneC"), c("cell1", "cell2", "cell3", "cell4")
    )),
    sparse = TRUE
  )
  result <- SpaMTP:::.pn_matrix_expression(
    matrix, c("GENEB", "missing"), c("cell4", "cell2")
  )

  expect_equal(result[[1]], c(11, 5))
  expect_null(result[[2]])
})
