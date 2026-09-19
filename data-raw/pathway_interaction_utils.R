# Source labels, never factor level numbers, define the shared style codes.
.pi_reaction_codes <- function(labels) {
  labels <- as.character(labels)
  known <- c(
    "Process(indirect)" = 1L, "Control(indirect)" = 1L,
    "Binding" = 2L, "Process(missing)" = 3L,
    "Process(inhibition)" = 4L, "Process(phosphorylation)" = 5L,
    "Process(activation)" = 6L, "Process(expression)" = 7L,
    "Process(missing interaction)" = 8L, "Process(indirect effect)" = 9L,
    "Process(binding/association)" = 10L, "Process(ubiquitination)" = 11L,
    "Process(dephosphorylation)" = 12L, "Process(dissociation)" = 13L,
    "Process(methylation)" = 14L, "Process(repression)" = 15L,
    "Process" = 16L, "Process(BiochemicalReaction)" = 16L,
    "Control(In)" = 16L, "Control(Out)" = 16L,
    "Process(glycosylation)" = 17L, "Process(state change)" = 18L
  )
  result <- unname(known[labels])
  result[grepl("^Control\\((In|Out): INHIBITION of ", labels)] <- 4L
  result[grepl("^Control\\((In|Out): ACTIVATION of ", labels)] <- 6L
  unknown <- !is.na(labels) & is.na(result)
  if (any(unknown)) stop("Unmapped source interaction: ", paste(unique(labels[unknown]), collapse = ", "))
  result
}

.pi_direction_codes <- function(labels) {
  labels <- as.character(labels)
  result <- match(labels, c("directed", "undirected"))
  if (any(!is.na(labels) & is.na(result))) stop("Unknown source edge direction.")
  result
}

.pi_convert_edges <- function(edge, identifiers, legacy = FALSE) {
  namespaces <- c(
    ENZYME = "EN", ENTREZID = "entrez", UNIPROT = "uniprot",
    ENSEMBL = "ensembl", SYMBOL = "gene_symbol", PFAM = "PF",
    CHEBI = "chebi", REFSEQ = "REFSQ", HMDB = "hmdb", KEGGCOMP = "kegg",
    CAS = "CAS", PUBCHEM = "pubchem", LIPIDMAPS = "LIPIDMAPS",
    KEGGGLYCAN = "kegg_glycan", CHEMSPIDER = "chemspider", KEGGDRUG = "keggdrug"
  )
  endpoint <- function(side) {
    prefix <- unname(namespaces[as.character(edge[[paste0(side, "_type")]])])
    keys <- paste0(prefix, ":", edge[[side]])
    unname(unlist(mget(keys, identifiers, ifnotfound = list(NA_character_)), use.names = FALSE))
  }
  if (!nrow(edge)) {
    # The historical resources used an all-NA sentinel for empty tables.
    out <- data.frame(src = NA_character_, dest = NA_character_,
                      directed = NA_integer_, reaction_type = NA_integer_)
    if (!legacy) out$source_reaction_type <- NA_character_
    return(out)
  }
  if (legacy) {
    # Reproduce the historical bug ONLY to verify the patch's input exactly.
    direction <- rep(as.integer(edge$direction)[1L], nrow(edge))
    reaction <- rep(as.integer(edge$type)[1L], nrow(edge))
  } else {
    direction <- .pi_direction_codes(edge$direction)
    reaction <- .pi_reaction_codes(edge$type)
  }
  out <- data.frame(src = endpoint("src"), dest = endpoint("dest"),
                    directed = direction, reaction_type = reaction)
  if (!legacy) out$source_reaction_type <- as.character(edge$type)
  out
}

.pi_edge_key <- function(edge) paste(edge$src, edge$dest, sep = "\r")

.pi_make_patch <- function(original, corrected) {
  rows <- match(.pi_edge_key(corrected), .pi_edge_key(original))
  if (anyNA(rows)) stop("A corrected edge has no endpoint match in the original graph.")
  list(rows = rows, directed = corrected$directed,
       reaction_type = corrected$reaction_type,
       source_reaction_type = corrected$source_reaction_type)
}
