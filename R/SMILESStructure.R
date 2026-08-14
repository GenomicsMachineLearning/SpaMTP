#### Lightweight SMILES decomposition and ionisation priors #################

.empty_smiles_graph <- function(valid = FALSE, message = NA_character_) {
  list(
    atoms = data.frame(
      element = character(), aromatic = logical(), bracket = logical(),
      explicit_h = integer(), charge = integer(), stringsAsFactors = FALSE
    ),
    edges = data.frame(
      from = integer(), to = integer(), order = numeric(),
      stringsAsFactors = FALSE
    ),
    valid = valid,
    message = message
  )
}

.parse_bracket_atom <- function(token) {
  content <- substring(token, 2L, nchar(token) - 1L)
  hit <- regexec("^[0-9]*([A-Z][a-z]?|[bcnops])", content, perl = TRUE)
  value <- regmatches(content, hit)[[1]]
  if (length(value) < 2L) return(NULL)
  raw_element <- value[[2]]
  aromatic <- raw_element %in% c("b", "c", "n", "o", "p", "s")
  element <- paste0(toupper(substring(raw_element, 1L, 1L)),
                    substring(raw_element, 2L))

  explicit_h <- 0L
  if (element != "H") {
    h_hit <- regexec("H([0-9]*)", content, perl = TRUE)
    h_value <- regmatches(content, h_hit)[[1]]
    if (length(h_value)) {
      explicit_h <- if (length(h_value) < 2L || !nzchar(h_value[[2]])) {
        1L
      } else {
        as.integer(h_value[[2]])
      }
    }
  }

  charge <- 0L
  charge_hits <- gregexpr("[+-][0-9]+|[+-]+", content, perl = TRUE)
  charge_tokens <- regmatches(content, charge_hits)[[1]]
  if (length(charge_tokens) && !identical(charge_tokens, character(0))) {
    charge <- sum(vapply(charge_tokens, function(x) {
      sign <- if (substring(x, 1L, 1L) == "+") 1L else -1L
      digits <- gsub("[^0-9]", "", x)
      magnitude <- if (nzchar(digits)) as.integer(digits) else nchar(x)
      sign * magnitude
    }, integer(1)))
  }
  list(
    element = element, aromatic = aromatic, bracket = TRUE,
    explicit_h = explicit_h, charge = charge
  )
}

.parse_smiles_graph <- function(smiles) {
  if (length(smiles) != 1L || is.na(smiles) || !nzchar(trimws(smiles))) {
    return(.empty_smiles_graph(message = "missing SMILES"))
  }
  smiles <- trimws(as.character(smiles))
  elements <- character()
  aromatic <- logical()
  bracket <- logical()
  explicit_h <- integer()
  charge <- integer()
  edge_from <- integer()
  edge_to <- integer()
  edge_order <- numeric()
  current <- NA_integer_
  branch_stack <- integer()
  rings <- list()
  pending_bond <- NA_real_
  position <- 1L
  n_chars <- nchar(smiles)
  error_message <- NULL

  add_atom <- function(atom) {
    new_id <- length(elements) + 1L
    elements <<- c(elements, atom$element)
    aromatic <<- c(aromatic, atom$aromatic)
    bracket <<- c(bracket, atom$bracket)
    explicit_h <<- c(explicit_h, atom$explicit_h)
    charge <<- c(charge, atom$charge)
    if (!is.na(current)) {
      order <- pending_bond
      if (is.na(order)) {
        order <- if (aromatic[[current]] && atom$aromatic) 1.5 else 1
      }
      edge_from <<- c(edge_from, current)
      edge_to <<- c(edge_to, new_id)
      edge_order <<- c(edge_order, order)
    }
    current <<- new_id
    pending_bond <<- NA_real_
  }

  while (position <= n_chars && is.null(error_message)) {
    char <- substring(smiles, position, position)
    two_chars <- if (position < n_chars) {
      substring(smiles, position, position + 1L)
    } else ""

    if (char == "[") {
      remainder <- substring(smiles, position)
      close <- regexpr("]", remainder, fixed = TRUE)[[1]]
      if (close < 0L) {
        error_message <- "unclosed bracket atom"
        next
      }
      token <- substring(smiles, position, position + close - 1L)
      atom <- .parse_bracket_atom(token)
      if (is.null(atom)) {
        error_message <- paste0("unsupported bracket atom ", token)
        next
      }
      add_atom(atom)
      position <- position + close
      next
    }

    if (two_chars %in% c("Cl", "Br")) {
      add_atom(list(
        element = two_chars, aromatic = FALSE, bracket = FALSE,
        explicit_h = 0L, charge = 0L
      ))
      position <- position + 2L
      next
    }

    if (char %in% c("B", "C", "N", "O", "P", "S", "F", "I",
                    "b", "c", "n", "o", "p", "s")) {
      is_aromatic <- char %in% c("b", "c", "n", "o", "p", "s")
      add_atom(list(
        element = toupper(char), aromatic = is_aromatic, bracket = FALSE,
        explicit_h = 0L, charge = 0L
      ))
      position <- position + 1L
      next
    }

    if (char == "(") {
      if (is.na(current)) error_message <- "branch without a preceding atom" else {
        branch_stack <- c(branch_stack, current)
      }
      position <- position + 1L
      next
    }
    if (char == ")") {
      if (!length(branch_stack)) error_message <- "unmatched closing branch" else {
        current <- branch_stack[[length(branch_stack)]]
        branch_stack <- branch_stack[-length(branch_stack)]
      }
      position <- position + 1L
      next
    }
    if (char %in% c("-", "=", "#", ":")) {
      pending_bond <- c("-" = 1, "=" = 2, "#" = 3, ":" = 1.5)[[char]]
      position <- position + 1L
      next
    }
    if (char %in% c("/", "\\", "@")) {
      position <- position + 1L
      next
    }
    if (char == ".") {
      current <- NA_integer_
      pending_bond <- NA_real_
      position <- position + 1L
      next
    }

    if (grepl("[0-9]", char) || char == "%") {
      if (is.na(current)) {
        error_message <- "ring closure without a preceding atom"
        next
      }
      if (char == "%") {
        if (position + 2L > n_chars) {
          error_message <- "incomplete two-digit ring closure"
          next
        }
        ring_label <- substring(smiles, position + 1L, position + 2L)
        if (!grepl("^[0-9]{2}$", ring_label)) {
          error_message <- "invalid two-digit ring closure"
          next
        }
        position <- position + 3L
      } else {
        ring_label <- char
        position <- position + 1L
      }
      previous <- rings[[ring_label]]
      if (is.null(previous)) {
        rings[[ring_label]] <- list(atom = current, order = pending_bond)
      } else {
        order <- pending_bond
        if (is.na(order)) order <- previous$order
        if (is.na(order)) {
          order <- if (aromatic[[current]] && aromatic[[previous$atom]]) 1.5 else 1
        }
        edge_from <- c(edge_from, previous$atom)
        edge_to <- c(edge_to, current)
        edge_order <- c(edge_order, order)
        rings[[ring_label]] <- NULL
      }
      pending_bond <- NA_real_
      next
    }

    error_message <- paste0("unsupported SMILES token '", char, "'")
  }

  if (is.null(error_message) && length(branch_stack)) {
    error_message <- "unclosed branch"
  }
  if (is.null(error_message) && length(rings)) {
    error_message <- "unclosed ring"
  }
  if (!length(elements) && is.null(error_message)) error_message <- "no atoms"
  if (!is.null(error_message)) {
    return(.empty_smiles_graph(message = error_message))
  }

  list(
    atoms = data.frame(
      element = elements, aromatic = aromatic, bracket = bracket,
      explicit_h = explicit_h, charge = charge, stringsAsFactors = FALSE
    ),
    edges = data.frame(
      from = edge_from, to = edge_to, order = edge_order,
      stringsAsFactors = FALSE
    ),
    valid = TRUE,
    message = NA_character_
  )
}

.smiles_neighbours <- function(graph, atom) {
  edges <- graph$edges
  hit_from <- which(edges$from == atom)
  hit_to <- which(edges$to == atom)
  list(
    atom = c(edges$to[hit_from], edges$from[hit_to]),
    order = c(edges$order[hit_from], edges$order[hit_to])
  )
}

.smiles_kekule_aromatic_atoms <- function(graph) {
  n_atoms <- nrow(graph$atoms)
  perceived <- rep(FALSE, n_atoms)
  if (n_atoms < 6L || sum(graph$edges$order == 2) < 3L) {
    return(perceived)
  }

  allowed <- graph$atoms$element %in% c("C", "N")
  adjacent_atoms <- vector("list", n_atoms)
  adjacent_orders <- vector("list", n_atoms)
  for (i in seq_len(nrow(graph$edges))) {
    from <- graph$edges$from[[i]]
    to <- graph$edges$to[[i]]
    order <- graph$edges$order[[i]]
    adjacent_atoms[[from]] <- c(adjacent_atoms[[from]], to)
    adjacent_orders[[from]] <- c(adjacent_orders[[from]], order)
    adjacent_atoms[[to]] <- c(adjacent_atoms[[to]], from)
    adjacent_orders[[to]] <- c(adjacent_orders[[to]], order)
  }

  in_alternating_six_ring <- function(start) {
    walk <- function(current, path, orders) {
      if (length(path) == 6L) {
        closing_index <- match(start, adjacent_atoms[[current]])
        if (is.na(closing_index)) return(FALSE)
        complete_orders <- c(
          orders,
          adjacent_orders[[current]][[closing_index]]
        )
        shifted_orders <- c(complete_orders[-1L], complete_orders[[1L]])
        return(
          all(complete_orders %in% c(1, 2)) &&
            sum(complete_orders == 2) == 3L &&
            all(complete_orders != shifted_orders)
        )
      }

      candidates <- adjacent_atoms[[current]]
      candidate_orders <- adjacent_orders[[current]]
      keep <- allowed[candidates] & !candidates %in% path &
        candidate_orders %in% c(1, 2)
      if (length(orders)) {
        keep <- keep & candidate_orders != orders[[length(orders)]]
      }
      candidates <- candidates[keep]
      candidate_orders <- candidate_orders[keep]
      if (!length(candidates)) return(FALSE)

      any(vapply(seq_along(candidates), function(i) {
        walk(
          candidates[[i]],
          c(path, candidates[[i]]),
          c(orders, candidate_orders[[i]])
        )
      }, logical(1)))
    }

    walk(start, start, numeric())
  }

  candidates <- which(allowed)
  perceived[candidates] <- vapply(
    candidates,
    in_alternating_six_ring,
    logical(1)
  )
  perceived
}

.smiles_hydrogen_count <- function(graph, atom) {
  atom_info <- graph$atoms[atom, , drop = FALSE]
  neighbours <- .smiles_neighbours(graph, atom)
  attached_h <- sum(graph$atoms$element[neighbours$atom] == "H")
  explicit <- atom_info$explicit_h[[1]] + attached_h
  if (atom_info$bracket[[1]]) return(explicit)

  heavy <- graph$atoms$element[neighbours$atom] != "H"
  heavy_orders <- neighbours$order[heavy]
  bond_sum <- sum(heavy_orders)
  degree <- sum(heavy)
  element <- atom_info$element[[1]]
  implicit <- 0
  if (element == "C") {
    implicit <- if (atom_info$aromatic[[1]]) {
      if (degree < 3L) 1 else 0
    } else max(0, round(4 - bond_sum))
  } else if (element == "N") {
    implicit <- if (atom_info$aromatic[[1]]) 0 else {
      target <- if (atom_info$charge[[1]] > 0) 4 else 3
      max(0, round(target - bond_sum))
    }
  } else if (element %in% c("O", "S")) {
    implicit <- if (atom_info$charge[[1]] < 0) 0 else {
      max(0, round(2 - bond_sum))
    }
  }
  explicit + implicit
}

.deconvolve_smiles_one <- function(smiles) {
  graph <- .parse_smiles_graph(smiles)
  empty <- list(
    structure_valid = FALSE,
    structure_message = graph$message,
    carboxyl_sites = NA_real_, hydroxyl_sites = NA_real_,
    alcohol_hydroxyl_sites = NA_real_, phenolic_hydroxyl_sites = NA_real_,
    ketone_sites = NA_real_, aldehyde_sites = NA_real_,
    carbonyl_sites = NA_real_, primary_amine_sites = NA_real_,
    secondary_amine_sites = NA_real_, tertiary_amine_sites = NA_real_,
    aromatic_nitrogen_sites = NA_real_, amide_sites = NA_real_,
    phosphate_acid_sites = NA_real_, sulfonic_acid_sites = NA_real_,
    thiol_sites = NA_real_, catechol_sites = NA_real_,
    acidic_sites = NA_real_, weak_acidic_sites = NA_real_,
    basic_sites = NA_real_, proton_acceptor_sites = NA_real_,
    alkali_binding_sites = NA_real_, alkali_chelation_motifs = NA_real_,
    max_exchangeable_protons = NA_real_, fmp10_reactive_sites = NA_real_,
    positive_mode_score = NA_real_, negative_mode_score = NA_real_,
    alkali_affinity_score = NA_real_, neutral_mass_score = NA_real_,
    positive_site_atoms = NA_character_, negative_site_atoms = NA_character_,
    alkali_site_atoms = NA_character_,
    structure_evidence = NA_character_
  )
  if (!isTRUE(graph$valid)) return(empty)

  atoms <- graph$atoms
  atom_ids <- seq_len(nrow(atoms))
  perceived_aromatic <- atoms$aromatic |
    .smiles_kekule_aromatic_atoms(graph)
  hydrogen_count <- vapply(
    atom_ids, function(i) .smiles_hydrogen_count(graph, i), numeric(1)
  )
  neighbours <- lapply(atom_ids, function(i) .smiles_neighbours(graph, i))
  heavy_neighbours <- lapply(neighbours, function(x) {
    keep <- atoms$element[x$atom] != "H"
    list(atom = x$atom[keep], order = x$order[keep])
  })

  carbonyl_centres <- atom_ids[vapply(atom_ids, function(i) {
    if (atoms$element[[i]] != "C") return(FALSE)
    n <- heavy_neighbours[[i]]
    any(atoms$element[n$atom] == "O" & n$order >= 1.9)
  }, logical(1))]

  carboxyl_centres <- carbonyl_centres[vapply(carbonyl_centres, function(i) {
    n <- heavy_neighbours[[i]]
    oxygen <- n$atom[atoms$element[n$atom] == "O" & n$order < 1.5]
    any(vapply(oxygen, function(o) {
      on <- heavy_neighbours[[o]]
      length(on$atom) == 1L &&
        (hydrogen_count[[o]] > 0 || atoms$charge[[o]] < 0)
    }, logical(1)))
  }, logical(1))]

  amide_centres <- carbonyl_centres[vapply(carbonyl_centres, function(i) {
    n <- heavy_neighbours[[i]]
    any(atoms$element[n$atom] == "N" & n$order < 1.5)
  }, logical(1))]

  ketone_centres <- setdiff(carbonyl_centres, c(carboxyl_centres, amide_centres))
  ketone_centres <- ketone_centres[vapply(ketone_centres, function(i) {
    n <- heavy_neighbours[[i]]
    substituents <- n$atom[!(atoms$element[n$atom] == "O" & n$order >= 1.9)]
    sum(atoms$element[substituents] == "C") == 2L
  }, logical(1))]
  aldehyde_centres <- setdiff(carbonyl_centres, c(carboxyl_centres, amide_centres))
  aldehyde_centres <- aldehyde_centres[vapply(aldehyde_centres, function(i) {
    n <- heavy_neighbours[[i]]
    substituents <- n$atom[!(atoms$element[n$atom] == "O" & n$order >= 1.9)]
    sum(atoms$element[substituents] == "C") == 1L && hydrogen_count[[i]] > 0
  }, logical(1))]

  hydroxyl_oxygen <- atom_ids[vapply(atom_ids, function(i) {
    if (atoms$element[[i]] != "O" || atoms$charge[[i]] != 0) return(FALSE)
    n <- heavy_neighbours[[i]]
    length(n$atom) == 1L && n$order[[1]] < 1.5 && hydrogen_count[[i]] > 0
  }, logical(1))]
  carboxyl_oxygen <- unique(unlist(lapply(carboxyl_centres, function(i) {
    n <- heavy_neighbours[[i]]
    n$atom[atoms$element[n$atom] == "O" & n$order < 1.5]
  }), use.names = FALSE))
  phenolic_oxygen <- hydroxyl_oxygen[vapply(hydroxyl_oxygen, function(i) {
    attached <- heavy_neighbours[[i]]$atom
    length(attached) == 1L && perceived_aromatic[[attached]] &&
      atoms$element[[attached]] == "C"
  }, logical(1))]
  alcohol_oxygen <- setdiff(hydroxyl_oxygen, c(carboxyl_oxygen, phenolic_oxygen))

  amide_nitrogen <- unique(unlist(lapply(amide_centres, function(i) {
    n <- heavy_neighbours[[i]]
    n$atom[atoms$element[n$atom] == "N" & n$order < 1.5]
  }), use.names = FALSE))
  amine_nitrogen <- atom_ids[
    atoms$element == "N" & !perceived_aromatic & atoms$charge >= 0 &
      !atom_ids %in% amide_nitrogen
  ]
  amine_degree <- vapply(amine_nitrogen, function(i) {
    sum(atoms$element[heavy_neighbours[[i]]$atom] != "H")
  }, integer(1))
  primary_amine <- sum(amine_degree == 1L & hydrogen_count[amine_nitrogen] >= 1)
  secondary_amine <- sum(amine_degree == 2L & hydrogen_count[amine_nitrogen] >= 1)
  tertiary_amine <- sum(amine_degree >= 3L)
  aromatic_n <- sum(
    atoms$element == "N" & perceived_aromatic & atoms$charge <= 0 &
      hydrogen_count == 0
  )

  acidic_hetero_sites <- function(element, centre, minimum_double_o) {
    centres <- atom_ids[atoms$element == centre]
    sum(vapply(centres, function(i) {
      n <- heavy_neighbours[[i]]
      double_o <- sum(atoms$element[n$atom] == "O" & n$order >= 1.9)
      single_o <- n$atom[atoms$element[n$atom] == "O" & n$order < 1.5]
      if (double_o < minimum_double_o) return(0L)
      sum(hydrogen_count[single_o] > 0)
    }, integer(1)))
  }
  phosphate_acid <- acidic_hetero_sites("O", "P", 1L)
  sulfonic_acid <- acidic_hetero_sites("O", "S", 2L)
  thiol <- sum(atoms$element == "S" & hydrogen_count > 0)

  phenolic_carbon <- unique(vapply(phenolic_oxygen, function(i) {
    heavy_neighbours[[i]]$atom[[1]]
  }, integer(1)))
  catechol <- 0L
  if (length(phenolic_carbon) >= 2L) {
    pairs <- utils::combn(phenolic_carbon, 2L)
    catechol <- sum(apply(pairs, 2L, function(pair) {
      any((graph$edges$from == pair[[1]] & graph$edges$to == pair[[2]]) |
            (graph$edges$from == pair[[2]] & graph$edges$to == pair[[1]]))
    }))
  }

  neutral_n <- atom_ids[atoms$element == "N" & atoms$charge <= 0]
  donor_atoms <- unique(c(
    atom_ids[atoms$element == "O" & atoms$charge <= 0],
    setdiff(neutral_n, amide_nitrogen),
    atom_ids[atoms$element == "S" & atoms$charge <= 0]
  ))
  donor_pairs <- 0L
  if (length(donor_atoms) >= 2L) {
    pairs <- utils::combn(donor_atoms, 2L)
    donor_pairs <- sum(apply(pairs, 2L, function(pair) {
      n1 <- heavy_neighbours[[pair[[1]]]]$atom
      n2 <- heavy_neighbours[[pair[[2]]]]$atom
      length(intersect(n1, n2)) > 0L || pair[[2]] %in% n1
    }))
  }
  chelation_motifs <- max(
    length(carboxyl_centres) + catechol + phosphate_acid + sulfonic_acid,
    donor_pairs
  )

  strong_acid <- length(carboxyl_centres) + phosphate_acid + sulfonic_acid
  medium_acid <- length(phenolic_oxygen) + thiol
  weak_acid <- length(alcohol_oxygen)
  basic_sites <- primary_amine + secondary_amine + tertiary_amine + aromatic_n
  carbonyl_oxygen <- unique(unlist(lapply(carbonyl_centres, function(i) {
    n <- heavy_neighbours[[i]]
    n$atom[atoms$element[n$atom] == "O" & n$order >= 1.9]
  }), use.names = FALSE))
  proton_acceptors <- length(unique(c(
    carbonyl_oxygen,
    atom_ids[atoms$element == "O" & hydrogen_count == 0 & atoms$charge <= 0],
    setdiff(neutral_n, amide_nitrogen)
  )))

  positive_score <- if (basic_sites > 0L) {
    min(1, 0.78 + 0.07 * basic_sites)
  } else if (proton_acceptors > 0L) {
    min(0.72, 0.42 + 0.05 * proton_acceptors)
  } else if (any(perceived_aromatic)) 0.3 else 0.2
  negative_score <- if (strong_acid > 0L) {
    min(1, 0.82 + 0.06 * strong_acid)
  } else if (medium_acid > 0L) {
    min(0.78, 0.58 + 0.05 * medium_acid)
  } else if (weak_acid > 0L) {
    min(0.5, 0.3 + 0.03 * weak_acid)
  } else 0.15
  alkali_score <- if (!length(donor_atoms)) 0.15 else {
    min(1, 0.35 + 0.08 * min(length(donor_atoms), 5L) +
          0.12 * min(chelation_motifs, 2L))
  }
  neutral_score <- if (sum(atoms$charge) == 0L) 1 else 0.4
  fmp10_sites <- primary_amine + secondary_amine + length(phenolic_oxygen)

  atom_labels <- function(ids) {
    ids <- sort(unique(as.integer(ids)))
    if (!length(ids)) return("none")
    paste0(atoms$element[ids], ids, collapse = ",")
  }
  positive_atoms <- unique(c(
    setdiff(neutral_n, amide_nitrogen), carbonyl_oxygen,
    atom_ids[atoms$element == "O" & hydrogen_count == 0 & atoms$charge <= 0]
  ))
  negative_atoms <- unique(c(
    carboxyl_oxygen, phenolic_oxygen,
    atom_ids[atoms$element == "S" & hydrogen_count > 0]
  ))
  positive_site_atoms <- atom_labels(positive_atoms)
  negative_site_atoms <- atom_labels(negative_atoms)
  alkali_site_atoms <- atom_labels(donor_atoms)

  evidence <- paste0(
    "carboxyl=", length(carboxyl_centres),
    "; alcohol-OH=", length(alcohol_oxygen),
    "; phenolic-OH=", length(phenolic_oxygen),
    "; ketone=", length(ketone_centres),
    "; amine=", primary_amine + secondary_amine + tertiary_amine,
    "; acidic-sites=", strong_acid + medium_acid,
    "; proton-acceptors=", proton_acceptors,
    "; alkali-donors=", length(donor_atoms),
    "; chelation-motifs=", chelation_motifs,
    "; positive-atoms=", positive_site_atoms,
    "; negative-atoms=", negative_site_atoms,
    "; alkali-atoms=", alkali_site_atoms
  )

  list(
    structure_valid = TRUE,
    structure_message = "parsed",
    carboxyl_sites = length(carboxyl_centres),
    hydroxyl_sites = length(hydroxyl_oxygen),
    alcohol_hydroxyl_sites = length(alcohol_oxygen),
    phenolic_hydroxyl_sites = length(phenolic_oxygen),
    ketone_sites = length(ketone_centres),
    aldehyde_sites = length(aldehyde_centres),
    carbonyl_sites = length(carbonyl_centres),
    primary_amine_sites = primary_amine,
    secondary_amine_sites = secondary_amine,
    tertiary_amine_sites = tertiary_amine,
    aromatic_nitrogen_sites = aromatic_n,
    amide_sites = length(amide_centres),
    phosphate_acid_sites = phosphate_acid,
    sulfonic_acid_sites = sulfonic_acid,
    thiol_sites = thiol,
    catechol_sites = catechol,
    acidic_sites = strong_acid + medium_acid,
    weak_acidic_sites = weak_acid,
    basic_sites = basic_sites,
    proton_acceptor_sites = proton_acceptors,
    alkali_binding_sites = length(donor_atoms),
    alkali_chelation_motifs = chelation_motifs,
    max_exchangeable_protons = strong_acid + medium_acid,
    fmp10_reactive_sites = fmp10_sites,
    positive_mode_score = positive_score,
    negative_mode_score = negative_score,
    alkali_affinity_score = alkali_score,
    neutral_mass_score = neutral_score,
    positive_site_atoms = positive_site_atoms,
    negative_site_atoms = negative_site_atoms,
    alkali_site_atoms = alkali_site_atoms,
    structure_evidence = evidence
  )
}

#' Decompose explicit functional groups in SMILES strings
#'
#' A lightweight, dependency-free SMILES graph parser identifies common
#' metabolite functional groups and derives interpretable positive-ion,
#' negative-ion, alkali-binding, and neutral-mass priors. The parser is intended
#' for candidate prioritisation rather than pKa, proton-affinity, or quantum
#' chemistry calculations. Stereochemistry is accepted but does not change the
#' counts.
#'
#' @param smiles Character vector of SMILES strings.
#' @param backend Currently `"auto"` and `"native"` both use SpaMTP's native
#'   parser. The argument reserves a stable extension point for optional
#'   chemistry toolkits.
#' @param strict Stop when any SMILES cannot be parsed. When `FALSE`, invalid
#'   rows are returned with `structure_valid = FALSE` and missing features.
#' @param workers Number of forked workers used for unique SMILES on platforms
#'   supporting `parallel::mclapply()`. The default is controlled by
#'   `options(SpaMTP.smiles_workers = 1)`.
#'
#' @return A data frame with functional-group counts, ion-mode scores, and a
#'   human-readable `structure_evidence` field.
#' @examples
#' utils::str(formals(DeconvolveSMILES))
#' @export
DeconvolveSMILES <- function(smiles, backend = c("auto", "native"),
                             strict = FALSE,
                             workers = getOption("SpaMTP.smiles_workers", 1L)) {
  backend <- match.arg(backend)
  workers <- suppressWarnings(as.integer(workers))
  if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    stop("workers must be a positive integer.")
  }
  smiles <- as.character(smiles)
  if (!length(smiles)) {
    return(data.frame(
      smiles = character(), structure_backend = character(),
      structure_valid = logical(), stringsAsFactors = FALSE
    ))
  }
  unique_smiles <- unique(smiles)
  parsed <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(
      unique_smiles, .deconvolve_smiles_one, mc.cores = workers,
      mc.preschedule = TRUE
    )
  } else {
    lapply(unique_smiles, .deconvolve_smiles_one)
  }
  valid <- vapply(parsed, `[[`, logical(1), "structure_valid")
  if (isTRUE(strict) && any(!valid)) {
    bad <- unique_smiles[!valid]
    messages <- vapply(parsed[!valid], `[[`, character(1), "structure_message")
    stop(
      "Could not parse SMILES: ",
      paste(paste0(bad, " (", messages, ")"), collapse = "; ")
    )
  }
  rows <- lapply(seq_along(parsed), function(i) {
    data.frame(
      smiles = unique_smiles[[i]],
      structure_backend = "SpaMTP-native",
      as.data.frame(parsed[[i]], stringsAsFactors = FALSE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result <- result[match(smiles, result$smiles), , drop = FALSE]
  rownames(result) <- NULL
  result
}

#' Add SMILES-derived structural features to a metabolite database
#'
#' @param db A metabolite data frame.
#' @param smiles_column Column containing SMILES. When `NULL`, SpaMTP detects
#'   `iso_smiles`, `canonical_smiles`, or `smiles`.
#' @param backend Parser backend passed to [DeconvolveSMILES()].
#' @param overwrite Replace existing feature values. By default only absent or
#'   missing values are filled.
#' @param strict Passed to [DeconvolveSMILES()].
#' @param workers Parallel workers passed to [DeconvolveSMILES()].
#'
#' @return `db` with structure-derived columns appended or completed.
#' @examples
#' utils::str(formals(AnnotateSMILESStructure))
#' @export
AnnotateSMILESStructure <- function(db, smiles_column = NULL,
                                    backend = c("auto", "native"),
                                    overwrite = FALSE, strict = FALSE,
                                    workers = getOption("SpaMTP.smiles_workers", 1L)) {
  if (!is.data.frame(db)) stop("db must be a data.frame or data.table.")
  existing_bound <- if ("max_exchangeable_protons" %in% names(db)) {
    is.finite(suppressWarnings(as.numeric(db$max_exchangeable_protons)))
  } else rep(FALSE, nrow(db))
  existing_bound_source <- if ("proton_bound_source" %in% names(db)) {
    as.character(db$proton_bound_source)
  } else rep(NA_character_, nrow(db))
  backend <- match.arg(backend)
  if (is.null(smiles_column)) {
    choices <- c("iso_smiles", "canonical_smiles", "smiles", "SMILES")
    detected <- choices[choices %in% names(db)]
    smiles_column <- if (length(detected)) detected[[1]] else NULL
  }
  if (is.null(smiles_column) || !smiles_column %in% names(db)) {
    stop("No SMILES column was found in db.")
  }
  features <- DeconvolveSMILES(db[[smiles_column]], backend = backend,
                              strict = strict, workers = workers)
  feature_names <- setdiff(names(features), "smiles")
  for (column in feature_names) {
    if (!column %in% names(db) || isTRUE(overwrite)) {
      db[[column]] <- features[[column]]
    } else {
      missing <- is.na(db[[column]])
      if (is.character(db[[column]])) missing <- missing | !nzchar(db[[column]])
      db[[column]][missing] <- features[[column]][missing]
    }
  }
  db$structure_smiles <- as.character(db[[smiles_column]])
  db$proton_bound_source <- ifelse(
    existing_bound & (is.na(existing_bound_source) |
                        existing_bound_source != "SMILES-inferred"),
    "provided",
    ifelse(features$structure_valid, "SMILES-inferred", "unavailable")
  )
  db
}
