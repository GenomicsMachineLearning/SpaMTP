#### Indexed and chemically constrained m/z annotation ########################

.spamtp_mass_constants <- function() {
  electron <- 0.000548579909065
  hydrogen <- 1.00782503223
  oxygen <- 15.99491461957
  carbon <- 12
  nitrogen <- 14.00307400443
  sulfur <- 31.9720711744
  fluorine <- 18.99840316273

  list(
    electron = electron,
    proton = 1.007276466621,
    hydrogen = hydrogen,
    sodium_ion = 22.9897692820 - electron,
    potassium_ion = 38.9637064864 - electron,
    chloride = 34.968852682 + electron,
    bromide = 78.9183376 + electron,
    ammonium = nitrogen + 4 * hydrogen - electron,
    water = 2 * hydrogen + oxygen,
    acetonitrile = 2 * carbon + 3 * hydrogen + nitrogen,
    methanol = carbon + 4 * hydrogen + oxygen,
    isopropanol = 3 * carbon + 8 * hydrogen + oxygen,
    dmso = 2 * carbon + 6 * hydrogen + oxygen + sulfur,
    formate = carbon + hydrogen + 2 * oxygen + electron,
    acetate = 2 * carbon + 3 * hydrogen + 2 * oxygen + electron,
    trifluoroacetate = 2 * carbon + 3 * fluorine + 2 * oxygen + electron,
    carbon13_spacing = 1.00335483507,
    chlorine_m2_spacing = 1.99704991,
    bromine_m2_spacing = 1.9979535
  )
}

.make_adduct_rule <- function(name, polarity, n_molecules, charge,
                              mass_shift, loss_h = 0L,
                              added_charge = charge + loss_h,
                              complexity = "simple", base_adduct = name,
                              contains_halogen = FALSE,
                              contains_metal = FALSE, prior = 1) {
  data.frame(
    name = name,
    notation = paste0(
      "[", name, "]",
      ifelse(abs(charge) == 1, "", abs(charge)),
      ifelse(charge > 0, "+", "-")
    ),
    polarity = polarity,
    n_molecules = as.integer(n_molecules),
    charge = as.integer(charge),
    mass_shift = as.numeric(mass_shift),
    loss_h = as.integer(loss_h),
    added_charge = as.integer(added_charge),
    complexity = complexity,
    base_adduct = base_adduct,
    contains_halogen = contains_halogen,
    contains_metal = contains_metal,
    prior = as.numeric(prior),
    stringsAsFactors = FALSE
  )
}

#' Return the validated SpaMTP adduct rule table
#'
#' The rules use neutral monoisotopic mass and explicit ion masses. Each row
#' stores the molecular multiplier, total numerator mass shift, charge, number
#' of protons removed from the analyte, and the charge contributed by added
#' species. This allows rules to be checked before candidate generation.
#'
#' @param polarity One of `"both"`, `"positive"`, `"negative"`, or
#'   `"neutral"`.
#' @param include_complex Include multimers, solvent clusters, metal-exchange,
#'   and multiply charged ions.
#'
#' @return A data frame containing validated adduct rules.
#' @export
AdductRules <- function(polarity = c("both", "positive", "negative", "neutral"),
                        include_complex = TRUE) {
  polarity <- match.arg(polarity)
  m <- .spamtp_mass_constants()

  rules <- do.call(rbind, list(
    .make_adduct_rule("M+H", "positive", 1, 1, m$proton),
    .make_adduct_rule("M+NH4", "positive", 1, 1, m$ammonium,
                      complexity = "salt", base_adduct = "M+H", prior = 0.9),
    .make_adduct_rule("M+Na", "positive", 1, 1, m$sodium_ion,
                      complexity = "salt", contains_metal = TRUE, prior = 0.9),
    .make_adduct_rule("M+K", "positive", 1, 1, m$potassium_ion,
                      complexity = "salt", contains_metal = TRUE, prior = 0.8),
    .make_adduct_rule("M+2H", "positive", 1, 2, 2 * m$proton,
                      complexity = "multicharged", base_adduct = "M+H", prior = 0.8),
    .make_adduct_rule("M+H+Na", "positive", 1, 2,
                      m$proton + m$sodium_ion,
                      complexity = "multicharged", base_adduct = "M+Na",
                      contains_metal = TRUE, prior = 0.7),
    .make_adduct_rule("M+2Na", "positive", 1, 2, 2 * m$sodium_ion,
                      complexity = "multicharged", base_adduct = "M+Na",
                      contains_metal = TRUE, prior = 0.65),
    .make_adduct_rule("M+H+K", "positive", 1, 2,
                      m$proton + m$potassium_ion,
                      complexity = "multicharged", base_adduct = "M+K",
                      contains_metal = TRUE, prior = 0.65),
    .make_adduct_rule("M+3H", "positive", 1, 3, 3 * m$proton,
                      complexity = "multicharged", base_adduct = "M+H", prior = 0.65),
    .make_adduct_rule("M+2H+Na", "positive", 1, 3,
                      2 * m$proton + m$sodium_ion,
                      complexity = "multicharged", base_adduct = "M+Na",
                      contains_metal = TRUE, prior = 0.55),
    .make_adduct_rule("M+H+2Na", "positive", 1, 3,
                      m$proton + 2 * m$sodium_ion,
                      complexity = "multicharged", base_adduct = "M+Na",
                      contains_metal = TRUE, prior = 0.5),
    .make_adduct_rule("M+3Na", "positive", 1, 3, 3 * m$sodium_ion,
                      complexity = "multicharged", base_adduct = "M+Na",
                      contains_metal = TRUE, prior = 0.45),
    .make_adduct_rule("M+2Na-H", "positive", 1, 1,
                      2 * m$sodium_ion - m$proton, loss_h = 1,
                      added_charge = 2, complexity = "metal_exchange",
                      base_adduct = "M+Na", contains_metal = TRUE, prior = 0.65),
    .make_adduct_rule("M+2K-H", "positive", 1, 1,
                      2 * m$potassium_ion - m$proton, loss_h = 1,
                      added_charge = 2, complexity = "metal_exchange",
                      base_adduct = "M+K", contains_metal = TRUE, prior = 0.5),
    .make_adduct_rule("M+ACN+H", "positive", 1, 1,
                      m$acetonitrile + m$proton, complexity = "solvent",
                      base_adduct = "M+H", prior = 0.5),
    .make_adduct_rule("M+ACN+Na", "positive", 1, 1,
                      m$acetonitrile + m$sodium_ion, complexity = "solvent",
                      base_adduct = "M+Na", contains_metal = TRUE, prior = 0.45),
    .make_adduct_rule("M+2ACN+H", "positive", 1, 1,
                      2 * m$acetonitrile + m$proton, complexity = "solvent",
                      base_adduct = "M+H", prior = 0.35),
    .make_adduct_rule("M+CH3OH+H", "positive", 1, 1,
                      m$methanol + m$proton, complexity = "solvent",
                      base_adduct = "M+H", prior = 0.45),
    .make_adduct_rule("M+IsoProp+H", "positive", 1, 1,
                      m$isopropanol + m$proton, complexity = "solvent",
                      base_adduct = "M+H", prior = 0.35),
    .make_adduct_rule("M+DMSO+H", "positive", 1, 1,
                      m$dmso + m$proton, complexity = "solvent",
                      base_adduct = "M+H", prior = 0.35),
    .make_adduct_rule("2M+H", "positive", 2, 1, m$proton,
                      complexity = "multimer", base_adduct = "M+H", prior = 0.55),
    .make_adduct_rule("2M+NH4", "positive", 2, 1, m$ammonium,
                      complexity = "multimer", base_adduct = "M+H", prior = 0.45),
    .make_adduct_rule("2M+Na", "positive", 2, 1, m$sodium_ion,
                      complexity = "multimer", base_adduct = "M+Na",
                      contains_metal = TRUE, prior = 0.5),
    .make_adduct_rule("2M+K", "positive", 2, 1, m$potassium_ion,
                      complexity = "multimer", base_adduct = "M+K",
                      contains_metal = TRUE, prior = 0.4),
    .make_adduct_rule("2M+ACN+H", "positive", 2, 1,
                      m$acetonitrile + m$proton, complexity = "multimer",
                      base_adduct = "M+H", prior = 0.3),
    .make_adduct_rule("2M+ACN+Na", "positive", 2, 1,
                      m$acetonitrile + m$sodium_ion, complexity = "multimer",
                      base_adduct = "M+Na", contains_metal = TRUE, prior = 0.3),

    .make_adduct_rule("M-H", "negative", 1, -1, -m$proton,
                      loss_h = 1, added_charge = 0),
    .make_adduct_rule("M-2H", "negative", 1, -2, -2 * m$proton,
                      loss_h = 2, added_charge = 0,
                      complexity = "multicharged", base_adduct = "M-H", prior = 0.8),
    .make_adduct_rule("M-3H", "negative", 1, -3, -3 * m$proton,
                      loss_h = 3, added_charge = 0,
                      complexity = "multicharged", base_adduct = "M-H", prior = 0.6),
    .make_adduct_rule("M+Cl", "negative", 1, -1, m$chloride,
                      added_charge = -1, complexity = "salt",
                      base_adduct = "M-H", contains_halogen = TRUE, prior = 0.65),
    .make_adduct_rule("M+Br", "negative", 1, -1, m$bromide,
                      added_charge = -1, complexity = "salt",
                      base_adduct = "M-H", contains_halogen = TRUE, prior = 0.55),
    .make_adduct_rule("M+Na-2H", "negative", 1, -1,
                      m$sodium_ion - 2 * m$proton, loss_h = 2,
                      added_charge = 1, complexity = "metal_exchange",
                      base_adduct = "M-H", contains_metal = TRUE, prior = 0.6),
    .make_adduct_rule("M+K-2H", "negative", 1, -1,
                      m$potassium_ion - 2 * m$proton, loss_h = 2,
                      added_charge = 1, complexity = "metal_exchange",
                      base_adduct = "M-H", contains_metal = TRUE, prior = 0.5),
    .make_adduct_rule("M+FA-H", "negative", 1, -1, m$formate,
                      added_charge = -1, complexity = "solvent",
                      base_adduct = "M-H", prior = 0.55),
    .make_adduct_rule("M+Hac-H", "negative", 1, -1, m$acetate,
                      added_charge = -1, complexity = "solvent",
                      base_adduct = "M-H", prior = 0.5),
    .make_adduct_rule("M+TFA-H", "negative", 1, -1,
                      m$trifluoroacetate, added_charge = -1,
                      complexity = "solvent", base_adduct = "M-H",
                      contains_halogen = TRUE, prior = 0.35),
    .make_adduct_rule("M-H2O-H", "negative", 1, -1,
                      -m$water - m$proton, loss_h = 1,
                      added_charge = 0, complexity = "neutral_loss",
                      base_adduct = "M-H", prior = 0.55),
    .make_adduct_rule("2M-H", "negative", 2, -1, -m$proton,
                      loss_h = 1, added_charge = 0,
                      complexity = "multimer", base_adduct = "M-H", prior = 0.55),
    .make_adduct_rule("2M+FA-H", "negative", 2, -1, m$formate,
                      added_charge = -1, complexity = "multimer",
                      base_adduct = "M-H", prior = 0.4),
    .make_adduct_rule("2M+Hac-H", "negative", 2, -1, m$acetate,
                      added_charge = -1, complexity = "multimer",
                      base_adduct = "M-H", prior = 0.35),
    .make_adduct_rule("3M-H", "negative", 3, -1, -m$proton,
                      loss_h = 1, added_charge = 0,
                      complexity = "multimer", base_adduct = "M-H", prior = 0.3),

    .make_adduct_rule("M", "neutral", 1, 1, 0,
                      added_charge = 0, complexity = "neutral",
                      base_adduct = "M", prior = 1)
  ))

  rownames(rules) <- NULL
  if (polarity != "both") {
    rules <- rules[rules$polarity == polarity, , drop = FALSE]
  }
  if (!isTRUE(include_complex)) {
    rules <- rules[rules$complexity %in% c("simple", "salt", "neutral"), , drop = FALSE]
  }
  rules
}

.normalise_adduct_name <- function(x) {
  x <- gsub("\\s+", "", as.character(x))
  x <- sub("^\\[", "", x)
  x <- sub("\\](?:[0-9]*[+-]|[+-][0-9]*)$", "", x)
  aliases <- c("M+2K+H" = "M+2K-H")
  hit <- x %in% names(aliases)
  x[hit] <- unname(aliases[x[hit]])
  x
}

.validate_adduct_rules <- function(rules) {
  required <- c(
    "name", "polarity", "n_molecules", "charge", "mass_shift",
    "loss_h", "added_charge", "complexity", "base_adduct",
    "contains_halogen", "contains_metal", "prior"
  )
  missing_columns <- setdiff(required, names(rules))
  if (length(missing_columns)) {
    stop("Adduct rule table is missing: ", paste(missing_columns, collapse = ", "))
  }
  if (any(!is.finite(rules$n_molecules) | rules$n_molecules < 1)) {
    stop("Every adduct rule must have n_molecules >= 1.")
  }
  if (any(!is.finite(rules$charge) | rules$charge == 0)) {
    stop("Every adduct rule must use a non-zero observed charge.")
  }
  non_neutral <- rules$polarity != "neutral"
  balanced <- rules$added_charge - rules$loss_h == rules$charge
  if (any(non_neutral & !balanced)) {
    bad <- rules$name[non_neutral & !balanced]
    stop("Unbalanced adduct rule(s): ", paste(bad, collapse = ", "),
         ". added_charge - loss_h must equal charge.")
  }
  if (any(!is.finite(rules$mass_shift)) || any(!is.finite(rules$prior))) {
    stop("Adduct mass shifts and priors must be finite.")
  }
  rules
}

.first_existing_column <- function(x, choices, required = TRUE) {
  hit <- choices[choices %in% names(x)]
  if (length(hit)) {
    return(hit[[1]])
  }
  if (required) {
    stop("Database must contain one of these columns: ",
         paste(choices, collapse = ", "))
  }
  NULL
}

.paste_unique <- function(x) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (length(x)) paste(x, collapse = "; ") else NA_character_
}

.normalise_metabolite_db <- function(db, collapse_isomers = TRUE) {
  if (!is.data.frame(db)) {
    stop("db must be a data.frame or data.table.")
  }
  formula_col <- .first_existing_column(db, c("formula", "mol_formula"))
  mass_col <- .first_existing_column(
    db, c("exactmass", "monoisotop_mass", "monoisotopic_mass", "neutral_mass")
  )
  id_col <- .first_existing_column(
    db, c("isomers", "chem_source_id", "compound_id", "id", "ramp_id"),
    required = FALSE
  )
  name_col <- .first_existing_column(
    db, c("isomers_names", "common_name", "name"), required = FALSE
  )
  inchikey_col <- .first_existing_column(
    db, c("isomers_inchikey", "inchi_key", "inchikey"), required = FALSE
  )
  ramp_col <- .first_existing_column(db, c("ramp_id", "rampId"), required = FALSE)
  proton_col <- .first_existing_column(
    db,
    c("max_exchangeable_protons", "max_protons", "exchangeable_protons"),
    required = FALSE
  )

  n <- nrow(db)
  value_or_na <- function(column) {
    if (is.null(column)) rep(NA_character_, n) else as.character(db[[column]])
  }
  compounds <- data.frame(
    formula = as.character(db[[formula_col]]),
    exactmass = suppressWarnings(as.numeric(db[[mass_col]])),
    isomers = value_or_na(id_col),
    isomers_inchikey = value_or_na(inchikey_col),
    isomers_names = value_or_na(name_col),
    ramp_ids = value_or_na(ramp_col),
    max_exchangeable_protons = if (is.null(proton_col)) {
      rep(NA_real_, n)
    } else {
      suppressWarnings(as.numeric(db[[proton_col]]))
    },
    stringsAsFactors = FALSE
  )
  keep <- is.finite(compounds$exactmass) & compounds$exactmass > 0 &
    !is.na(compounds$formula) & nzchar(compounds$formula)
  compounds <- compounds[keep, , drop = FALSE]

  if (isTRUE(collapse_isomers) && nrow(compounds)) {
    group_key <- paste(
      compounds$formula,
      format(compounds$exactmass, digits = 17, scientific = FALSE, trim = TRUE),
      ifelse(is.na(compounds$max_exchangeable_protons), "NA",
             compounds$max_exchangeable_protons),
      sep = "\r"
    )
    if (anyDuplicated(group_key)) {
      grouped_rows <- split(seq_len(nrow(compounds)), group_key)
      first_row <- vapply(grouped_rows, function(i) i[[1]], integer(1))
      collapse_column <- function(column) {
        vapply(grouped_rows, function(i) .paste_unique(column[i]), character(1))
      }
      compounds <- data.frame(
        formula = compounds$formula[first_row],
        exactmass = compounds$exactmass[first_row],
        isomers = collapse_column(compounds$isomers),
        isomers_inchikey = collapse_column(compounds$isomers_inchikey),
        isomers_names = collapse_column(compounds$isomers_names),
        ramp_ids = collapse_column(compounds$ramp_ids),
        max_exchangeable_protons =
          compounds$max_exchangeable_protons[first_row],
        stringsAsFactors = FALSE
      )
    }
  }
  rownames(compounds) <- NULL
  compounds
}

#' Build a reusable indexed metabolite annotation search space
#'
#' Expected adduct m/z values are generated once, sorted, and queried by binary
#' interval search. This is the one-dimensional equivalent of an interval tree
#' and avoids scanning every database row for every observed peak.
#'
#' @param db Metabolite database. Both the legacy SpaMTP five-column/wide
#'   database and the RaMP `chem_props` schema are supported.
#' @param polarity `"positive"`, `"negative"`, or `"neutral"`.
#' @param adducts Optional character vector of adduct names or bracketed
#'   notations.
#' @param rules Optional custom rule table. See [AdductRules()].
#' @param collapse_isomers Collapse records sharing formula, exact mass, and
#'   proton bound before indexing.
#'
#' @return An object of class `spamtp_mz_index`.
#' @export
BuildMZAnnotationIndex <- function(db, polarity = c("positive", "negative", "neutral"),
                                   adducts = NULL, rules = NULL,
                                   collapse_isomers = TRUE) {
  polarity <- match.arg(polarity)
  if (is.null(rules)) {
    rules <- AdductRules(polarity = polarity)
  } else {
    rules <- as.data.frame(rules, stringsAsFactors = FALSE)
    rules <- rules[rules$polarity == polarity, , drop = FALSE]
  }
  rules <- .validate_adduct_rules(rules)
  family_rules <- rules

  if (!is.null(adducts)) {
    requested <- .normalise_adduct_name(adducts)
    available <- .normalise_adduct_name(rules$name)
    missing_rules <- setdiff(requested, available)
    if (length(missing_rules)) {
      stop("Unknown or unavailable adduct(s) for ", polarity, " mode: ",
           paste(missing_rules, collapse = ", "))
    }
    rules <- rules[match(requested, available), , drop = FALSE]
  }
  if (!nrow(rules)) {
    stop("No adduct rules remain after filtering.")
  }

  compounds <- .normalise_metabolite_db(db, collapse_isomers = collapse_isomers)
  if (!nrow(compounds)) {
    stop("No valid formula/monoisotopic-mass records were found in db.")
  }

  expected_parts <- vector("list", nrow(rules))
  compound_parts <- vector("list", nrow(rules))
  rule_parts <- vector("list", nrow(rules))
  max_h <- compounds$max_exchangeable_protons

  for (j in seq_len(nrow(rules))) {
    rule <- rules[j, , drop = FALSE]
    expected <- (rule$n_molecules * compounds$exactmass + rule$mass_shift) /
      abs(rule$charge)
    valid <- is.finite(expected) & expected > 0
    have_bound <- !is.na(max_h)
    valid[have_bound] <- valid[have_bound] &
      rule$loss_h <= rule$n_molecules * max_h[have_bound]
    idx <- which(valid)
    expected_parts[[j]] <- expected[idx]
    compound_parts[[j]] <- idx
    rule_parts[[j]] <- rep.int(j, length(idx))
  }

  expected_mz <- unlist(expected_parts, use.names = FALSE)
  compound_idx <- as.integer(unlist(compound_parts, use.names = FALSE))
  rule_idx <- as.integer(unlist(rule_parts, use.names = FALSE))
  ord <- order(expected_mz, method = "radix")

  structure(
    list(
      expected_mz = expected_mz[ord],
      compound_idx = compound_idx[ord],
      rule_idx = rule_idx[ord],
      compounds = compounds,
      rules = rules,
      family_rules = family_rules,
      polarity = polarity
    ),
    class = "spamtp_mz_index"
  )
}

#' Print an indexed metabolite annotation search space
#'
#' @param x A `spamtp_mz_index` returned by [BuildMZAnnotationIndex()].
#' @param ... Additional arguments reserved for print-method compatibility.
#'
#' @return `x`, invisibly.
#' @method print spamtp_mz_index
#' @export
print.spamtp_mz_index <- function(x, ...) {
  cat("SpaMTP m/z annotation index\n")
  cat("  polarity: ", x$polarity, "\n", sep = "")
  cat("  metabolites: ", nrow(x$compounds), "\n", sep = "")
  cat("  adduct rules: ", nrow(x$rules), "\n", sep = "")
  cat("  expected ions: ", length(x$expected_mz), "\n", sep = "")
  invisible(x)
}

.normalise_ms1_spectrum <- function(ms1_spectrum) {
  if (is.null(ms1_spectrum)) {
    return(NULL)
  }
  if (is.numeric(ms1_spectrum) && is.null(dim(ms1_spectrum))) {
    return(data.frame(mz = as.numeric(ms1_spectrum), intensity = 1))
  }
  spectrum <- as.data.frame(ms1_spectrum)
  mz_col <- .first_existing_column(
    spectrum, c("mz", "m.z", "mass", "raw_mz", "observed_mz")
  )
  intensity_col <- .first_existing_column(
    spectrum, c("intensity", "abundance", "counts", "signal"), required = FALSE
  )
  intensity <- if (is.null(intensity_col)) rep(1, nrow(spectrum)) else {
    suppressWarnings(as.numeric(spectrum[[intensity_col]]))
  }
  out <- data.frame(
    mz = suppressWarnings(as.numeric(spectrum[[mz_col]])),
    intensity = intensity
  )
  out <- out[is.finite(out$mz) & out$mz > 0, , drop = FALSE]
  out[order(out$mz), , drop = FALSE]
}

.find_spectrum_peak <- function(spectrum, target_mz, ppm) {
  if (is.null(spectrum) || !nrow(spectrum) || !is.finite(target_mz)) {
    return(c(mz = NA_real_, intensity = NA_real_, ppm = NA_real_))
  }
  pos <- findInterval(target_mz, spectrum$mz)
  idx <- unique(pmax(1L, pmin(nrow(spectrum), c(pos, pos + 1L))))
  errors <- abs(spectrum$mz[idx] - target_mz) / target_mz * 1e6
  best <- idx[which.min(errors)]
  if (!length(best) || errors[which.min(errors)] > ppm) {
    return(c(mz = NA_real_, intensity = NA_real_, ppm = NA_real_))
  }
  c(mz = spectrum$mz[best], intensity = spectrum$intensity[best],
    ppm = errors[which.min(errors)])
}

.formula_elements <- function(formula) {
  unique(unlist(regmatches(formula, gregexpr("[A-Z][a-z]?", formula))))
}

.element_count <- function(formula, element) {
  pattern <- paste0("(?:^|[^a-z])", element, "([0-9]*)")
  hit <- regexec(pattern, formula, perl = TRUE)
  value <- regmatches(formula, hit)[[1]]
  if (!length(value)) return(0L)
  if (length(value) < 2 || !nzchar(value[[2]])) 1L else as.integer(value[[2]])
}

.mass_defect_score <- function(observed_mz, formula, rule) {
  signed_defect <- observed_mz - round(observed_mz)
  elements <- .formula_elements(formula)
  pure_cho <- length(elements) > 0 && all(elements %in% c("C", "H", "O"))
  suspicious <- pure_cho && signed_defect >= -0.15 && signed_defect <= -0.01 &&
    !isTRUE(rule$contains_halogen) && !isTRUE(rule$contains_metal)
  if (suspicious) 0.2 else 1
}

.ratio_score <- function(observed, expected) {
  if (!is.finite(observed) || observed <= 0 || !is.finite(expected) || expected <= 0) {
    return(0.5)
  }
  max(0.05, exp(-abs(log(observed / expected))))
}

.isotope_score <- function(expected_mz, formula, rule, spectrum, ppm) {
  if (is.null(spectrum)) return(NA_real_)
  m <- .spamtp_mass_constants()
  charge <- abs(rule$charge)
  main <- .find_spectrum_peak(spectrum, expected_mz, ppm)
  score <- 1

  carbon_count <- .element_count(formula, "C")
  if (carbon_count > 0) {
    isotope <- .find_spectrum_peak(
      spectrum, expected_mz + m$carbon13_spacing / charge, ppm
    )
    if (!is.finite(isotope[["mz"]])) {
      score <- score * 0.5
    } else if (is.finite(main[["intensity"]]) && main[["intensity"]] > 0) {
      ratio <- isotope[["intensity"]] / main[["intensity"]]
      score <- score * .ratio_score(ratio, min(0.95, carbon_count * 0.0107))
    }
  }

  n_cl <- .element_count(formula, "Cl") + as.integer(isTRUE(rule$name == "M+Cl"))
  n_br <- .element_count(formula, "Br") + as.integer(isTRUE(rule$name == "M+Br"))
  if (n_cl > 0 || n_br > 0) {
    spacing <- if (n_br > 0) m$bromine_m2_spacing else m$chlorine_m2_spacing
    expected_ratio <- if (n_br > 0) n_br * 0.97 else n_cl * 0.32
    isotope_m2 <- .find_spectrum_peak(
      spectrum, expected_mz + spacing / charge, ppm
    )
    if (!is.finite(isotope_m2[["mz"]])) {
      score <- score * 0.2
    } else if (is.finite(main[["intensity"]]) && main[["intensity"]] > 0) {
      score <- score * .ratio_score(
        isotope_m2[["intensity"]] / main[["intensity"]], expected_ratio
      )
    }
  }
  score
}

.adduct_family_score <- function(neutral_mass, rule, all_rules, spectrum, ppm) {
  complex_types <- c("multimer", "solvent", "metal_exchange")
  if (is.null(spectrum) || !rule$complexity %in% complex_types) {
    return(NA_real_)
  }
  base_names <- unique(c(
    rule$base_adduct,
    if (rule$polarity == "positive") c("M+H", "M+Na") else "M-H"
  ))
  base_rules <- all_rules[all_rules$name %in% base_names &
                            all_rules$n_molecules == 1, , drop = FALSE]
  if (!nrow(base_rules)) return(NA_real_)
  base_mz <- (base_rules$n_molecules * neutral_mass + base_rules$mass_shift) /
    abs(base_rules$charge)
  found <- vapply(base_mz, function(x) {
    is.finite(.find_spectrum_peak(spectrum, x, ppm)[["mz"]])
  }, logical(1))
  if (any(found)) 1 else 0.1
}

.empty_annotation_result <- function() {
  data.frame(
    observed_mz = numeric(), expected_mz = numeric(), ppm_error = numeric(),
    score = numeric(), mass_score = numeric(), rule_prior = numeric(),
    chemical_score = numeric(), isotope_score = numeric(),
    adduct_network_score = numeric(), adduct = character(), charge = integer(),
    formula = character(), neutral_mass = numeric(), metabolite_ids = character(),
    inchikeys = character(), metabolite_names = character(), ramp_ids = character(),
    stringsAsFactors = FALSE
  )
}

#' Query an indexed metabolite annotation search space
#'
#' Candidates are pruned by ppm tolerance and proton/charge validity during
#' index construction, then ranked using mass error, rule prior, mass-defect,
#' isotope, and contextual adduct-family evidence.
#'
#' @param observed_mz Numeric vector of observed m/z values.
#' @param index A `spamtp_mz_index` from [BuildMZAnnotationIndex()].
#' @param ppm Mass tolerance in parts per million.
#' @param ms1_spectrum Optional contextual spectrum with `mz` and `intensity`
#'   columns. It should represent the same retention-time window or spatial
#'   pixel/region as the queried peaks.
#' @param use_mass_defect Apply the CHO negative mass-defect penalty.
#' @param check_isotopes Score carbon-13 and, when relevant, chlorine/bromine
#'   M+2 patterns.
#' @param check_adduct_network Downweight complex ions whose base monomer
#'   family is absent from `ms1_spectrum`.
#' @param min_score Minimum final score to retain.
#'
#' @return A ranked candidate data frame.
#' @export
QueryMZAnnotationIndex <- function(observed_mz, index, ppm = 5,
                                   ms1_spectrum = NULL,
                                   use_mass_defect = TRUE,
                                   check_isotopes = TRUE,
                                   check_adduct_network = TRUE,
                                   min_score = 0) {
  if (!inherits(index, "spamtp_mz_index")) {
    stop("index must be created by BuildMZAnnotationIndex().")
  }
  observed_mz <- suppressWarnings(as.numeric(observed_mz))
  if (!length(observed_mz) || any(!is.finite(observed_mz) | observed_mz <= 0)) {
    stop("observed_mz must contain positive finite values.")
  }
  if (length(ppm) != 1 || !is.finite(ppm) || ppm < 0) {
    stop("ppm must be a single non-negative number.")
  }
  spectrum <- .normalise_ms1_spectrum(ms1_spectrum)
  tolerance <- ppm / 1e6
  results <- vector("list", length(observed_mz))

  for (i in seq_along(observed_mz)) {
    obs <- observed_mz[[i]]
    lower <- obs / (1 + tolerance)
    upper <- if (tolerance < 1) obs / (1 - tolerance) else Inf
    lower_open <- lower - .Machine$double.eps * max(1, abs(lower))
    lo <- findInterval(lower_open, index$expected_mz)
    hi <- findInterval(upper, index$expected_mz)
    if (hi <= lo) next
    hit <- seq.int(lo + 1L, hi)
    expected <- index$expected_mz[hit]
    error <- abs(obs - expected) / expected * 1e6
    keep <- is.finite(error) & error <= ppm
    if (!any(keep)) next
    hit <- hit[keep]
    expected <- expected[keep]
    error <- error[keep]
    compounds <- index$compounds[index$compound_idx[hit], , drop = FALSE]
    rules <- index$rules[index$rule_idx[hit], , drop = FALSE]

    mass_score <- if (ppm == 0) rep(1, length(error)) else {
      exp(-0.5 * (error / (ppm / 3))^2)
    }
    chemical_score <- if (isTRUE(use_mass_defect)) {
      vapply(seq_along(error), function(j) {
        .mass_defect_score(obs, compounds$formula[[j]], rules[j, , drop = FALSE])
      }, numeric(1))
    } else rep(1, length(error))
    isotope_score <- if (isTRUE(check_isotopes)) {
      vapply(seq_along(error), function(j) {
        .isotope_score(expected[[j]], compounds$formula[[j]],
                       rules[j, , drop = FALSE], spectrum, ppm)
      }, numeric(1))
    } else rep(NA_real_, length(error))
    network_score <- if (isTRUE(check_adduct_network)) {
      vapply(seq_along(error), function(j) {
        .adduct_family_score(compounds$exactmass[[j]], rules[j, , drop = FALSE],
                             index$family_rules, spectrum, ppm)
      }, numeric(1))
    } else rep(NA_real_, length(error))

    evidence_isotope <- ifelse(is.na(isotope_score), 1, isotope_score)
    evidence_network <- ifelse(is.na(network_score), 1, network_score)
    final_score <- mass_score * rules$prior * chemical_score *
      evidence_isotope * evidence_network

    out <- data.frame(
      observed_mz = obs,
      expected_mz = expected,
      ppm_error = error,
      score = final_score,
      mass_score = mass_score,
      rule_prior = rules$prior,
      chemical_score = chemical_score,
      isotope_score = isotope_score,
      adduct_network_score = network_score,
      adduct = rules$name,
      charge = rules$charge,
      formula = compounds$formula,
      neutral_mass = compounds$exactmass,
      metabolite_ids = compounds$isomers,
      inchikeys = compounds$isomers_inchikey,
      metabolite_names = compounds$isomers_names,
      ramp_ids = compounds$ramp_ids,
      stringsAsFactors = FALSE
    )
    results[[i]] <- out[out$score >= min_score, , drop = FALSE]
  }

  result <- do.call(rbind, results)
  if (is.null(result) || !nrow(result)) return(.empty_annotation_result())
  result <- result[order(result$observed_mz, -result$score, result$ppm_error), , drop = FALSE]
  rownames(result) <- NULL
  result
}

#' Annotate one or more observed m/z values
#'
#' Convenience wrapper that builds an index when one is not supplied. Reuse a
#' pre-built index for repeated or large searches.
#'
#' @param observed_mz Numeric vector of observed m/z values.
#' @param db Metabolite database used when `index` is `NULL`.
#' @param index Optional pre-built `spamtp_mz_index`.
#' @param polarity Ion mode.
#' @param adducts Optional adduct subset.
#' @param rules Optional custom rule table.
#' @param ppm Mass tolerance in ppm.
#' @param ms1_spectrum Optional contextual spectrum.
#' @param ... Additional arguments passed to [QueryMZAnnotationIndex()].
#'
#' @return A ranked candidate data frame.
#' @export
AnnotateMZ <- function(observed_mz, db = NULL, index = NULL,
                       polarity = c("positive", "negative", "neutral"),
                       adducts = NULL, rules = NULL, ppm = 5,
                       ms1_spectrum = NULL, ...) {
  polarity <- match.arg(polarity)
  if (is.null(index)) {
    if (is.null(db)) stop("Supply either db or a pre-built index.")
    index <- BuildMZAnnotationIndex(
      db = db, polarity = polarity, adducts = adducts, rules = rules
    )
  }
  QueryMZAnnotationIndex(
    observed_mz = observed_mz, index = index, ppm = ppm,
    ms1_spectrum = ms1_spectrum, ...
  )
}
