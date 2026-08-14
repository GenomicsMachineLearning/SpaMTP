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
    carbon = carbon,
    hydrogen = hydrogen,
    nitrogen = nitrogen,
    oxygen = oxygen,
    fluorine = fluorine,
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

#' List MALDI matrix and on-tissue derivatization profiles
#'
#' The registry separates conventional MALDI matrices, reactive matrices, and
#' derivatization reagents. A registered profile does not necessarily imply
#' that a universal mass-shift rule is available: some reactions require a
#' coupling reagent, produce several products, or remain strongly dependent on
#' the analyte structure and acquisition conditions.
#'
#' @param matrix Optional matrix/reagent name or alias. When `NULL`, return the
#'   complete registry.
#'
#' @return A data frame describing supported matrix profiles and their current
#'   automatic-rule status.
#' @export
MALDIMatrixProfiles <- function(matrix = NULL) {
  profiles <- data.frame(
    matrix = c(
      "none", "dhb", "chca", "9aa", "dan", "norharmane",
      "fmp10", "fmp8", "fmp9", "dpp_tfb", "tmp_tfb", "n_mepyba",
      "dnph", "coniferyl_aldehyde", "dhba", "dhap", "girard_t",
      "girard_p", "2_picolylamine", "tmpa", "ampp_hatu", "tahs"
    ),
    display_name = c(
      "No specified MALDI matrix", "2,5-Dihydroxybenzoic acid (DHB)",
      "alpha-Cyano-4-hydroxycinnamic acid (CHCA)",
      "9-Aminoacridine (9-AA)", "1,5-Diaminonaphthalene (DAN)",
      "Norharmane", "FMP-10", "FMP-8", "FMP-9",
      "2,4-Diphenylpyrylium tetrafluoroborate (DPP-TFB)",
      "2,4,6-Trimethylpyrylium tetrafluoroborate (TMP-TFB)",
      "N-Methylpyridinium boronic acid (N-MePyBA)",
      "2,4-Dinitrophenylhydrazine (DNPH)", "Coniferyl aldehyde (CA)",
      "2,4-Dihydroxybenzaldehyde (DHBA)",
      "2,5-Dihydroxyacetophenone (DHAP)", "Girard reagent T",
      "Girard reagent P", "2-Picolylamine (2-PA)",
      "N,N,N-Trimethyl-2-(piperazin-1-yl)ethanaminium (TMPA)",
      "AMPP/HATU", "TAHS"
    ),
    category = c(
      "unspecified", rep("conventional_matrix", 5),
      rep("reactive_matrix", 10), rep("otcd_reagent", 6)
    ),
    default_polarity = c(
      "positive", "positive", "positive", "negative", "negative", "both",
      "positive", "positive", "positive", "positive", "positive",
      "positive", "positive", "positive", "positive", "positive",
      "positive", "positive", "positive", "positive", "positive", "positive"
    ),
    target_groups = c(
      "none", "none", "none", "none", "none", "none",
      "primary/secondary amine; phenolic hydroxyl",
      "primary/secondary amine; phenolic hydroxyl",
      "primary/secondary amine; phenolic hydroxyl", "primary amine",
      "primary amine", "catechol/1,2-diol", "aldehyde/ketone",
      "primary amine", "primary amine", "primary amine",
      "aldehyde/ketone", "aldehyde/ketone", "carboxylic acid",
      "carboxylic acid", "carboxylic acid; aldehyde", "catecholamine"
    ),
    automatic_rules = c(
      "standard", "validated_matrix_adduct", "validated_matrix_adduct",
      "standard", "standard", "standard", "validated_reactive_product",
      rep("profile_only", 15)
    ),
    reference = c(
      NA_character_, "10.1021/acs.analchem.0c04720",
      "10.1021/acs.analchem.0c04720", "10.1002/rcm.750",
      NA_character_, NA_character_, "10.1038/s41592-019-0551-3",
      "10.1038/s41592-019-0551-3", "10.1038/s41592-019-0551-3",
      "10.1007/s13361-015-1119-9", NA_character_,
      "10.1021/acs.analchem.8b03746", NA_character_, NA_character_, NA_character_,
      NA_character_, NA_character_, NA_character_,
      "10.1021/acs.analchem.6b01021", "10.1021/acs.analchem.0c02303",
      "10.1021/jasms.2c00336", NA_character_
    ),
    stringsAsFactors = FALSE
  )

  if (is.null(matrix)) return(profiles)
  canonical <- vapply(matrix, .normalise_maldi_matrix, character(1))
  profiles[match(canonical, profiles$matrix), , drop = FALSE]
}

.normalise_maldi_matrix <- function(matrix) {
  if (is.null(matrix) || !length(matrix)) return("none")
  if (length(matrix) != 1L || is.na(matrix) || !nzchar(trimws(matrix))) {
    stop("maldi_matrix must be one non-empty matrix/reagent name.")
  }
  key <- tolower(gsub("[^[:alnum:]]", "", matrix))
  aliases <- c(
    none = "none", unspecified = "none", nomatrix = "none",
    dhb = "dhb", `25dhb` = "dhb", dihydroxybenzoicacid = "dhb",
    chca = "chca", hcca = "chca", alphacyano4hydroxycinnamicacid = "chca",
    `9aa` = "9aa", `9aminoacridine` = "9aa",
    dan = "dan", `15dan` = "dan", diaminonaphthalene = "dan",
    norharmane = "norharmane", fmp10 = "fmp10", fmp8 = "fmp8",
    fmp9 = "fmp9", dpptfb = "dpp_tfb", dpp = "dpp_tfb",
    tmptfb = "tmp_tfb", tmp = "tmp_tfb", nmepyba = "n_mepyba",
    dnph = "dnph", coniferylaldehyde = "coniferyl_aldehyde",
    ca = "coniferyl_aldehyde", dhba = "dhba", dhap = "dhap",
    girardt = "girard_t", gt = "girard_t", girardp = "girard_p",
    gp = "girard_p", `2picolylamine` = "2_picolylamine",
    `2pa` = "2_picolylamine", tmpa = "tmpa", ampphatu = "ampp_hatu",
    ampp = "ampp_hatu", tahs = "tahs"
  )
  if (!key %in% names(aliases)) {
    stop(
      "Unknown MALDI matrix/reagent '", matrix, "'. See MALDIMatrixProfiles()."
    )
  }
  unname(aliases[[key]])
}

.resolve_maldi_polarity <- function(polarity = NULL, maldi_matrix = NULL) {
  if (!is.null(polarity)) {
    return(match.arg(polarity, c("positive", "negative", "neutral")))
  }
  if (is.null(maldi_matrix)) return("positive")
  profile <- MALDIMatrixProfiles(maldi_matrix)
  selected <- profile$default_polarity[[1]]
  if (identical(selected, "both")) {
    stop(
      profile$display_name[[1]],
      " supports both ion modes; supply polarity explicitly."
    )
  }
  selected
}

.decorate_matrix_rules <- function(rules, matrix, rule_class, source,
                                   reactive_group = "none",
                                   min_reactive_sites = 0L,
                                   site_count_column = NA_character_) {
  n <- nrow(rules)
  rules$maldi_matrix <- rep_len(matrix, n)
  rules$rule_class <- rep_len(rule_class, n)
  rules$rule_source <- rep_len(source, n)
  rules$reactive_group <- rep_len(reactive_group, n)
  rules$min_reactive_sites <- rep_len(as.integer(min_reactive_sites), n)
  rules$site_count_column <- rep_len(site_count_column, n)
  rules
}

.matrix_specific_rules <- function(matrix, polarity) {
  m <- .spamtp_mass_constants()
  empty <- .decorate_matrix_rules(
    AdductRules(polarity)[0, , drop = FALSE], matrix,
    "matrix_profile", NA_character_
  )
  if (polarity != "positive") return(empty)

  if (matrix == "fmp10") {
    retained_tag <- 20 * m$carbon + 14 * m$hydrogen + m$nitrogen
    rules <- rbind(
      .make_adduct_rule(
        "M+FMP10", "positive", 1, 1, retained_tag,
        complexity = "matrix_derivative", prior = 0.8
      ),
      .make_adduct_rule(
        "M+2FMP10a", "positive", 1, 1,
        2 * retained_tag - (m$carbon + 3 * m$hydrogen),
        complexity = "matrix_derivative", prior = 0.75
      ),
      .make_adduct_rule(
        "M+2FMP10b", "positive", 1, 1,
        2 * retained_tag - m$proton,
        complexity = "matrix_derivative", prior = 0.7
      )
    )
    rules <- .decorate_matrix_rules(
      rules, matrix, "reactive_product",
      "10.1038/s41592-019-0551-3; 10.1002/cmtd.202500062",
      "FMP-10-reactive site", c(1L, 2L, 2L), "fmp10_reactive_sites"
    )
    return(rules)
  }

  if (matrix == "dhb") {
    dhb <- 7 * m$carbon + 6 * m$hydrogen + 4 * m$oxygen
    dehydrated_dhb <- dhb - m$water
    rules <- rbind(
      .make_adduct_rule(
        "M+(DHB-H2O)+H", "positive", 1, 1,
        dehydrated_dhb + m$proton, complexity = "matrix_adduct",
        base_adduct = "M+H", prior = 0.4
      ),
      .make_adduct_rule(
        "M+(DHB-H2O)+Na", "positive", 1, 1,
        dehydrated_dhb + m$sodium_ion, complexity = "matrix_adduct",
        base_adduct = "M+Na", contains_metal = TRUE, prior = 0.3
      ),
      .make_adduct_rule(
        "M+2(DHB-H2O)+H", "positive", 1, 1,
        2 * dehydrated_dhb + m$proton, complexity = "matrix_adduct",
        base_adduct = "M+H", prior = 0.2
      )
    )
    return(.decorate_matrix_rules(
      rules, matrix, "matrix_adduct", "10.1021/acs.analchem.0c04720",
      "amine-enriched; parent-ion support required"
    ))
  }

  if (matrix == "chca") {
    chca <- 10 * m$carbon + 7 * m$hydrogen + m$nitrogen + 3 * m$oxygen
    rule <- .make_adduct_rule(
      "M+CHCA+Na", "positive", 1, 1, chca + m$sodium_ion,
      complexity = "matrix_adduct", base_adduct = "M+Na",
      contains_metal = TRUE, prior = 0.25
    )
    return(.decorate_matrix_rules(
      rule, matrix, "matrix_adduct", "10.1021/acs.analchem.0c04720",
      "amine-enriched; parent-ion support required"
    ))
  }
  empty
}

#' Select annotation rules from the MALDI matrix/reagent profile
#'
#' Conventional matrices select ordinary ion rules and, where validated,
#' matrix-metabolite adducts. Reactive matrices additionally select covalent
#' product rules. Profiles without a verified universal net transformation are
#' retained in the registry but do not generate speculative mass shifts.
#'
#' @param maldi_matrix Matrix or derivatization reagent name; aliases such as
#'   `"2,5-DHB"`, `"HCCA"`, `"9-AA"`, and `"FMP-10"` are accepted.
#' @param polarity `"positive"`, `"negative"`, or `"neutral"`. If `NULL`, use
#'   the profile default; profiles supporting both modes require an explicit
#'   choice.
#' @param include_standard Include the validated ordinary adduct rules.
#' @param include_matrix_products Include validated matrix-adduct or reactive
#'   product rules.
#'
#' @return A validated adduct/reaction rule data frame accepted by
#'   [BuildMZAnnotationIndex()].
#' @export
MALDIMatrixRules <- function(maldi_matrix, polarity = NULL,
                             include_standard = TRUE,
                             include_matrix_products = TRUE) {
  matrix <- .normalise_maldi_matrix(maldi_matrix)
  profile <- MALDIMatrixProfiles(matrix)
  polarity <- .resolve_maldi_polarity(polarity, matrix)

  standard <- .decorate_matrix_rules(
    AdductRules(polarity)[if (isTRUE(include_standard)) TRUE else FALSE, , drop = FALSE],
    matrix, "standard_adduct", "SpaMTP validated general adduct rules"
  )
  specific <- if (isTRUE(include_matrix_products)) {
    .matrix_specific_rules(matrix, polarity)
  } else {
    standard[0, , drop = FALSE]
  }

  if (isTRUE(include_matrix_products) &&
      profile$automatic_rules[[1]] == "profile_only") {
    warning(
      "The ", profile$display_name[[1]],
      " profile is registered, but SpaMTP does not yet bundle a universal ",
      "validated net product-mass rule. Standard adduct rules were selected; ",
      "supply a study-specific custom rule table when appropriate.",
      call. = FALSE
    )
  }
  .validate_adduct_rules(rbind(standard, specific))
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
  optional_defaults <- list(
    maldi_matrix = "unspecified",
    rule_class = "custom_adduct",
    rule_source = "user-supplied rule",
    reactive_group = "none",
    min_reactive_sites = 0L,
    site_count_column = NA_character_
  )
  for (column in names(optional_defaults)) {
    if (!column %in% names(rules)) {
      rules[[column]] <- rep_len(optional_defaults[[column]], nrow(rules))
    }
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
  rules$min_reactive_sites <- suppressWarnings(
    as.integer(rules$min_reactive_sites)
  )
  if (any(is.na(rules$min_reactive_sites) | rules$min_reactive_sites < 0L)) {
    stop("min_reactive_sites must contain non-negative integers.")
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
  site_columns <- c(
    "fmp10_reactive_sites", "primary_amine_sites", "secondary_amine_sites",
    "phenolic_hydroxyl_sites", "catechol_sites", "carbonyl_sites",
    "carboxyl_sites"
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
  for (column in site_columns) {
    compounds[[column]] <- if (column %in% names(db)) {
      suppressWarnings(as.numeric(db[[column]]))
    } else {
      rep(NA_real_, n)
    }
  }
  keep <- is.finite(compounds$exactmass) & compounds$exactmass > 0 &
    !is.na(compounds$formula) & nzchar(compounds$formula)
  compounds <- compounds[keep, , drop = FALSE]

  if (isTRUE(collapse_isomers) && nrow(compounds)) {
    group_key <- paste(
      compounds$formula,
      format(compounds$exactmass, digits = 17, scientific = FALSE, trim = TRUE),
      ifelse(is.na(compounds$max_exchangeable_protons), "NA",
             compounds$max_exchangeable_protons),
      do.call(
        paste,
        c(
          lapply(compounds[site_columns], function(x) {
            ifelse(is.na(x), "NA", format(x, trim = TRUE))
          }),
          sep = ":"
        )
      ),
      sep = "\r"
    )
    if (anyDuplicated(group_key)) {
      grouped_rows <- split(seq_len(nrow(compounds)), group_key)
      first_row <- vapply(grouped_rows, function(i) i[[1]], integer(1))
      collapsed_site_values <- lapply(
        compounds[site_columns], function(x) x[first_row]
      )
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
      for (column in site_columns) {
        compounds[[column]] <- collapsed_site_values[[column]]
      }
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
#' @param polarity `"positive"`, `"negative"`, or `"neutral"`. When `NULL`,
#'   use the matrix-profile default, or positive mode when no profile is given.
#' @param adducts Optional character vector of adduct names or bracketed
#'   notations. When `NULL`, use the complete automatically selected rule
#'   space; this is a filter, not a compulsory input.
#' @param rules Optional custom rule table. See [AdductRules()].
#' @param maldi_matrix Optional MALDI matrix or derivatization reagent profile.
#'   When supplied and `rules` is `NULL`, [MALDIMatrixRules()] automatically
#'   selects standard and validated matrix-specific rules. `adducts` remains an
#'   optional filter on that selected rule space.
#' @param collapse_isomers Collapse records sharing formula, exact mass, and
#'   proton bound before indexing.
#'
#' @return An object of class `spamtp_mz_index`.
#' @export
BuildMZAnnotationIndex <- function(db, polarity = NULL,
                                   adducts = NULL, rules = NULL,
                                   maldi_matrix = NULL,
                                   collapse_isomers = TRUE) {
  polarity <- .resolve_maldi_polarity(polarity, maldi_matrix)
  if (is.null(rules)) {
    rules <- if (is.null(maldi_matrix)) {
      AdductRules(polarity = polarity)
    } else {
      MALDIMatrixRules(maldi_matrix, polarity = polarity)
    }
  } else {
    if (!is.null(maldi_matrix)) {
      stop(
        "Supply either maldi_matrix for automatic selection or an explicit ",
        "rules table, not both."
      )
    }
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

  reactive_columns <- unique(rules$site_count_column[
    rules$min_reactive_sites > 0L & !is.na(rules$site_count_column)
  ])
  missing_reactive_columns <- reactive_columns[vapply(
    reactive_columns,
    function(column) {
      !column %in% names(compounds) || !any(is.finite(compounds[[column]]))
    },
    logical(1)
  )]
  if (length(missing_reactive_columns)) {
    warning(
      "No usable structure-derived reaction-site counts were found in: ",
      paste(missing_reactive_columns, collapse = ", "), ". Reactive-product ",
      "candidates will be retained with status 'unknown' and downweighted; ",
      "add these columns to db for chemical eligibility pruning.",
      call. = FALSE
    )
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
    required_sites <- rule$min_reactive_sites[[1]]
    site_column <- rule$site_count_column[[1]]
    if (required_sites > 0L && !is.na(site_column) &&
        site_column %in% names(compounds)) {
      site_count <- compounds[[site_column]]
      known_site_count <- is.finite(site_count)
      valid[known_site_count] <- valid[known_site_count] &
        site_count[known_site_count] >= required_sites
    }
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
      polarity = polarity,
      maldi_matrix = if (is.null(maldi_matrix)) {
        "unspecified"
      } else {
        .normalise_maldi_matrix(maldi_matrix)
      }
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
  matrix_profile <- if (is.null(x$maldi_matrix)) "unspecified" else x$maldi_matrix
  cat("  MALDI matrix profile: ", matrix_profile, "\n", sep = "")
  cat("  metabolites: ", nrow(x$compounds), "\n", sep = "")
  cat("  ion/reaction rules: ", nrow(x$rules), "\n", sep = "")
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
  complex_types <- c("multimer", "solvent", "metal_exchange", "matrix_adduct")
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
    chemical_score = numeric(), reactive_site_score = numeric(),
    reactive_site_status = character(), isotope_score = numeric(),
    adduct_network_score = numeric(), adduct = character(), charge = integer(),
    maldi_matrix = character(), rule_class = character(),
    rule_source = character(), reactive_group = character(),
    min_reactive_sites = integer(),
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
    site_status <- vapply(seq_along(error), function(j) {
      required_sites <- rules$min_reactive_sites[[j]]
      site_column <- rules$site_count_column[[j]]
      if (required_sites <= 0L) return("not_required")
      if (is.na(site_column) || !site_column %in% names(compounds)) {
        return("unknown")
      }
      count <- compounds[[site_column]][[j]]
      if (!is.finite(count)) return("unknown")
      if (count >= required_sites) "verified" else "insufficient"
    }, character(1))
    reactive_site_score <- ifelse(site_status == "unknown", 0.25, 1)
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
      reactive_site_score *
      evidence_isotope * evidence_network

    out <- data.frame(
      observed_mz = obs,
      expected_mz = expected,
      ppm_error = error,
      score = final_score,
      mass_score = mass_score,
      rule_prior = rules$prior,
      chemical_score = chemical_score,
      reactive_site_score = reactive_site_score,
      reactive_site_status = site_status,
      isotope_score = isotope_score,
      adduct_network_score = network_score,
      adduct = rules$name,
      charge = rules$charge,
      maldi_matrix = rules$maldi_matrix,
      rule_class = rules$rule_class,
      rule_source = rules$rule_source,
      reactive_group = rules$reactive_group,
      min_reactive_sites = rules$min_reactive_sites,
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
#' @param polarity Ion mode. When `NULL`, infer it from `index` or the matrix
#'   profile, falling back to positive mode.
#' @param adducts Optional adduct subset. When `NULL`, retain the complete rule
#'   space selected by the matrix profile or polarity.
#' @param rules Optional custom rule table.
#' @param maldi_matrix Optional MALDI matrix/reagent profile used to select
#'   rules automatically when `rules` is `NULL`.
#' @param ppm Mass tolerance in ppm.
#' @param ms1_spectrum Optional contextual spectrum.
#' @param database_version SpaMTPdb/RaMP version used when `db = NULL`.
#' @param database_source Database source used when `db = NULL`; see
#'   [LoadSpaMTPDatabase()].
#' @param database_local_dir Optional staged SpaMTPdb resource directory.
#' @param ... Additional arguments passed to [QueryMZAnnotationIndex()].
#'
#' @return A ranked candidate data frame.
#' @export
AnnotateMZ <- function(observed_mz, db = NULL, index = NULL,
                       polarity = NULL,
                       adducts = NULL, rules = NULL, maldi_matrix = NULL,
                       ppm = 5,
                       ms1_spectrum = NULL,
                       database_version = "latest",
                       database_source = c("auto", "spamtpdb", "bundled"),
                       database_local_dir = NULL, ...) {
  if (!is.null(index) && !inherits(index, "spamtp_mz_index")) {
    stop("index must be created by BuildMZAnnotationIndex().")
  }
  polarity <- if (is.null(polarity) && inherits(index, "spamtp_mz_index")) {
    index$polarity
  } else {
    .resolve_maldi_polarity(polarity, maldi_matrix)
  }
  if (is.null(index)) {
    if (is.null(db)) {
      db <- .spamtp_db_resource(
        "chem_props",
        version = database_version,
        source = match.arg(database_source),
        local_dir = database_local_dir
      )
    }
    index <- BuildMZAnnotationIndex(
      db = db, polarity = polarity, adducts = adducts, rules = rules,
      maldi_matrix = maldi_matrix
    )
  } else if (!identical(index$polarity, polarity)) {
    stop("The supplied index polarity does not match polarity.")
  } else if (!is.null(maldi_matrix) &&
             !identical(index$maldi_matrix, .normalise_maldi_matrix(maldi_matrix))) {
    stop("The supplied index MALDI matrix profile does not match maldi_matrix.")
  }
  QueryMZAnnotationIndex(
    observed_mz = observed_mz, index = index, ppm = ppm,
    ms1_spectrum = ms1_spectrum, ...
  )
}
