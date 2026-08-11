test_that("validated adduct rules are charge balanced", {
  rules <- AdductRules("both")
  charged <- rules$polarity != "neutral"
  expect_true(all(rules$added_charge[charged] - rules$loss_h[charged] ==
                    rules$charge[charged]))
  expect_true(all(rules$n_molecules >= 1))
  expect_true(all(rules$charge != 0))
})

test_that("2HG positive adduct masses use explicit ion masses", {
  target <- data.frame(
    formula = "C5H8O5",
    exactmass = 148.037173366,
    isomers = "HMDB0000606",
    isomers_inchikey = "HWXBTNAVRSUOJR-GSVOUGTGSA-N",
    isomers_names = "D-2-Hydroxyglutaric acid",
    max_exchangeable_protons = 2
  )
  index <- BuildMZAnnotationIndex(
    target,
    adducts = c("[M+H]+", "M+Na", "M+2Na-H", "2M+H")
  )
  expected <- c(
    "M+H" = 149.044449832621,
    "M+Na" = 171.026394068091,
    "M+2Na-H" = 193.008338303561,
    "2M+H" = 297.081623198621
  )
  actual <- setNames(
    index$expected_mz,
    index$rules$name[index$rule_idx]
  )
  expect_equal(actual[names(expected)], expected, tolerance = 1e-10)

  result <- QueryMZAnnotationIndex(unname(actual[names(expected)]), index, ppm = 0)
  expect_equal(nrow(result), 4)
  expect_true(all(result$ppm_error == 0))
})

test_that("exchangeable proton bounds prune impossible ion states", {
  target <- data.frame(
    formula = "C5H8O5",
    exactmass = 148.037173366,
    id = "2HG",
    name = "2HG",
    max_exchangeable_protons = 2
  )
  index <- BuildMZAnnotationIndex(
    target,
    polarity = "negative",
    adducts = c("M-H", "M-2H", "M-3H")
  )
  indexed_rules <- unique(index$rules$name[index$rule_idx])
  expect_setequal(indexed_rules, c("M-H", "M-2H"))
  expect_false("M-3H" %in% indexed_rules)
})

test_that("unbalanced custom rules fail before candidate generation", {
  target <- data.frame(formula = "C5H8O5", exactmass = 148.037173366)
  bad_rule <- AdductRules("positive")[1, ]
  bad_rule$added_charge <- 2L
  expect_error(
    BuildMZAnnotationIndex(target, rules = bad_rule),
    "Unbalanced adduct"
  )
})

test_that("binary interval query enforces ppm tolerance", {
  target <- data.frame(formula = "C5H8O5", exactmass = 148.037173366)
  index <- BuildMZAnnotationIndex(target, adducts = "M+H")
  expected <- index$expected_mz[[1]]
  inside <- expected * (1 + 4.9e-6)
  outside <- expected * (1 + 5.1e-6)
  expect_equal(nrow(QueryMZAnnotationIndex(inside, index, ppm = 5)), 1)
  expect_equal(nrow(QueryMZAnnotationIndex(outside, index, ppm = 5)), 0)
})

test_that("mass defect is a score penalty rather than an early rejection", {
  target <- data.frame(
    formula = "C2H2O2",
    exactmass = 99.942723533379,
    id = "synthetic"
  )
  index <- BuildMZAnnotationIndex(target, adducts = "M+H")
  result <- QueryMZAnnotationIndex(index$expected_mz, index, ppm = 0)
  expect_equal(nrow(result), 1)
  expect_equal(result$chemical_score, 0.2)
  expect_gt(result$score, 0)
})

test_that("isotope and adduct-family evidence changes candidate scores", {
  target <- data.frame(
    formula = "C5H8O5",
    exactmass = 148.037173366,
    id = "2HG",
    max_exchangeable_protons = 2
  )
  monomer_index <- BuildMZAnnotationIndex(target, adducts = "M+H")
  monomer_mz <- monomer_index$expected_mz[[1]]
  isotope_mz <- monomer_mz + 1.00335483507
  spectrum <- data.frame(
    mz = c(monomer_mz, isotope_mz),
    intensity = c(1000, 53.5)
  )
  with_isotope <- QueryMZAnnotationIndex(
    monomer_mz, monomer_index, ppm = 5, ms1_spectrum = spectrum
  )
  without_isotope <- QueryMZAnnotationIndex(
    monomer_mz, monomer_index, ppm = 5,
    ms1_spectrum = spectrum[1, , drop = FALSE]
  )
  expect_gt(with_isotope$isotope_score, without_isotope$isotope_score)

  complex_index <- BuildMZAnnotationIndex(target, adducts = "2M+H")
  complex_mz <- complex_index$expected_mz[[1]]
  isolated <- QueryMZAnnotationIndex(
    complex_mz, complex_index, ppm = 5,
    ms1_spectrum = data.frame(mz = complex_mz, intensity = 100)
  )
  supported <- QueryMZAnnotationIndex(
    complex_mz, complex_index, ppm = 5,
    ms1_spectrum = data.frame(
      mz = c(complex_mz, monomer_mz), intensity = c(100, 500)
    )
  )
  expect_equal(isolated$adduct_network_score, 0.1)
  expect_equal(supported$adduct_network_score, 1)
  expect_gt(supported$score, isolated$score)
})

test_that("legacy and RaMP chemical-property schemas are both accepted", {
  legacy <- data.frame(
    formula = "C5H8O5",
    exactmass = 148.037173366,
    isomers = "HMDB0000606",
    isomers_inchikey = "HWXBTNAVRSUOJR-GSVOUGTGSA-N",
    isomers_names = "D-2-Hydroxyglutaric acid"
  )
  ramp <- data.frame(
    ramp_id = "RAMP_C_000219574",
    chem_data_source = "hmdb",
    chem_source_id = "hmdb:HMDB0000606",
    inchi_key = "HWXBTNAVRSUOJR-GSVOUGTGSA-N",
    monoisotop_mass = 148.037173366,
    common_name = "D-2-Hydroxyglutaric acid",
    mol_formula = "C5H8O5"
  )
  legacy_result <- AnnotateMZ(149.044449832621, legacy, adducts = "M+H", ppm = 1e-6)
  ramp_result <- AnnotateMZ(149.044449832621, ramp, adducts = "M+H", ppm = 1e-6)
  expect_equal(legacy_result$metabolite_names, ramp_result$metabolite_names)
  expect_equal(ramp_result$ramp_ids, "RAMP_C_000219574")
})

test_that("legacy wrappers reuse an index and return zero-row misses safely", {
  target <- data.frame(
    formula = "C5H8O5",
    exactmass = 148.037173366,
    isomers = "HMDB0000606",
    isomers_inchikey = "HWXBTNAVRSUOJR-GSVOUGTGSA-N",
    isomers_names = "D-2-Hydroxyglutaric acid"
  )
  index <- BuildMZAnnotationIndex(target, adducts = "M+H")

  hit <- AnnotateBigData(
    index$expected_mz[[1]], index = index, ppm_error = 0, verbose = FALSE
  )
  miss <- AnnotateBigData(
    999999, index = index, ppm_error = 5, verbose = FALSE
  )

  expect_equal(nrow(hit), 1)
  expect_equal(nrow(miss), 0)
  expect_true(all(c("all_Scores", "all_Ramp_IDs", "mz_names") %in% names(miss)))
})

test_that("MALDI matrix aliases select rules without compulsory adduct input", {
  profiles <- MALDIMatrixProfiles()
  expect_true(all(c(
    "dhb", "chca", "9aa", "fmp10", "dpp_tfb", "girard_t", "ampp_hatu"
  ) %in% profiles$matrix))
  expect_equal(MALDIMatrixProfiles("2,5-DHB")$matrix, "dhb")
  expect_equal(MALDIMatrixProfiles("HCCA")$matrix, "chca")
  expect_equal(MALDIMatrixProfiles("9-AA")$default_polarity, "negative")

  target <- data.frame(
    formula = "C8H11NO2", exactmass = 153.078978601,
    name = "dopamine", fmp10_reactive_sites = 3
  )
  index <- BuildMZAnnotationIndex(
    target, polarity = "positive", maldi_matrix = "FMP-10"
  )
  expect_equal(index$maldi_matrix, "fmp10")
  expect_true(all(c("M+FMP10", "M+2FMP10a", "M+2FMP10b") %in%
                    index$rules$name))

  batch_hit <- AnnotateBigData(
    index$expected_mz[index$rules$name[index$rule_idx] == "M+2FMP10a"][[1]],
    db = target, maldi_matrix = "FMP-10", ppm_error = 0, verbose = FALSE
  )
  expect_match(batch_hit$all_Adducts, "M\\+2FMP10a")

  aa_index <- BuildMZAnnotationIndex(target, maldi_matrix = "9-AA")
  expect_equal(aa_index$polarity, "negative")
  expect_true("M-H" %in% aa_index$rules$name)
  expect_error(
    AnnotateMZ(100, index = aa_index, polarity = "positive"),
    "polarity does not match"
  )
})

test_that("FMP10 reaction sites prune impossible products and flag unknowns", {
  target <- data.frame(
    formula = c("C8H11NO2", "C8H11NO2"),
    exactmass = c(153.078978601, 153.078978601),
    name = c("eligible", "ineligible"),
    fmp10_reactive_sites = c(3, 1)
  )
  index <- BuildMZAnnotationIndex(
    target, polarity = "positive", maldi_matrix = "FMP-10",
    adducts = "M+2FMP10a", collapse_isomers = FALSE
  )
  expect_equal(length(index$expected_mz), 1)
  result <- QueryMZAnnotationIndex(index$expected_mz, index, ppm = 0)
  expect_equal(result$metabolite_names, "eligible")
  expect_equal(result$reactive_site_status, "verified")
  expect_equal(result$reactive_site_score, 1)

  unknown <- target[1, setdiff(names(target), "fmp10_reactive_sites")]
  expect_warning(
    unknown_index <- BuildMZAnnotationIndex(
      unknown, polarity = "positive", maldi_matrix = "FMP-10",
      adducts = "M+2FMP10a"
    ),
    "reaction-site counts"
  )
  unknown_result <- QueryMZAnnotationIndex(
    unknown_index$expected_mz, unknown_index, ppm = 0
  )
  expect_equal(unknown_result$reactive_site_status, "unknown")
  expect_equal(unknown_result$reactive_site_score, 0.25)
})

test_that("DHB and CHCA matrix adducts require parent-family support", {
  target <- data.frame(
    formula = "C10H16O4", exactmass = 200,
    id = "synthetic", name = "synthetic"
  )
  dhb_index <- BuildMZAnnotationIndex(
    target, polarity = "positive", maldi_matrix = "DHB",
    adducts = "M+(DHB-H2O)+H"
  )
  matrix_mz <- dhb_index$expected_mz[[1]]
  parent_mz <- 200 + .spamtp_mass_constants()$proton
  isolated <- QueryMZAnnotationIndex(
    matrix_mz, dhb_index, ppm = 1,
    ms1_spectrum = data.frame(mz = matrix_mz, intensity = 100)
  )
  supported <- QueryMZAnnotationIndex(
    matrix_mz, dhb_index, ppm = 1,
    ms1_spectrum = data.frame(
      mz = c(parent_mz, matrix_mz), intensity = c(500, 100)
    )
  )
  expect_equal(isolated$adduct_network_score, 0.1)
  expect_equal(supported$adduct_network_score, 1)
  expect_gt(supported$score, isolated$score)

  chca_rules <- MALDIMatrixRules("CHCA", polarity = "positive")
  expect_true("M+CHCA+Na" %in% chca_rules$name)
  aa_rules <- MALDIMatrixRules("9-AA")
  expect_false(any(aa_rules$rule_class == "matrix_adduct"))
  expect_true("M-H" %in% aa_rules$name)
})

test_that("profile-only reactive matrices do not invent product masses", {
  expect_warning(
    rules <- MALDIMatrixRules("DPP-TFB", polarity = "positive"),
    "does not yet bundle a universal"
  )
  expect_true(all(rules$rule_class == "standard_adduct"))
  expect_false(any(grepl("DPP", rules$name, fixed = TRUE)))
})
