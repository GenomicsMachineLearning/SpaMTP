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
