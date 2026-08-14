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

test_that("SMILES decomposition recognises explicit metabolite functional groups", {
  structures <- DeconvolveSMILES(c(
    lactic_acid = "CC(O)C(=O)O",
    pyruvic_acid = "CC(=O)C(=O)O",
    acetone = "CC(=O)C",
    dopamine = "NCCc1ccc(O)c(O)c1",
    methyl_acetate = "CC(=O)OC",
    two_hg = "O[C@H](CCC(O)=O)C(O)=O"
  ), strict = TRUE)

  expect_true(all(structures$structure_valid))
  expect_equal(structures$carboxyl_sites[c(1, 2, 6)], c(1, 1, 2))
  expect_equal(structures$ketone_sites[c(2, 3)], c(1, 1))
  expect_equal(structures$carboxyl_sites[[5]], 0)
  expect_equal(structures$hydroxyl_sites[[5]], 0)
  expect_equal(structures$primary_amine_sites[[4]], 1)
  expect_equal(structures$phenolic_hydroxyl_sites[[4]], 2)
  expect_equal(structures$catechol_sites[[4]], 1)
  expect_equal(structures$fmp10_reactive_sites[[4]], 3)
  expect_match(structures$structure_evidence[[6]], "carboxyl=2")
})

test_that("ion-mode and alkali priors are inferred before mass matching", {
  target <- data.frame(
    formula = c("C5H8O5", "C3H6O"),
    exactmass = c(148.037173366, 58.041864812),
    name = c("2HG", "acetone"),
    iso_smiles = c("O[C@H](CCC(O)=O)C(O)=O", "CC(=O)C")
  )
  positive <- BuildMZAnnotationIndex(
    target, polarity = "positive", adducts = c("M+H", "M+Na"),
    collapse_isomers = FALSE
  )
  two_hg <- positive$compounds$isomers_names == "2HG"
  expect_true(positive$compounds$alkali_affinity_score[two_hg] >
                positive$compounds$positive_mode_score[two_hg])

  result <- QueryMZAnnotationIndex(positive$expected_mz, positive, ppm = 0)
  result_2hg <- result[result$metabolite_names == "2HG", ]
  expect_gt(
    result_2hg$structure_score[result_2hg$adduct == "M+Na"],
    result_2hg$structure_score[result_2hg$adduct == "M+H"]
  )
  expect_true(all(result$structure_status %in%
                    c("structure-supported", "structure-weak")))
  expect_true(all(nzchar(result$structure_rule_evidence)))

  negative <- BuildMZAnnotationIndex(
    target, polarity = "negative", adducts = c("M-H", "M-2H", "M-3H"),
    collapse_isomers = FALSE
  )
  indexed <- data.frame(
    name = negative$compounds$isomers_names[negative$compound_idx],
    adduct = negative$rules$name[negative$rule_idx]
  )
  expect_setequal(indexed$adduct[indexed$name == "2HG"], c("M-H", "M-2H"))
  expect_false(any(indexed$name == "acetone"))

  neutral <- BuildMZAnnotationIndex(
    target, polarity = "neutral", collapse_isomers = FALSE
  )
  neutral_result <- QueryMZAnnotationIndex(
    neutral$expected_mz, neutral, ppm = 0
  )
  expect_true(all(neutral_result$structure_score == 1))
})

test_that("SMILES-derived reactive sites automate FMP10 rule eligibility", {
  dopamine <- data.frame(
    formula = "C8H11NO2", exactmass = 153.078978601,
    name = "dopamine", iso_smiles = "NCCc1ccc(O)c(O)c1"
  )
  index <- BuildMZAnnotationIndex(
    dopamine, polarity = "positive", maldi_matrix = "FMP-10",
    adducts = "M+2FMP10a"
  )
  expect_equal(index$compounds$fmp10_reactive_sites, 3)
  expect_equal(length(index$expected_mz), 1)
  result <- QueryMZAnnotationIndex(index$expected_mz, index, ppm = 0)
  expect_equal(result$reactive_site_status, "verified")
  expect_match(result$structure_evidence, "phenolic-OH=2")
})

test_that("structure-only adduct prediction is mode-specific and auditable", {
  two_hg <- "O[C@H](CCC(O)=O)C(O)=O"
  positive <- PredictAdductsFromSMILES(two_hg, polarity = "positive")
  negative <- PredictAdductsFromSMILES(two_hg, polarity = "negative")
  neutral <- PredictAdductsFromSMILES(two_hg, polarity = "neutral")

  expect_true(positive$retained[positive$adduct == "M+Na"])
  expect_gt(
    positive$structure_score[positive$adduct == "M+Na"],
    positive$structure_score[positive$adduct == "M+H"]
  )
  expect_true(negative$retained[negative$adduct == "M-H"])
  expect_true(negative$retained[negative$adduct == "M-2H"])
  expect_false(negative$retained[negative$adduct == "M-3H"])
  expect_equal(neutral$adduct, "M")
  expect_equal(neutral$structure_score, 1)

  acetone <- PredictAdductsFromSMILES("CC(=O)C", polarity = "negative")
  expect_false(acetone$retained[acetone$adduct == "M-H"])
  expect_equal(
    acetone$proton_bound_status[acetone$adduct == "M-H"], "insufficient"
  )

  dopamine <- PredictAdductsFromSMILES(
    "NCCc1ccc(O)c(O)c1", maldi_matrix = "FMP-10"
  )
  expect_true(dopamine$retained[dopamine$adduct == "M+2FMP10a"])
  expect_equal(
    dopamine$reactive_site_status[dopamine$adduct == "M+2FMP10a"],
    "verified"
  )
})

test_that("weak alcohol deprotonation is scored without overriding curated bounds", {
  glucose <- data.frame(
    formula = "C6H12O6", exactmass = 180.063388,
    smiles = "OCC(O)C(O)C(O)CO"
  )
  inferred <- AnnotateSMILESStructure(glucose)
  expect_equal(inferred$max_exchangeable_protons, 0)
  expect_equal(inferred$proton_bound_source, "SMILES-inferred")
  expect_length(
    BuildMZAnnotationIndex(
      inferred, polarity = "negative", adducts = "M-H"
    )$expected_mz,
    1
  )

  curated <- glucose
  curated$max_protons <- 0
  expect_length(
    BuildMZAnnotationIndex(
      curated, polarity = "negative", adducts = "M-H"
    )$expected_mz,
    0
  )
})

test_that("registered reactive matrices report SMILES target compatibility", {
  dopamine <- "NCCc1ccc(O)c(O)c1"
  acetone <- "CC(=O)C"
  two_hg <- "O[C@H](CCC(O)=O)C(O)=O"

  dpp <- suppressWarnings(PredictAdductsFromSMILES(
    dopamine, maldi_matrix = "DPP-TFB"
  ))
  dnph <- suppressWarnings(PredictAdductsFromSMILES(
    acetone, maldi_matrix = "DNPH"
  ))
  picolyl <- suppressWarnings(PredictAdductsFromSMILES(
    two_hg, maldi_matrix = "2-picolylamine"
  ))
  tahs <- suppressWarnings(PredictAdductsFromSMILES(
    dopamine, maldi_matrix = "TAHS"
  ))
  conventional <- PredictAdductsFromSMILES(
    two_hg, maldi_matrix = "9-AA"
  )

  expect_true(all(dpp$matrix_target_status == "compatible-target"))
  expect_true(all(dnph$matrix_target_status == "compatible-target"))
  expect_true(all(picolyl$matrix_target_status == "compatible-target"))
  expect_true(all(tahs$matrix_target_status == "compatible-target"))
  expect_true(all(conventional$matrix_target_status == "not_reactive"))
})
