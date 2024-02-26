
# Test that owas_mixed_effects works ----------------
test_that("owas_mixed_effects works", {
  
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[grep("feature_",
                                                   colnames(example_data))][1:5]
  
  
  # Continuous exp -------------
  ## with covars and conf int------------
  owas_mixed_effects_out <- owas_mixed_effects(
    df = example_data,
    var = "cc_status",
    omics = colnames_omic_fts[1],
    random_effects = "1|case_control_set",
    covars = c("age"),
    var_exposure_or_outcome = "exposure",
    family = "gaussian",
    conf_int = TRUE)
  
  ## without covars and with conf int------------
  owas_mixed_effects_out <- owas_mixed_effects(
    df = example_data,
    var = "cc_status",
    omics = colnames_omic_fts[1],
    random_effects = "1|case_control_set",
    var_exposure_or_outcome = "exposure",
    family = "gaussian",
    conf_int = TRUE)
  
  ## with covars, without conf int------------
  owas_mixed_effects_out <- owas_mixed_effects(
    df = example_data,
    var = "diseasecat1",
    omics = colnames_omic_fts[1],
    random_effects = "1|case_control_set",
    covars = c("age", "sex"),
    var_exposure_or_outcome = "outcome",
    family = "gaussian",
    ref_group = "healthy",
    conf_int = FALSE)
  
  ## without covars, without conf int------------
  # Dont change this one
  owas_mixed_effects_out <- owas_mixed_effects(
    df = example_data,
    var = "cc_status",
    omics = colnames_omic_fts[1],
    random_effects = "1|case_control_set",
    var_exposure_or_outcome = "exposure",
    family = "gaussian",
    conf_int = FALSE)
  
  
  # Test continuous analysis
  ### Run single model ---------------
  fit <- lmerTest::lmer(feature_V1 ~  cc_status + (1|case_control_set), 
                        data = example_data, REML = TRUE) |>
    summary() |> 
    coef()
  
  
  res_manual <- fit[nrow(fit),] |> t() |> as.data.frame()
  
  est_manual <- res_manual$Estimate
  
  ### Test that single model and mixed model give same estimates
  est_fxn <- owas_mixed_effects_out[
    owas_mixed_effects_out$feature_name == "feature_V1", 
  ]$estimate
  
  testthat::expect_equal(object = est_fxn, 
                         expected = est_manual)
  
  # Categorical exp -------------------------
  ## without covars and conf int ---------------
  owas_mixed_effects_out <- owas_mixed_effects(
    df = example_data,
    var = "disease1",
    omics = colnames_omic_fts[1:2],
    random_effects = "1|case_control_set",
    # covars = c("age", "sex"),
    var_exposure_or_outcome = "outcome", 
    family = "binomial",
    conf_int = TRUE)
  
  ### Run single model ---------------
  fit <- lme4::glmer(disease1 ~ feature_V1 + (1| case_control_set), 
                     data = example_data, 
                     family = binomial(link = "logit")) |>
    summary() |> 
    coef()
  
  res_manual <- fit[nrow(fit),] |> t() |> as.data.frame()
  
  est_manual <- res_manual$Estimate
  
  ### Test that single model and mixed model give same estimates
  est_fxn <- owas_mixed_effects_out[
    owas_mixed_effects_out$feature_name == "feature_V1", ]$estimate
  
  testthat::expect_equal(object = est_fxn, 
                         expected = est_manual)
  
  ## with covars and with conf int ---------------
  owas_mixed_effects_out <- owas_mixed_effects(
    df = example_data,
    var = "disease1",
    omics = colnames_omic_fts[1],
    random_effects = "1|case_control_set",
    covars = c("age", "sex"),
    var_exposure_or_outcome = "outcome", 
    family = "binomial",
    conf_int = TRUE)
  
  
  ## without covars, without conf int ---------------
  owas_mixed_effects_out <- owas_mixed_effects(
    df = example_data,
    var = "disease1",
    omics = colnames_omic_fts[1],
    random_effects = "1|case_control_set",
    # covars = c("age", "sex"),
    var_exposure_or_outcome = "outcome", 
    family = "binomial",
    conf_int = FALSE)
  ## without covars, with conf int ---------------
  owas_mixed_effects_out <- owas_mixed_effects(
    df = example_data,
    var = "disease1",
    omics = colnames_omic_fts[1],
    random_effects = "1|case_control_set",
    # covars = c("age", "sex"),
    var_exposure_or_outcome = "outcome", 
    family = "binomial",
    conf_int = TRUE)
  
  
})



# fxn gives correct errors ----------------
test_that("owas_ gives correct errors", {
  # Load Example Data
  data("example_data")
  example_data$disease_3cat <- c(example_data$diseasecat1[-1], "other")
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[
    grep("feature_",
         colnames(example_data))][1:5]
  
  # Run function with multiple exposures as the variable of interest
  disease_names <- c("diseasecat1", "diseasecat2", "diseasecat3") 
  
  
  
  ## Error that ref_group must be specified ----
  error_message <- testthat::capture_error(
    owas_mixed_effects(
      df = example_data,
      var = c("test"),
      omics = colnames_omic_fts[1],
      random_effects = "1|case_control_set",
      var_exposure_or_outcome = "outcome", 
      family = "binomial",
      conf_int = TRUE)
  )
  
  # Test that function returns expected dimensions
  testthat::expect_equal(
    object = error_message$message, 
    expected = 
      'Variable \'test\' not found in data. Check data.')
  
  
  # Error: two categories
  error_message <- testthat::capture_error(
    owas_mixed_effects(
      df = example_data,
      var = c("diseasecat1", "disease1"),
      omics = colnames_omic_fts[1],
      random_effects = "1|case_control_set",
      var_exposure_or_outcome = "outcome", 
      family = "binomial",
      conf_int = TRUE)
  )
  
  # Test that function returns expected dimensions
  testthat::expect_equal(
    object = error_message$message, 
    expected = 'All variables in \'var\' must be the same type')
  
  # Error: ref group
  error_message <- testthat::capture_error(
    owas_mixed_effects(
      df = example_data,
      var = c("diseasecat1"),
      omics = colnames_omic_fts[1],
      random_effects = "1|case_control_set",
      var_exposure_or_outcome = "outcome", 
      family = "binomial",
      conf_int = TRUE)
  )
  
  # Test that function returns correct error
  testthat::expect_equal(
    object = error_message$message, 
    expected = 'If var is character or factor, ref_group must be specified')
  
  ## Var not found in data  ----
  error_message <- testthat::capture_error(
    owas_mixed_effects(
      df = example_data,
      var = c("disease_3cat"),
      omics = colnames_omic_fts[1],
      random_effects = "1|case_control_set",
      var_exposure_or_outcome = "outcome", 
      ref_group = "disease",
      family = "binomial",
      conf_int = TRUE)
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 'Var can only contain a maximum of two unique categories')
  
  
  
  ## Not all omics variables are found in the data ----
  error_message <- testthat::capture_error(
    owas_mixed_effects(
      df = example_data,
      var = c("disease1"),
      omics = c("other"),
      random_effects = "1|case_control_set",
      var_exposure_or_outcome = "outcome", 
      ref_group = "disease",
      family = "binomial",
      conf_int = TRUE)
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 
      "Not all omics vars are in the data. Check omics column names."
  )
  
  ## Not all covars are found in the data ----
  error_message <- testthat::capture_error(
    owas_mixed_effects(
      df = example_data,
      var = c("disease1"),
      omics = c(colnames_omic_fts[1]),
      covars = c("test"),
      random_effects = "1|case_control_set",
      var_exposure_or_outcome = "outcome", 
      ref_group = "disease",
      family = "binomial",
      conf_int = TRUE))
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 
      "Not all covariates are in the data. Check covariate column names.") 
  
  ## Check that random effects exist  ----
  error_message <- testthat::capture_error(
    owas_mixed_effects(
      df = example_data,
      var = c("disease1"),
      omics = c(colnames_omic_fts[1]),
      random_effects = "1|test",
      var_exposure_or_outcome = "outcome", 
      ref_group = "disease",
      family = "binomial",
      conf_int = TRUE)
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 'Not all random effect vars are found in data')
  
  ## Check that family is specified ----
  error_message <- testthat::capture_error(
    owas_mixed_effects(
      df = example_data,
      var = c("disease1"),
      omics = c(colnames_omic_fts[1]),
      random_effects = "1|case_control_set",
      var_exposure_or_outcome = "outcome", 
      ref_group = "disease",
      family = "test",
      conf_int = TRUE)
  )
  
  testthat::expect_equal(
    object = error_message$message, 
    expected = 'family must be either \"gaussian\" or \"binomial\" ')
  
  ## Check that family is specified ----
  error_message <- testthat::capture_error(
    owas_mixed_effects(
      df = example_data,
      var = c("disease1"),
      omics = c(colnames_omic_fts[1]),
      random_effects = "1|case_control_set",
      var_exposure_or_outcome = "test", 
      family = "gaussian",
      conf_int = TRUE)
  )
  
  testthat::expect_equal(
    object = error_message$message, 
    expected = 'var_exposure_or_outcome must be either "exposure" or "outcome" ')
  
})
