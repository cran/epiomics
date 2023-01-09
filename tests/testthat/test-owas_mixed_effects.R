
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
    omics = colnames_omic_fts,
    random_effects = "1|case_control_set",
    covars = c("age"),
    var_exposure_or_outcome = "exposure",
    family = "gaussian",
    conf_int = TRUE)
  
  ## without covars and with conf int------------
  owas_mixed_effects_out <- owas_mixed_effects(
    df = example_data,
    var = "cc_status",
    omics = colnames_omic_fts,
    random_effects = "1|case_control_set",
    var_exposure_or_outcome = "exposure",
    family = "gaussian",
    conf_int = TRUE)

  ## with covars, without conf int------------
  owas_mixed_effects_out <- owas_mixed_effects(
    df = example_data,
    var = "cc_status",
    omics = colnames_omic_fts,
    random_effects = "1|case_control_set",
    covars = c("age", "sex"),
    var_exposure_or_outcome = "exposure",
    family = "gaussian",
    conf_int = FALSE)
  
  ## without covars, without conf int------------
  owas_mixed_effects_out <- owas_mixed_effects(
    df = example_data,
    var = "cc_status",
    omics = colnames_omic_fts,
    random_effects = "1|case_control_set",
    var_exposure_or_outcome = "exposure",
    family = "gaussian",
    conf_int = FALSE, 
    REML = TRUE)
  
  
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
  ## with covars and conf int ---------------
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
  
})