# Test that meet in middle works ----------------
test_that("meet_in_middle works", {
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[grep("feature_",
                                                   colnames(example_data))]
  
  # Run function with continuous exposure as the variable of interest
  cont_meet_in_middle <- meet_in_middle(df = example_data,
                                        exposure = "exposure1", 
                                        outcome = "disease1", 
                                        omics = colnames_omic_fts,
                                        # covars = c("age", "sex"), 
                                        outcome_family = "binomial", 
                                        confidence_level = 0.8)
  # Test that function returns expected dimensions ----
  testthat::expect_equal(object = ncol(cont_meet_in_middle$overlap), 
                         expected = 15)
  testthat::expect_equal(object = nrow(cont_meet_in_middle$overlap), 
                         expected = 4)
  
  #Test that this grabs the correct estimates: 
  testthat::expect_equal(
    object = cont_meet_in_middle$overlap$estimate_exp_omic[1],
    expected = as.numeric(coef(lm(feature_V47 ~ exposure1, 
                                  data = example_data))[2]))  
  
  #Test that this grabs the correct estimates: 
  testthat::expect_equal(
    object = cont_meet_in_middle$overlap$estimate_omic_out[1],
    expected = as.numeric(
      coef(glm(disease1 ~ feature_V47, 
               data = example_data, 
               family = binomial(link = "logit")))[2]))  
  
  
  # Run function with no overlapping features
  cont_meet_in_middle <- meet_in_middle(df = example_data,
                                        exposure = "exposure1", 
                                        outcome = "disease1", 
                                        omics = colnames_omic_fts[1:10],
                                        outcome_family = "binomial", 
                                        confidence_level = 0.8)
  
  testthat::expect_null(cont_meet_in_middle$overlap)
  
})

# Test meet in middle errors
test_that("meet_in_middle errors", {
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[grep("feature_",
                                                   colnames(example_data))]
  
  # Run function with continuous exposure as the variable of interest
  error_message <- testthat::capture_error(
    meet_in_middle(df = example_data,
                   exposure = c("exposure1", 
                                "exposure2"), 
                   outcome = "disease1", 
                   omics = colnames_omic_fts,
                   outcome_family = "binomial", 
                   confidence_level = 0.8)
  )
  # Test that function returns expected dimensions ----
  testthat::expect_equal(
    object = error_message$message, 
    expected = "More than one exposure or outcome is not currently supported.")
})