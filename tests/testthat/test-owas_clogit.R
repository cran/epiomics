
# Test owas_clogit ----------------
test_that("owas_clogit works", {
  
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[grep("feature_",
                                                colnames(example_data))]
  
  # No covars
  owas_clogit_out <- owas_clogit(
    df = example_data,
    cc_status = "cc_status",
    cc_set = "case_control_set", 
    omics = colnames_omic_fts, 
    covars = NULL,
    confidence_level = 0.95, 
    conf_int = FALSE)
  
  # Test that function returns expected dimensions
  testthat::expect_equal(object = ncol(owas_clogit_out), 
                         expected = 8)
  testthat::expect_equal(object = nrow(owas_clogit_out), 
                         expected = length(colnames_omic_fts))
  
  
  
  # One Covar
  owas_clogit_out <- owas_clogit(
    df = example_data,
    cc_status = "cc_status",
    cc_set = "case_control_set", 
    omics = colnames_omic_fts, 
    covars = "age")
  
  # With Covariates
  owas_clogit_out <- owas_clogit(
    df = example_data,
    cc_status = "cc_status",
    cc_set = "case_control_set", 
    omics = colnames_omic_fts, 
    covars = "age", 
    conf_int = TRUE)
  
})




# Test errors  ----------------
testthat::test_that("owas_clogit errors", {
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[
    grep("feature_",
         colnames(example_data))][1:10]
  
  # Names of exposures in mixture
  exposure_names <- c("exposure1", "exposure2", "exposure3")
  
  ## Exposure found in data ----
  error_message <- testthat::capture_error(
    owas_clogit(df = example_data,
                cc_status = "cc_status_wrong",
                cc_set = "case_control_set", 
                omics = colnames_omic_fts, 
                covars = NULL,
                confidence_level = 0.95, 
                conf_int = FALSE)
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 'Variable \'cc_status_wrong\' not found in data. Check data.')
  
  
  ## Not all omics variables are found in the data ----
  error_message <- testthat::capture_error(
    owas_clogit(df = example_data,
                cc_status = "cc_status",
                cc_set = "case_control_set", 
                omics = c(colnames_omic_fts, "other"),
                covars = NULL,
                confidence_level = 0.95, 
                conf_int = FALSE)
    )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 
      "Not all omics vars are found in the data. Check omics column names."
  )
  
  ## Not all covars  are found in the data ----
  error_message <- testthat::capture_error(
    owas_clogit(df = example_data,
                cc_status = "cc_status",
                cc_set = "case_control_set", 
                omics = colnames_omic_fts,
                covars = c("wrong_covar"),
                confidence_level = 0.95, 
                conf_int = FALSE)
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 
      "Not all covars are found in the data. Check covar column names."
  ) 
  
})
