
# owas with multiple outcomes ----------------
test_that("owas works with multiple traits", {
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[
    grep("feature_",
         colnames(example_data))][1:5]
  
  # Run function with multiple exposures as the variable of interest
  expnms <- c("exposure1", "exposure2", "exposure3")
  
  # Continuous exposure, covars
  cont_owas_out <- owas(
    df = example_data,
    var = expnms,
    omics = colnames_omic_fts,
    covars = c("age", "sex"),
    var_exposure_or_outcome = "exposure",
    family = "gaussian")
  
  # Test that function returns expected dimensions
  testthat::expect_equal(object = ncol(cont_owas_out), 
                         expected = 8)
  testthat::expect_equal(object = nrow(cont_owas_out), 
                         expected = length(colnames_omic_fts)*length(expnms))
  
  # Continuous exposure, no covars
  owas(df = example_data,
       var = expnms,
       omics = colnames_omic_fts,
       # covars = c("age", "sex"),
       var_exposure_or_outcome = "exposure",
       family = "gaussian")
  
  # Continuous outcome, with covars
  owas(df = example_data,
       var = expnms,
       omics = colnames_omic_fts,
       covars = c("age", "sex"),
       var_exposure_or_outcome = "outcome",
       family = "gaussian")
  
  # Continuous outcome, no covars
  owas(df = example_data,
       var = expnms,
       omics = colnames_omic_fts,
       # covars = c("age", "sex"),
       var_exposure_or_outcome = "outcome",
       family = "gaussian")
  
})

# owas with single outcome ----------------
test_that("owas works with single continuous outcome", {
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[
    grep("feature_",
         colnames(example_data))][1:5]
  
  # Run function with continuous exposure as the variable of interest
  cont_owas_out <- owas(df = example_data,
                        var = "exposure1",
                        omics = colnames_omic_fts,
                        covars = c("age", "sex"),
                        var_exposure_or_outcome = "exposure",
                        family = "gaussian", 
                        conf_int = TRUE)
  
  # Test that function returns expected dimensions
  testthat::expect_equal(object = ncol(cont_owas_out), 
                         expected = 10)
  testthat::expect_equal(object = nrow(cont_owas_out), 
                         expected = 5)
  
  # Run function with dichotomous outcome as the variable of interest
  owas(df = example_data,
       var = "disease1",
       omics = colnames_omic_fts[1],
       covars = c("age", "sex"),
       var_exposure_or_outcome = "outcome",
       family = "binomial", 
       conf_int = TRUE)
  
  # Run function with dichotomous exposure as the variable of interest
  owas(df = example_data,
       var = "exposure4",
       omics = colnames_omic_fts[1],
       covars = c("age", "sex"),
       var_exposure_or_outcome = "outcome",
       ref_group = c("low"),
       family = "binomial") 
  
  
})



# owas gives correct errors ----------------
test_that("owas gives correct errors", {
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[
    grep("feature_",
         colnames(example_data))][1:5]
  
  # Run function with multiple exposures as the variable of interest
  disease_names <- c("diseasecat1", "diseasecat2", "diseasecat3") 
  
  ## Error that ref_group must be specified ----
  error_message <- testthat::capture_error(
    cont_owas_out <- owas(df = example_data,
                          var = disease_names,
                          omics = colnames_omic_fts,
                          covars = c("age", "sex"),
                          var_exposure_or_outcome = "exposure",
                          family = "gaussian")
  )
  
  # Test that function returns expected dimensions
  testthat::expect_equal(
    object = error_message$message, 
    expected = 
      'If var is character or factor, ref_group must be specified')
  
  
  # Error: two categories
  error_message <- testthat::capture_error(
    cont_owas_out <- owas(df = example_data,
                          var = c(disease_names),
                          omics = colnames_omic_fts,
                          covars = c("age", "sex"),
                          var_exposure_or_outcome = "exposure",
                          family = "gaussian", 
                          ref_group = "healthy")
  )
  
  # Test that function returns expected dimensions
  testthat::expect_equal(
    object = error_message$message, 
    expected = 'Var can only contain a maximum of two unique categories')
  
  # Error: same type
  error_message <- testthat::capture_error(
    owas(df = example_data,
         var = c(disease_names, "disease1"),
         omics = colnames_omic_fts,
         covars = c("age", "sex"),
         var_exposure_or_outcome = "exposure",
         family = "gaussian", 
         ref_group = "healthy")
  )
  
  # Test that function returns correct error
  testthat::expect_equal(
    object = error_message$message, 
    expected = 'All variables in \'var\' must be the same type')
  
  ## Var not found in data  ----
  error_message <- testthat::capture_error(
    owas(df = example_data,
         var = c(disease_names, "test"),
         omics = colnames_omic_fts,
         covars = c("age", "sex"),
         var_exposure_or_outcome = "exposure",
         family = "gaussian", 
         ref_group = "healthy")
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 'Variable \'test\' not found in data. Check data.')
  
  
  
  ## Not all omics variables are found in the data ----
  error_message <- testthat::capture_error(
    owas(df = example_data,
         var = c(disease_names[1:2]),
         omics = c(colnames_omic_fts, "test"),
         covars = c("age", "sex"),
         var_exposure_or_outcome = "exposure",
         family = "gaussian", 
         ref_group = "healthy")
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 
      "Not all omics vars are found in the data. Check omics column names."
  )
  
  ## Not all covars are found in the data ----
  error_message <- testthat::capture_error(
    owas(df = example_data,
         var = c(disease_names[1:2]),
         omics = c(colnames_omic_fts),
         covars = c("age", "test"),
         var_exposure_or_outcome = "exposure",
         family = "gaussian", 
         ref_group = "healthy")
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 
      "Not all covars are found in the data. Check covar column names."
  ) 
  
  ## Check that var_exposure_or_outcome is specified  ----
  error_message <- testthat::capture_error(
    owas(df = example_data,
         var = c("exposure1"),
         omics = colnames_omic_fts,
         covars = c("age", "sex"),
         var_exposure_or_outcome = "other",
         family = "gaussian", 
         ref_group = "healthy")
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 'var_exposure_or_outcome must be either "exposure" or "outcome" ')
  
  ## Check that family is specified  ----
  error_message <- testthat::capture_error(
    owas(df = example_data,
         var = c("exposure1"),
         omics = colnames_omic_fts,
         covars = c("age", "sex"),
         var_exposure_or_outcome = "exposure",
         family = "test", 
         ref_group = "healthy")
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 'family must be either "gaussian" or "binomial" ')
  
  
})
