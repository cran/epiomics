
# Test that owas_qgcomp  ----------------
test_that("owas_qgcomp works", {
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[
    grep("feature_",
         colnames(example_data))][1:10]
  
  # Get names of categorical omics
  colnames_cat_out <- colnames(example_data)[
   grep("disease",
   colnames(example_data))][1]
  
  # Names of exposures in mixture
  exposure_names <- c("exposure1", "exposure2", "exposure3")
  
  # Run function without covariates
  out <- owas_qgcomp(df = example_data,
                     expnms = exposure_names,
                     omics = colnames_omic_fts,
                     q = 4,
                     confidence_level = 0.95)
  
  
  # Continuous omics, with covars ----
  out <- owas_qgcomp(df = example_data,
                     expnms = c("exposure1", "exposure2", "exposure3"),
                     covars = c("weight", "age", "sex"),
                     omics = colnames_omic_fts,
                     q = 4,
                     confidence_level = 0.95)
  
  # Test that function returns expected dimensions
  testthat::expect_equal(object = ncol(out), 
                         expected = 12)
  testthat::expect_equal(object = nrow(out), 
                         expected = length(colnames_omic_fts))
  
  
  # Categorical omics, with covars ----
  out_boot <- owas_qgcomp(df = example_data,
                     expnms = c("exposure1", "exposure2", "exposure3"),
                     covars = c("weight", "age", "sex"),
                     omics = colnames_cat_out,
                     q = 4,
                     confidence_level = 0.95, 
                     family = "binomial", 
                     run.qgcomp.boot = TRUE)
  
  out_noboot <- owas_qgcomp(df = example_data,
                     expnms = c("exposure1", "exposure2", "exposure3"),
                     covars = c("weight", "age", "sex"),
                     omics = colnames_cat_out,
                     q = 4,
                     confidence_level = 0.95, 
                     family = "binomial", 
                     run.qgcomp.boot = FALSE)
  
})






# Test that owas_qgcomp errors  ----------------
testthat::test_that("owas_qgcomp errors", {
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
    owas_qgcomp(df = example_data,
                expnms = c("PFAS"),
                omics = colnames_omic_fts,
                q = 4,
                confidence_level = 0.95)
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 'Variable \'PFAS\' not found in data. Check data.')
  
  
  ## Different Var Types found in data ----
  error_message <- testthat::capture_error(
    owas_qgcomp(df = example_data,
                expnms = c("exposure1",
                           "exposure4"),
                omics = colnames_omic_fts,
                q = 4,
                confidence_level = 0.95)
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 'All variables in \'expnms\' must be the same type')
  
  
  ## must be numeric ----
  error_message <- testthat::capture_error(
    owas_qgcomp(df = example_data,
                expnms = c("exposure4"),
                omics = colnames_omic_fts,
                q = 4,
                confidence_level = 0.95)
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 'Currently exposures must be numeric, consider reformatting')
  
  
  ## Not all omics variables are found in the data ----
  error_message <- testthat::capture_error(
    owas_qgcomp(df = example_data,
                expnms = c("exposure1"),
                omics = c(colnames_omic_fts, "other"),
                q = 4,
                confidence_level = 0.95)
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 
      "Not all omics vars are found in the data. Check omics column names."
    )
  
  ## Not all covars  are found in the data ----
  error_message <- testthat::capture_error(
    owas_qgcomp(df = example_data,
                expnms = c("exposure1"),
                omics = c(colnames_omic_fts),
                covars = c("other"),
                q = 4)
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 
      "Not all covars are found in the data. Check covar column names."
  ) 
  
  ## family must be gaussian or binomial  ----
  error_message <- testthat::capture_error(
    owas_qgcomp(df = example_data,
                expnms = c("exposure1"),
                omics = c(colnames_omic_fts),
                family = "other",
                q = 4)
  )
  
  # Test error in data
  testthat::expect_equal(
    object = error_message$message, 
    expected = 
      "family must be either \"gaussian\" or \"binomial\" "
  ) 
  

  })
