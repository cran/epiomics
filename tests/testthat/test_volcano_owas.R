
# Test that owas works with a single var ----------------
test_that("volcano_owas works with a single var", {
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[
    grep("feature_",
         colnames(example_data))][1:5]
  
  # Run function with continuous exposure as the variable of interest
  owas_out <- owas(df = example_data,
                        var = "exposure1",
                        omics = colnames_omic_fts,
                        covars = c("age", "sex"),
                        var_exposure_or_outcome = "exposure",
                        family = "gaussian")
  
  vp <- volcano_owas(owas_out)
  
  testthat::expect_equal(object = vp$data$threshold, 
                         expected = c("Significant", rep("Non-significant", 4)))
  
})

# Test that owas works with multiple vars ----------------
test_that("volcano_owas works with multiple vars", {
  
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[
    grep("feature_",
         colnames(example_data))][1:5]
  
  # Run function with multiple exposures as the variable of interest
  expnms <- c("exposure1", "exposure2", "exposure3")
  
  # Continuous exposure, covars
  owas_multiple_vars <- owas(
    df = example_data,
    var = expnms,
    omics = colnames_omic_fts,
    covars = c("age", "sex"),
    var_exposure_or_outcome = "exposure",
    family = "gaussian")
  
  # Volcano plot
  vp <- volcano_owas(owas_multiple_vars, 
                     highlight_adj_p = FALSE)
  
  testthat::expect_equal(object = vp$data$threshold, 
                         expected = c("Significant", 
                                      rep("Non-significant", 14)))
 
})
