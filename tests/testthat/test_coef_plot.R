
# Test that owas works with a single var ----------------
test_that("coef_plot works with a single var", {
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[
    grep("feature_",
         colnames(example_data))][1:5]
  
  # Run function with continuous exposure as the variable of interest
  coefplot <- owas(df = example_data,
                   var = "exposure1",
                   omics = colnames_omic_fts,
                   covars = c("age", "sex"),
                   var_exposure_or_outcome = "exposure",
                   family = "gaussian", 
                   conf_int = TRUE) |> 
    coef_plot_from_owas()
  
  testthat::expect_equal(object = coefplot$data$threshold, 
                         expected = c("Significant",
                                      rep("Non-significant", 4)))
  
  
  # Test error for no CI's
  owas(df = example_data,
       var = "exposure1",
       omics = colnames_omic_fts,
       covars = c("age", "sex"), 
       var_exposure_or_outcome = "exposure",
       family = "gaussian") |> 
    coef_plot_from_owas() |>
    testthat::capture_error() 
  
  # Check with multiple exposures
  coefplot <- owas(df = example_data,
       var = c("exposure1", "exposure2"),
       omics = colnames_omic_fts,
       covars = c("age", "sex"),
       var_exposure_or_outcome = "exposure",
       family = "gaussian", 
       conf_int = TRUE) |> 
    coef_plot_from_owas(main_cat_var = "feature_name") 
  
  
    testthat::expect_equal(object = coefplot$data$threshold, 
                           expected = c("Significant",
                                        rep("Non-significant", 9)))
  
})
