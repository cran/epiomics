# epiomics: Analysis of Omics Data in Observational Studies

<!-- badges: start -->

[![R-CMD-check](https://github.com/Goodrich-Lab/epiomics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Goodrich-Lab/epiomics/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/Goodrich-Lab/epiomics/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Goodrich-Lab/epiomics?branch=main)
<!-- badges: end -->

**epiomics** provides a collection of fast and flexible functions for the analysis of omics data in observational studies.

## Installation

You can install epiomics from [CRAN](https://cran.r-project.org/) with:

``` r
install.packages("epiomics")
```

You can download the developmental version of epiomics from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Goodrich-Lab/epiomics")

library(epiomics)
```

## Omics wide association study (owas)

The basis of many omics analysis in epidemiology begin with an omics wide association study. The function `owas()` implements an omics wide association study with the option of using the 'omics data as either the dependent variable (i.e., for performing an exposure --\> 'omics analysis) or using the 'omics as the independent variable (i.e., for performing an 'omics --\> outcome analysis). `owas()` provides the option to adjust for covariates, and allows for either continuous or dichotomous outcomes. `owas()` can also handle multiple variables of interest (ie, multiple exposures or multiple traits).

Start with loading example data:

``` r
# Load Example Data
data("example_data")

# Get names of omics
colnames_omic_fts <- colnames(example_data)[grep("feature_",
                                               colnames(example_data))][1:10]

# Get names of traits
trait_nms = c("disease1", "disease2")
```

### Run owas with continuous exposure as the variable of interest

``` r
owas(df = example_data, 
     var = "exposure1", 
     omics = colnames_omic_fts, 
     covars = c("age", "sex"), 
     var_exposure_or_outcome = "exposure", 
     family = "gaussian")
     
# Equivalent: 
owas(df = example_data, 
     var = "exposure1", 
     omics = colnames_omic_fts, 
     covars = c("age", "sex"), 
     var_exposure_or_outcome = "exposure")  
```

### Run owas with dichotomous outcome as the variable of interest

``` r
owas(df = example_data, 
     var = "disease1", 
     omics = colnames_omic_fts, 
     covars = c("age", "sex"), 
     var_exposure_or_outcome = "outcome", 
     family = "binomial")
```

### Run owas with multiple continuous exposures as the variable of interest

``` r
# Get names of exposures
expnms = c("exposure1", "exposure2", "exposure3")

owas(df = example_data, 
     var = expnms, 
     omics = colnames_omic_fts, 
     covars = c("age", "sex"), 
     var_exposure_or_outcome = "exposure", 
     family = "gaussian")
```

## Meet in the Middle

The function `meet_in_middle()` conducts meet in the middle screening between an exposure, omics, and an outcome, as described by Cadiou et al., 2021. This function provides the option to adjust for covariates, and allows for either continuous or dichotomous outcomes. Examples are based on the simulated data created above.

### Meet in the middle with a dichotomous outcome

``` r
res <- meet_in_middle(df = example_data,
                      exposure = "exposure1", 
                      outcome = "disease1", 
                      omics = colnames_omic_fts,
                      covars = c("age", "sex"), 
                      outcome_family = "binomial")
res
```

### Meet in the middle with a continuous outcome

``` r
res <- meet_in_middle(df = example_data,
                      exposure = "exposure1", 
                      outcome = "weight", 
                      omics = colnames_omic_fts,
                      covars = c("age", "sex"), 
                      outcome_family = "gaussian")
```

### Meet in the middle with a continuous outcome and no covariates

``` r
res <- meet_in_middle(df = example_data,
                      exposure = "exposure1", 
                      outcome = "weight", 
                      omics = colnames_omic_fts,
                      outcome_family = "gaussian")
```


## Omics wide environmental chemical mixtures association study

The `owas_qgcomp()` function implements an omics wide association study using quantile-based g-Computation (as described by Keil et al., (2019) <doi:10.1289/EHP5838>) to examine associations of exposure mixtures with each individual 'omics feature as an outcome 'omics data as either the dependent variable. This function allows for either continuous or dichotomous outcomes, and provides the option to adjust for covariates.

### Run owas with continuous exposure as the variable of interest

``` r
exposure_names = c("exposure1", "exposure2", "exposure3")

# Run function without covariates
out <- owas_qgcomp(df = example_data,
                   expnms = exposure_names,
                   omics = colnames_omic_fts,
                   q = 4, 
                   confidence_level = 0.95) 
```
