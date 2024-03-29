% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/owas_qgcomp.R
\name{owas_qgcomp}
\alias{owas_qgcomp}
\title{Perform omics wide association study using qgcomp}
\usage{
owas_qgcomp(
  df,
  expnms,
  omics,
  covars = NULL,
  q = 4,
  confidence_level = 0.95,
  family = "gaussian",
  rr = TRUE,
  run.qgcomp.boot = TRUE
)
}
\arguments{
\item{df}{Dataset}

\item{expnms}{Name of the exposures. Can be either continuous or
dichotomous. For dichotomous variables, must set \code{q} to
"NULL", and values must be either 0/1.}

\item{omics}{Names of all omics features in the dataset}

\item{covars}{Names of covariates (can be NULL)}

\item{q}{NULL or number of quantiles used to create quantile indicator
variables representing the exposure variables. Defaults to 4If NULL, then
qgcomp proceeds with un-transformed version of exposures in the input
datasets (useful if data are already transformed, or for performing standard
g-computation).}

\item{confidence_level}{Confidence level for marginal significance
(defaults to 0.95, or an alpha of 0.05)}

\item{family}{Currently only "gaussian" (default) for linear models (via lm)
or "binomial" for logistic. Default is "gaussian".}

\item{rr}{see \code{qgcomp()}}

\item{run.qgcomp.boot}{Should the model be fit with qgcomp.boot? See
package \link[qgcomp]{qgcomp.boot} for details. Default is TRUE.
Setting to FALSE decreases computational time.}
}
\value{
A data frame with the following columns:
feature: name of the omics feature
psi: the model estimate for the feature. For linear models, this is the
beta; for logistic models, this is the log odds.
lcl_psi: the lower confidence interval.
ucl_psi: the upper confidence interval.
p_value: p-value for the estimate
test_statistic: t-statistic for psi coefficient
adjusted_pval: FDR adjusted p-value
threshold: Marginal significance, based on unadjusted p-values
covariates: the names of covariates in the model, if any
coef_exposure: the individual coefficient of each exposure
}
\description{
Omics wide association study using quantile-based g-Computation
(as described by Keil et al., (2019) \doi{10.1289/EHP5838}) to examine
associations of exposure mixtures with each individual 'omics feature as an
outcome 'omics data as either the dependent variable. Allows for either
continuous or dichotomous outcomes, and provides the option to adjust for
covariates.
}
\examples{
# Load Example Data
data("example_data")

# Get names of omics
colnames_omic_fts <- colnames(example_data)[grep("feature_",
                                              colnames(example_data))][1:5]

# Names of exposures in mixture
 exposure_names = c("exposure1", "exposure2", "exposure3")

# Run function without covariates
out <- owas_qgcomp(df = example_data,
                   expnms = exposure_names,
                   omics = colnames_omic_fts,
                   q = 4, 
                   confidence_level = 0.95) 


# Run analysis with covariates
out <- owas_qgcomp(df = example_data,
                   expnms = c("exposure1", "exposure2", "exposure3"),
                   covars = c("weight", "age", "sex"),
                   omics = colnames_omic_fts,
                   q = 4, 
                   confidence_level = 0.95) 
 
}
