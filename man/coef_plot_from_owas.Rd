% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef_plot_from_owas.R
\name{coef_plot_from_owas}
\alias{coef_plot_from_owas}
\title{Create volcano plot using results from owas}
\usage{
coef_plot_from_owas(
  df,
  main_cat_var = NULL,
  order_effects = TRUE,
  highlight_adj_p = TRUE,
  highlight_adj_p_threshold = 0.05,
  effect_ratio = FALSE,
  flip_axis = FALSE,
  filter_p_less_than = 1
)
}
\arguments{
\item{df}{output from \code{owas} function call, using conf_int = TRUE.}

\item{main_cat_var}{Which variable should be the primary categorical
variable? Should be either var_name or feature_name. Only relevant if
both var_name and feature_name have more than one level. Default is NULL,
and the y-axis is chosen as the variable that has more levels.}

\item{order_effects}{Should features be ordered by the mean effect estimate?
Default is TRUE.}

\item{highlight_adj_p}{Should features which meet a specific adjusted p-value
threshold be highlighted? Default is TRUE.}

\item{highlight_adj_p_threshold}{If \code{highlight_adj_p} = TRUE, can set
annotation_adj_p_threshold to change the adjusted p-value threshold for
which features will be highlighted. Defaults to 0.05.}

\item{effect_ratio}{Are the effect estimates on the ratio scale (ie, should
the null effect line be centered at 1)? Defaults to FALSE.}

\item{flip_axis}{Flip the x and y axis? Default is FALSE, and the y-axis is
plotted with the features or variable names.}

\item{filter_p_less_than}{P-value threshold for which features/variables will
be included in the plot. Default is 1, and all features will be included.}
}
\value{
A ggplot figure
}
\description{
Creates a coefficient plot based on ggplot using the results from the
\code{owas} function.
}
\examples{
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
                 family = "gaussian", 
                 conf_int = TRUE)

coef_plot_from_owas(owas_out)

}
