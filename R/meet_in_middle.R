#' Perform 'omics wide association study
#' @description 
#' Implements a meet in the middle analysis for identifying omics associated 
#' with both exposures and outcomes, as described by Chadeau-Hyam et al., 2010.
#' 
#' @import data.table
#' @export 
#' @param df Dataframe
#' @param exposure Name of the exposure of interest. Can be either continuous 
#' or dichotomous. Currently, only a single exposure is supported.
#' @param outcome Name of the outcome of interest. Can be either continuous or 
#' dichotomous. For dichotomous variables, must set \code{outcome_family} to
#'  "logistic", and values must be either 0/1 or a factor with the first level 
#'  representing the reference group. Currently, only a single outcome is 
#'  supported.
#' @param omics Names of all omics features in the dataset 
#' @param covars Names of covariates (can be NULL)
#' @param outcome_family "gaussian" for linear models (via lm) or "binomial" 
#' for logistic (via glm) 
#' @param confidence_level Confidence level for marginal significance 
#' (defaults to 0.95)
#' @param conf_int Should Confidence intervals be generated for the estimates? 
#' Default is FALSE. Setting to TRUE will take longer. For logistic models, 
#' calculates Wald confidence intervals via \code{confint.default}.
#' @param ref_group_exposure Reference category if the exposure is a
#' character or factor. If not, can leave empty. 
#' @param ref_group_outcome Reference category if the outcome is a
#' character or factor.  If not, can leave empty. 
#' 
#' @returns 
#' A list of three dataframes, containing: 
#' 1) Results from the Exposure-Omics Wide Association Study  
#' 2) Results from the Omics-Outcome Wide Association Study
#' 3) Overlapping significant features from 1 and 2.     
#' For each omics wide association, results are provided in a data frame with 6 
#' columns:  
#' feature_name: name of the omics feature   
#' estimate: the model estimate for the feature. For linear models, this is the 
#' beta: for logistic models, this is the log odds. 
#' se: Standard error of the estimate
#' p_value: p-value for the estimate
#' adjusted_pval: FDR adjusted p-value
#' threshold: Marginal significance, based on unadjusted p-values 
#' 
#' @examples 
#' # Load Example Data
#' data("example_data")
#' 
#' # Get names of omics
#' colnames_omic_fts <- colnames(example_data)[grep("feature_",
#'                                               colnames(example_data))][1:10]
#' 
#' # Meet in the middle with a dichotomous outcome
#' res <- meet_in_middle(df = example_data,
#'                       exposure = "exposure1", 
#'                       outcome = "disease1", 
#'                       omics = colnames_omic_fts,
#'                       covars = c("age", "sex"), 
#'                       outcome_family = "binomial")
#' 
#' # Meet in the middle with a continuous outcome 
#' res <- meet_in_middle(df = example_data,
#'                       exposure = "exposure1", 
#'                       outcome = "weight", 
#'                       omics = colnames_omic_fts,
#'                       covars = c("age", "sex"), 
#'                       outcome_family = "gaussian")
#' 
#' # Meet in the middle with a continuous outcome and no covariates
#' res <- meet_in_middle(df = example_data,
#'                       exposure = "exposure1", 
#'                       outcome = "weight", 
#'                       omics = colnames_omic_fts,
#'                       outcome_family = "gaussian")
#' 
meet_in_middle <- compiler::cmpfun(
  function(df, 
           exposure,
           outcome,
           omics,
           covars = NULL,
           outcome_family = "gaussian",
           confidence_level = 0.95, 
           conf_int = FALSE, 
           ref_group_exposure = NULL, 
           ref_group_outcome = NULL){
    alpha <- 1-confidence_level
    
    # Check if more than one exposure
    if((length(exposure)+length(outcome))>2){ 
      stop("More than one exposure or outcome is not currently supported.") 
    }   
    
    
    df <- data.table::as.data.table(df)
    # exposure_omics_owas
    exposure_omics_owas <- epiomics::owas(df = df,
                                          var = exposure, 
                                          omics = omics, 
                                          covars = covars,
                                          var_exposure_or_outcome = "exposure", 
                                          family = "gaussian",
                                          confidence_level = confidence_level, 
                                          conf_int = conf_int, 
                                          ref_group = ref_group_exposure)
    
    # omics_outcome_owas
    omics_outcome_owas <- epiomics::owas(df = df,
                                         var = outcome, 
                                         omics = omics,
                                         covars = covars,
                                         var_exposure_or_outcome = "outcome", 
                                         family = outcome_family,
                                         confidence_level = confidence_level, 
                                         conf_int = conf_int, 
                                         ref_group = ref_group_outcome)
    
    
    # Find overlap
    x_o_fts <- exposure_omics_owas[
      exposure_omics_owas$p_value<alpha,]$feature_name
    o_y_fts <- omics_outcome_owas[
      omics_outcome_owas$p_value<alpha,]$feature_name
    overlap_fts <- intersect(x_o_fts, o_y_fts)
    
    # Create overlapping data frame
    if(length(overlap_fts)>0){
      # Subset only overlapping sig
      # Exposure-feature
      x_o <- exposure_omics_owas[
        exposure_omics_owas$feature_name %in% overlap_fts,]
      colnames(x_o) <- c("exposure_name", 
                         paste0(colnames(x_o)[-1],
                                c("",rep("_exp_omic", ncol(x_o)-2))))
      # feature-outcome
      o_y <- omics_outcome_owas[
        omics_outcome_owas$feature_name %in% overlap_fts,]
      colnames(o_y) <- c("outcome_name", 
                         paste0(colnames(o_y)[-1],
                                c("",rep("_omic_out", ncol(o_y)-2))))
      
      # Merge
      overlap <- merge(x_o, o_y, by = c("feature_name")) 
      # Reorder
      first_cols <- c("exposure_name", "outcome_name", "feature_name")
      overlap <- overlap[,c(first_cols, 
                            setdiff(colnames(overlap), first_cols))]
      
      
    } else {overlap <- NULL}
    
    final_results <- list(exposure_omics_owas = exposure_omics_owas, 
                          omics_outcome_owas = omics_outcome_owas, 
                          overlap = overlap)
    
    return(final_results)
    
  }
)
