#' Perform 'omics wide association study for matched case control studies
#' @description 
#' Implements an omics wide association study for matched case control studies 
#' using conditional logistic regression. For this function, the variable of 
#' of interest should be a dichotomous outcome, and the strata is the variable
#' indicating the matching. 
#' 
#' @importFrom stats coef glm lm p.adjust confint confint.default as.formula var
#' @importFrom survival clogit strata coxph Surv
#' @export 
#' @param df Dataset
#' @param cc_status Name of the variable indicating case control status.
#'  Must be either 0/1 or a factor with the first level representing the
#'  reference group. 
#' @param cc_set Name of the variable indicating the case control set.
#' @param omics Names of all omics features in the dataset 
#'  reference group. 
#' @param covars Names of covariates (can be NULL)
#' @param confidence_level Confidence level for marginal significance 
#' (defaults to 0.95, or an alpha of 0.05)
#' @param conf_int Should Confidence intervals be generated for the estimates? 
#' Default is FALSE. Setting to TRUE will take longer. For logistic models, 
#' calculates Wald confidence intervals via \code{confint.default}.
#' @param method method used the correct (exact) calculation in the
#'  conditional likelihood or one of the approximations. Default is "efron". 
#'  Passed to \code{clogit}.
#' @param test_data_quality If TRUE (default), then code will ensure that 
#' the variance of all variables in the analysis is greater than 0 after 
#' dropping any missing data.
#' 
#' @return
#' A data frame with 6 columns:  
#' feature_name: name of the omics feature   
#' estimate: the model estimate for the feature. For linear models, this is the 
#' beta; for logistic models, this is the log odds. 
#' se: Standard error of the estimate
#' test statistic: t-value
#' p_value: p-value for the estimate
#' adjusted_pval: FDR adjusted p-value
#' threshold: Marginal significance, based on unadjusted p-values 
#' 
owas_clogit <- function(df, 
                        cc_status,
                        cc_set, 
                        omics, 
                        covars = NULL,
                        confidence_level = 0.95, 
                        conf_int = FALSE, 
                        method = "efron", 
                        test_data_quality = TRUE){
  
  alpha <- 1-confidence_level
  
  # Check for issues in data 
  # Check if variable of interest is in data
  if(FALSE %in% (cc_status %in% colnames(df))){ 
    stop(paste0("Variable '", 
                paste0(cc_status[!(cc_status %in% colnames(df))],
                       collapse = ", "),
                "' not found in data. Check data.") ) 
  }    
  # Check if all omics features are in the data
  if(FALSE %in% (omics %in% colnames(df))){ 
    stop(
      "Not all omics vars are found in the data. Check omics column names."
    ) 
  }    
  # Check if covars are in data
  if(FALSE %in% (covars %in% colnames(df))){ 
    stop(
      "Not all covars are found in the data. Check covar column names."
    )  
  }    
  
  # Test whether any omics have a variance of zero -------
  if(test_data_quality == TRUE){
    # Change data frame to data table for speed
    df2 <- data.table(df)
    # Drop rows with any missing values in the selected columns
    df_complete <- df2[complete.cases(df2[, c(omics, covars), with = FALSE]), 
                       c(omics, covars), with = FALSE]
    
    # For any non-numeric variables, turn them into factors then numeric
    var_unique_summary <- df_complete[, lapply(.SD, function(x) {
      if (is.numeric(x)) var(x, na.rm = TRUE) else (length(unique(x))-1)
    }), .SDcols = c(omics, covars)]
    
    var_unique_summary <- as.matrix(var_unique_summary==0)
    # Check which variables have exactly zero variance
    zero_var_vars <- colnames(var_unique_summary)[var_unique_summary != 0]
    
    if(length(zero_var_vars) > 0){
      stop(paste0("The following variables have zero variance for complete cases: ", 
                  paste(zero_var_vars, collapse = ", "), 
                  ". Please remove before analysis."))
    }
  }
  
  
  # Set formula for model ------------------
  # depending on whether covariates are included 
  if(is.null(covars)){
    # Initialize data frame
    mod_formula_df <- data.frame(
      cc_status = rep(cc_status, length(omics)), 
      feature_name = omics, 
      cc_set = rep(cc_set, length(omics)))
    # Set formula
    mod_formula_df$formula <- paste0(mod_formula_df$cc_status, 
                                     "~", 
                                     mod_formula_df$feature_name, 
                                     "+strata(",
                                     mod_formula_df$cc_set,
                                     ")")
  } else{  
    # Initialize data frame
    mod_formula_df <- data.frame(
      cc_status = rep(cc_status, length(omics)), 
      feature_name = omics, 
      cc_set = rep(cc_set, length(omics)), 
      covars = rep(paste0(covars, collapse = "+")))
    # Set formula
    mod_formula_df$formula <- paste0(mod_formula_df$cc_status, 
                                     "~", 
                                     mod_formula_df$covars, 
                                     "+",
                                     mod_formula_df$feature_name, 
                                     "+strata(",
                                     mod_formula_df$cc_set,
                                     ")")
  }
  
  
  # Run models -------------------------
  ## If no confidence intervals are requested: -----------------
  if(!conf_int){
    # Get empty data frame
    res_out <- data.frame(feature_name = omics, 
                          estimate = vector("numeric", length(omics)), 
                          se_est = vector("numeric", length(omics)), 
                          test_statistic = vector("numeric", length(omics)), 
                          p_value = vector("numeric", length(omics)), 
                          formula = mod_formula_df$formula) 
    
    # Run for loop to get results
    for(i in seq_along(omics)){ 
      mod_formula <- as.formula(res_out$formula[i])
      fit <- clogit(mod_formula, data = df, method = method) |>
        summary() |> 
        coef()
      res_out[i,2:5] <- fit[nrow(fit),  # Select last row
                            c(1, 3, 4, 5)]
    }
    
  } else if(conf_int){
    ## If confidence intervals are requested: --------------------
    # Get empty data frame
    res_out <- data.frame(feature_name = omics, 
                          estimate = vector("numeric", length(omics)), 
                          se_est = vector("numeric", length(omics)), 
                          test_statistic = vector("numeric", length(omics)), 
                          p_value = vector("numeric", length(omics)), 
                          conf_low = vector("numeric", length(omics)), 
                          conf_high = vector("numeric", length(omics)),
                          formula = mod_formula_df$formula) 
    
    # Run for loop to get results
    for(i in seq_along(omics)){ 
      mod_formula <- as.formula(res_out$formula[i])
      fit <- clogit(mod_formula, data = df, method = method)
      mod_coef <- fit |> summary() |> coef()
      mod_ci <- fit |> confint(level = confidence_level)
      
      res_out[i,2:5] <- mod_coef[nrow(mod_coef),  # Select last row
                                 c(1, 3, 4, 5)]
      res_out[i, 6:7] <- mod_ci[nrow(mod_ci), ]  # Select last row
    }
  }
  
  # Separate ftr_var_group back to ftr name and var name ------------
  # Calculate adjusted p value
  res_out$adjusted_pval <- p.adjust(res_out$p_value,
                                    method = "fdr")
  
  res_out$threshold <- ifelse(res_out$p_value < alpha,
                              "Significant",
                              "Non-significant")
  
  
  # Reorder 
  final_col_names <- c(setdiff(colnames(res_out), 
                               c("threshold", 
                                 "formula")), 
                       "threshold", "formula")
  
  # Select columns
  final_out <- res_out[,final_col_names] 
  
  return(final_out)
  
}
