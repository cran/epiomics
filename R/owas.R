#' Perform 'omics wide association study
#' @description 
#' Implements an omics wide association study with the option of using the 
#' 'omics data as either the dependent variable (i.e., for performing an 
#' exposure --> 'omics analysis) or using the 'omics as the independent
#'  variable (i.e., for performing an 'omics --> outcome analysis). Allows for
#'  either continuous or dichotomous outcomes, and provides the option to 
#'  adjust for covariates. 
#' 
#' @import data.table
#' @importFrom stats binomial coef glm lm p.adjust confint confint.default na.omit
#' @export 
#' @param df Dataset
#' @param var Name of the variable or variables of interest- this is usually
#'  either an exposure variable or an outcome variable. Can be either 
#'  continuous or dichotomous. For dichotomous variables, must set \code{family}
#'  to "binomial", and values must be either 0/1 or a factor with the first 
#'  level representing the reference group. Can handle multiple variables, but 
#'  they must all be of the same \code{family}.
#' @param omics Names of all omics features in the dataset 
#' @param covars Names of covariates (can be NULL)
#' @param var_exposure_or_outcome Is the variable of interest an exposure 
#' (independent variable) or outcome (dependent variable)? Must be either
#' "exposure" or "outcome"
#' @param family "gaussian" (default) for linear models (via lm) or "binomial"
#'  for logistic (via glm) 
#' @param confidence_level Confidence level for marginal significance 
#' (defaults to 0.95, or an alpha of 0.05)
#' @param conf_int Should Confidence intervals be generated for the estimates? 
#' Default is FALSE. Setting to TRUE will take longer. For logistic models, 
#' calculates Wald confidence intervals via \code{confint.default}.
#' @param ref_group Reference category if the variable of interest is a
#' character or factor. If not, can leave empty. 
#' 
#' @return
#' A data frame with 6 columns:  
#' feature_name: name of the omics feature   
#' estimate: the model estimate for the feature. For linear models, this is the 
#' beta; for logistic models, this is the log odds. 
#' se: Standard error of the estimate
#' test_statistic: t-value
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
#' # Get names of exposures
#' expnms = c("exposure1", "exposure2", "exposure3")
#' 
#' # Run function with one continuous exposure as the variable of interest
#' owas(df = example_data, 
#'      var = "exposure1", 
#'      omics = colnames_omic_fts, 
#'      covars = c("age", "sex"), 
#'      var_exposure_or_outcome = "exposure", 
#'      family = "gaussian")
#'      
#' # Run function with multiple continuous exposures as the variable of interest
#' owas(df = example_data, 
#'      var = expnms, 
#'      omics = colnames_omic_fts, 
#'      covars = c("age", "sex"), 
#'      var_exposure_or_outcome = "exposure", 
#'      family = "gaussian")
#' 
#' # Run function with dichotomous outcome as the variable of interest
#' owas(df = example_data, 
#'      var = "disease1", 
#'      omics = colnames_omic_fts, 
#'      covars = c("age", "sex"), 
#'      var_exposure_or_outcome = "outcome", 
#'      family = "binomial")
#' 
owas <- compiler::cmpfun(
  function(df, 
           var,
           omics, 
           covars = NULL,
           var_exposure_or_outcome, 
           family = "gaussian", 
           confidence_level = 0.95, 
           conf_int = FALSE, 
           ref_group = NULL){
    df <- base::as.data.frame(df)
    final_col_names <- ftr_var_group <- NULL
    alpha <- 1-confidence_level
    # Check for issues in data ----
    # Get var variable types
    var_types <- df[,(colnames(df) %in% var)] |>
      lapply(function(x)class(x)) |> unlist() |> unique()
    
    ## Check if variable of interest is in data ----
    if(FALSE %in% (var %in% colnames(df))){ 
      stop(paste0("Variable '", 
                  paste0(var[!(var %in% colnames(df))],
                         collapse = ", "),
                  "' not found in data. Check data.") ) 
    }    
    ## Check if var has different types ----
    if(length(var_types) > 1){ 
      stop("All variables in \'var\' must be the same type")  
    }   
    ## Check if var is numeric or character/factor with max 2 levels ----
    if((var_types == "character" | var_types == "factor") ){ 
      if(is.null(ref_group)){ 
        stop("If var is character or factor, ref_group must be specified")  
      } 
      if((df[,(colnames(df) %in% var)] |>
          as.matrix() |>
          na.omit() |>
          as.character() |> 
          unique() |>
          length()) > 2){ 
        stop("Var can only contain a maximum of two unique categories")  
      }
    }  
    ## Check if all omics features are in the data ----
    if(FALSE %in% (omics %in% colnames(df))){ 
      stop(
        "Not all omics vars are found in the data. Check omics column names.")  
    }    
    ## Check if covars are in data ----
    if(FALSE %in% (covars %in% colnames(df))){ 
      stop(
        "Not all covars are found in the data. Check covar column names."
      ) 
    }    
    ## Check that var is specified
    if(!(var_exposure_or_outcome %in% c("exposure", "outcome"))){ 
      stop(
        "var_exposure_or_outcome must be either \"exposure\" or \"outcome\" "
      )
    }
    ## Check that family is specified
    if(!(family %in% c("gaussian", "binomial"))){ 
      stop("family must be either \"gaussian\" or \"binomial\" ")
    }
    
    
    # Change data frame to data table for speed
    df <- data.table(df)
    
    # Pivot longer on omics 
    dt_l <- melt.data.table(
      df,
      id.vars = c(var, covars),
      measure.vars = omics,
      variable.name = "feature_name",
      value.name = "feature_value"
    )
    
    # Pivot longer on variables 
    dt_l2 <- melt.data.table(
      dt_l,
      id.vars = c("feature_name", "feature_value", covars),
      measure.vars = var,
      variable.name = "var_name",
      value.name = "var_value")
    
    # Create feature_name_var_name variable
    dt_l2$ftr_var_group <- paste0(dt_l2$feature_name, "_", dt_l2$var_name)
    
    # relevel var_value to specify correct reference group
    if(var_types == "character" | var_types == "factor") { 
      dt_l2$var_value <- ifelse(dt_l2$var_value == ref_group, 0, 1)
    }
    
    # Set formula for model ------------------
    # depending on whether variable of interest is the exposure or the outcome 
    if(var_exposure_or_outcome == "exposure"){
      # If variable is exposure: 
      if(is.null(covars)){
        # mod_formula <- paste0("feature_value~", var)
        mod_formula <- "feature_value~var_value"
      } else {
        mod_formula <- paste0("feature_value~", 
                              paste0(covars, collapse = "+"),
                              "+var_value")
      } 
      
    } else if(var_exposure_or_outcome == "outcome"){
      # If variable is outcome: 
      if(is.null(covars)){
        # mod_formula <- paste0(var, "~feature_value")
        mod_formula <- "var_value~feature_value"
      } else{
        mod_formula <- paste0("var_value~", 
                              paste0(covars, collapse = "+"),
                              "+feature_value")
      }
      
    } 
    # Run models -------------------------
    ## If no confidence intervals are requested: -----------------
    if(!conf_int){
      if(family == "gaussian"){
        # Linear models:
        res <- dt_l2[, 
                     {fit <- lm(mod_formula, data = .SD) 
                     coef(summary(fit))[nrow(coef(summary(fit))), 
                                        c(1, 2, 3, 4)]
                     }, 
                     by = ftr_var_group]
        
        # Add column for estimate 
        res <- cbind(res, c("estimate", "se", "test_statistic", "p_value"))
        
        # Pivot wider
        final_results <- dcast(data = res, 
                               ftr_var_group ~ V2, 
                               value.var = "V1")[,c(1, 2, 4, 5, 3)]
        
      } else if(family == "binomial"){
        # Logistic models
        res <- dt_l2[, 
                     {fit <- glm(mod_formula, 
                                 data = .SD, 
                                 family=binomial(link='logit')) 
                     coef(summary(fit))[nrow(coef(summary(fit))), 
                                        c(1, 2, 3, 4)] 
                     }, 
                     by = ftr_var_group]
        
        # Add column for estimate 
        res <- cbind(res, c("estimate", "se", "test_statistic", "p_value"))
        
        # Pivot wider
        final_results <- dcast(data = res, 
                               ftr_var_group ~ V2, 
                               value.var = "V1")[,c(1, 2, 4, 5, 3)]
        
      }
      
    } else if(conf_int){
      ## If confidence intervals are requested: --------------------
      if(family == "gaussian"){
        # Linear models:
        res <- dt_l2[, 
                     {fit <- lm(mod_formula, data = .SD) 
                     out <- coef(summary(fit))[nrow(coef(summary(fit))), 
                                               c(1, 2, 3, 4)] 
                     ci <- confint(fit)[length(coef(fit)), ]
                     c(out,ci)
                     }, 
                     by = ftr_var_group]
        
        # Add column for estimate 
        res <- cbind(res, c("estimate", "se", "test_statistic",
                            "p_value", "conf_low", "conf_high"))
        
        # Pivot wider
        final_results <- dcast(data = res, 
                               ftr_var_group ~ V2, 
                               value.var = "V1")[,c(1, 4, 6, 7, 5, 3, 2)]
        
      } else if(family == "binomial"){
        # Logistic models
        res <- dt_l2[, 
                     {fit <- glm(mod_formula, 
                                 data = .SD, 
                                 family=binomial(link='logit')) 
                     out <- coef(summary(fit))[nrow(coef(summary(fit))), 
                                               c(1, 2, 3, 4)] 
                     ci <- confint.default(fit)[length(coef(fit)), ]
                     c(out,ci)
                     }, 
                     by = ftr_var_group]
        
        # Add column for estimate 
        res <- cbind(res, c("estimate", "se", "test_statistic",
                            "p_value", "conf_low", "conf_high"))
        
        # Pivot wider
        final_results <- dcast(data = res, 
                               ftr_var_group ~ V2, 
                               value.var = "V1")[,c(1, 4, 6, 7, 5, 3, 2)]
        
      } 
    }
    
    # Separate ftr_var_group back to ftr name and var name ------------
    # get all combinations of ftnme and var name
    unique_ftnme_varnme_combo <- unique(dt_l2[,c("feature_name",
                                                 "var_name",
                                                 "ftr_var_group")])
    # merge with final data
    final_results_2 <- merge(final_results, 
                             unique_ftnme_varnme_combo, 
                             by = "ftr_var_group")
    
    # Calculate adjusted p value
    final_results_2$adjusted_pval <- p.adjust(final_results_2$p_value,
                                              method = "fdr")
    
    final_results_2$threshold <- ifelse(final_results_2$p_value < alpha,
                                        "Significant",
                                        "Non-significant")
    
    # Reorder 
    final_col_names <- c("var_name", "feature_name", 
                         setdiff(colnames(final_results_2), 
                                 c("ftr_var_group", 
                                   "feature_name", 
                                   "var_name")))
    
    # Select columns
    reordered <- final_results_2[,final_col_names, with = FALSE] 
    
    
    return(as.data.frame(reordered))
    
  }
)
