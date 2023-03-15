#' Perform 'omics wide association study with linear or generalized mixed 
#' models
#' @description 
#' Implements an omics wide association study with the option of using the 
#' 'omics data as either the dependent variable (i.e., for performing an 
#' exposure --> 'omics analysis) or using the 'omics as the independent 
#' variable (i.e., for performing an 'omics --> outcome analysis). Allows for
#' either continuous or dichotomous outcomes, and provides the option to 
#' adjust for covariates. 
#' 
#' @import data.table 
#' @importFrom stats binomial coef glm lm p.adjust confint 
#' @importFrom lmerTest lmer 
#' @importFrom lme4 glmer confint.merMod
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
#' @param random_effects Random effects, formatted as specified by 
#' \link[lme4]{lmer} or \link[lme4]{glmer}
#' @param var_exposure_or_outcome Is the variable of interest an exposure 
#' (independent variable) or outcome (dependent variable)? Must be either
#' "exposure" or "outcome"
#' @param family "gaussian" (default) for linear models (via 
#' \link[lmerTest]{lmer}) or "binomial" for logistic (via \link[lme4]{glmer}) 
#' @param confidence_level Confidence level for marginal significance 
#' (defaults to 0.95, or an alpha of 0.05)
#' @param REML logical scalar - Should the estimates be chosen to optimize 
#' the REML criterion (as opposed to the log-likelihood)? Default is TRUE
#' @param conf_int Should Confidence intervals be generated for the estimates? 
#' Default is FALSE. Setting to TRUE will take longer. For logistic models, 
#' calculates Wald confidence intervals via \code{confint.default}.
#' @param ref_group Reference category if the variable of interest is a
#' character or factor. If not, can leave empty.
#' 
#' @returns 
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
owas_mixed_effects <- compiler::cmpfun(
  function(df, 
           var,
           omics, 
           random_effects, 
           covars = NULL,
           var_exposure_or_outcome, 
           family = "gaussian", 
           confidence_level = 0.95, 
           conf_int = FALSE, 
           REML = TRUE, 
           ref_group = NULL){
    df <- base::as.data.frame(df)
    final_col_names <- ftr_var_group <- NULL
    alpha <- 1-confidence_level
    
    # Check for issues in data 
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
          as.character() |> 
          unique() |>
          length()) > 2){ 
        stop("Var can only contain a maximum of two unique categories")  
      }
    }  
    # Check if all omics features are in the data
    if(FALSE %in% (omics %in% colnames(df))){ 
      stop("Not all omics vars are in the data. Check omics column names.")  
    }    
    # Check if covars are in data
    if(FALSE %in% (covars %in% colnames(df))){ 
      stop("Not all covariates are in the data. Check covariate column names.") 
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
    
    # Check if random effects are in data
    random_effect_vars_all <- gsub("\\|", ",", x = random_effects) |> 
      strsplit("\\,\\s|\\,|\\s") |> 
      unlist()
    
    random_effect_vars <- random_effect_vars_all[random_effect_vars_all %in% 
                                                   colnames(df)]
    if(length(random_effect_vars)==0){ 
      stop("Not all random effect vars are found in data")
    }
    
    # Change data frame to data table for speed ------------
    df <- data.table(df)
    
    # Pivot longer on omics 
    dt_l <- melt.data.table(
      df,
      id.vars = c(var, covars, random_effect_vars),
      measure.vars = omics,
      variable.name = "feature_name",
      value.name = "feature_value"
    )
    
    
    # Pivot longer on variables 
    dt_l2 <- melt.data.table(
      dt_l,
      id.vars = c("feature_name", "feature_value", covars, random_effect_vars),
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
        mod_formula <- paste0("feature_value~var_value+(",
                              random_effects, 
                              ")") |> 
          as.formula()
      } else {
        mod_formula <- paste0("feature_value~", 
                              paste0(covars, collapse = "+"),
                              "+var_value+(",
                              random_effects, ")") |>
          as.formula()
      } 
      
    } else if(var_exposure_or_outcome == "outcome"){
      # If variable is outcome: 
      if(is.null(covars)){
        mod_formula <- paste0("var_value~feature_value+(",
                              random_effects, 
                              ")") |> 
          as.formula()
      } else{
        mod_formula <- paste0("var_value~", 
                              paste0(covars, collapse = "+"),
                              "+feature_value+(",random_effects, ")") |> 
          as.formula()
      }
      
    } 
    
    # Run models -------------------------
    ## If no confidence intervals are requested: -----------------
    if(!conf_int){
      if(family == "gaussian"){
        # Linear models:
        res <- dt_l2[, 
                     {fit <- lmerTest::lmer(mod_formula, 
                                            data = .SD, 
                                            REML = REML) |>
                       summary() |> 
                       coef()
                     fit[nrow(fit), ]# Select last row of fixed effects
                     }, 
                     by = ftr_var_group]
        
        # Add column for estimate 
        res <- cbind(res, c("estimate",
                            "se", 
                            "df", 
                            "test_statistic",
                            "p_value"))
        
        # Pivot wider
        final_results <- dcast(data = res, 
                               ftr_var_group ~ V2, 
                               value.var = "V1")[,c(1,3,5,2,6,4)]
        
      } else if(family == "binomial"){
        # Logistic models
        res <- dt_l2[,
                     {fit <- glmer(mod_formula,
                                   data = .SD,
                                   family=binomial(link='logit')) |>
                       summary() |> 
                       coef()
                     fit[nrow(fit), ]# Select last row of fixed effects
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
                     {fit <- lmerTest::lmer(mod_formula, 
                                            data = .SD, 
                                            REML = REML)
                     coef_out <- fit |> summary() |> coef()
                     confint_out <- confint.merMod(fit, method = "Wald")
                     out <- coef_out[nrow(coef_out), ]# Select last row
                     
                     ci <- confint_out[nrow(confint_out), ]# Select last row
                     c(out,ci)
                     }, 
                     by = ftr_var_group]
        
        # Add column for estimate 
        res <- cbind(res, c("estimate", "se", "df", 
                            "test_statistic", "p_value", 
                            "conf_low", "conf_high"))
        
        # Pivot wider
        final_results <- dcast(data = res, 
                               ftr_var_group ~ V2, 
                               value.var = "V1")[,c(1,5,7,4,8,6,2,3)]
        
      } else if(family == "binomial"){
        # Logistic models
        res <- dt_l2[,
                     {fit <- lme4::glmer(mod_formula,
                                         data = .SD,
                                         family=binomial(link='logit'))
                     out_coef <- fit |> 
                       summary() |> 
                       coef()
                     out_ci <- confint.merMod(fit, method = "Wald")
                     
                     out <- out_coef[nrow(out_coef), ]# Select last row
                     
                     ci <- out_ci[nrow(out_ci), ]# Select last row
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
    reordered <- final_results_2[, final_col_names, with = FALSE] 
    
    
    return(as.data.frame(reordered))
    
  }
)
