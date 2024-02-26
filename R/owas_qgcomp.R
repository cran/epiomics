#' Perform omics wide association study using qgcomp
#' @description 
#' Omics wide association study using quantile-based g-Computation 
#' (as described by Keil et al., (2019) \doi{10.1289/EHP5838}) to examine
#' associations of exposure mixtures with each individual 'omics feature as an
#' outcome 'omics data as either the dependent variable. Allows for either 
#' continuous or dichotomous outcomes, and provides the option to adjust for 
#' covariates. 
#' @import data.table qgcomp
#' @importFrom stats binomial gaussian coef glm lm p.adjust  
#' @export 
#' @param df Dataset
#' @param expnms Name of the exposures. Can be either continuous or 
#' dichotomous. For dichotomous variables, must set \code{q} to
#'  "NULL", and values must be either 0/1.
#' @param omics Names of all omics features in the dataset 
#' @param covars Names of covariates (can be NULL)
#' @param q NULL or number of quantiles used to create quantile indicator
#' variables representing the exposure variables. Defaults to 4If NULL, then
#' qgcomp proceeds with un-transformed version of exposures in the input 
#' datasets (useful if data are already transformed, or for performing standard
#' g-computation). 
#' @param confidence_level Confidence level for marginal significance 
#' (defaults to 0.95, or an alpha of 0.05)
#' @param family Currently only "gaussian" (default) for linear models (via lm) 
#' or "binomial" for logistic. Default is "gaussian".
#' @param rr see \code{qgcomp()}
#' @param run.qgcomp.boot Should the model be fit with qgcomp.boot? See 
#' package \link[qgcomp]{qgcomp.boot} for details. Default is TRUE. 
#' Setting to FALSE decreases computational time.  
#' 
#' @return
#' A data frame with the following columns:  
#' feature: name of the omics feature   
#' psi: the model estimate for the feature. For linear models, this is the 
#' beta; for logistic models, this is the log odds. 
#' lcl_psi: the lower confidence interval. 
#' ucl_psi: the upper confidence interval. 
#' p_value: p-value for the estimate
#' test_statistic: t-statistic for psi coefficient
#' adjusted_pval: FDR adjusted p-value
#' threshold: Marginal significance, based on unadjusted p-values 
#' covariates: the names of covariates in the model, if any
#' coef_exposure: the individual coefficient of each exposure
#' 
#' @examples 
#' # Load Example Data
#' data("example_data")
#' 
#' # Get names of omics
#' colnames_omic_fts <- colnames(example_data)[grep("feature_",
#'                                               colnames(example_data))][1:5]
#' 
#' # Names of exposures in mixture
#'  exposure_names = c("exposure1", "exposure2", "exposure3")
#' 
#' # Run function without covariates
#' out <- owas_qgcomp(df = example_data,
#'                    expnms = exposure_names,
#'                    omics = colnames_omic_fts,
#'                    q = 4, 
#'                    confidence_level = 0.95) 
#' 
#' 
#' # Run analysis with covariates
#' out <- owas_qgcomp(df = example_data,
#'                    expnms = c("exposure1", "exposure2", "exposure3"),
#'                    covars = c("weight", "age", "sex"),
#'                    omics = colnames_omic_fts,
#'                    q = 4, 
#'                    confidence_level = 0.95) 
#'  
owas_qgcomp <- compiler::cmpfun(
  function(df, 
           expnms,
           omics, 
           covars = NULL,
           q = 4, 
           confidence_level = 0.95, 
           family = "gaussian", 
           rr = TRUE, 
           run.qgcomp.boot = TRUE){
    
    alpha <- 1-confidence_level
    # Get var variable types
    var_types <- df[,(colnames(df) %in% expnms)] |>
      lapply(function(x)class(x)) |> unlist() |> unique()
    
    ## Check if variable of interest is in data ----
    if(FALSE %in% (expnms %in% colnames(df))){ 
      stop(paste0("Variable '", 
                  paste0(expnms[!(expnms %in% colnames(df))],
                         collapse = ", "),
                  "' not found in data. Check data.") ) 
    }    
    ## Check if var has different types ----
    if(length(var_types) > 1){ 
      stop("All variables in \'expnms\' must be the same type")  
    }   
    ## Check if var is numeric or character/factor with max 2 levels ----
    if((var_types == "character" | var_types == "factor")){ 
      stop("Currently exposures must be numeric, consider reformatting")  
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
    ## Check that family is specified ----
    if(!(family %in% c("gaussian", "binomial"))){ 
      stop("family must be either \"gaussian\" or \"binomial\" ")
    }
    
    # Set family
    if(family == "gaussian"){family_type = gaussian()}
    if(family == "binomial"){family_type = binomial()}
    
    
    # Change data frame to data table for speed
    df <- data.table(df)
    
    # Pivot longer  
    dt_l <- melt.data.table(
      df,
      id.vars = c(expnms, covars),
      measure.vars = omics,
      variable.name = "feature_name",
      value.name = "feature_value"
    )
    
    # Run models -------------------------
    # Split to list
    dt_list <- split(dt_l, f = dt_l$feature_name)
    
    # Run qgcomp noboot
    if(run.qgcomp.boot){
      res <- lapply(dt_list, 
                    FUN = function(x){ 
                      dt_l <- x[,c(colnames(x) %in% 
                                     c(expnms, "feature_value", covars)),
                                with = FALSE]
                      fit <- qgcomp::qgcomp(feature_value~.,
                                            expnms = expnms, 
                                            q = q,
                                            alpha = 0.05,
                                            data = dt_l, 
                                            family = family_type, 
                                            rr = rr)
                      # Get names of coefficients
                      names(fit$fit$coefficients) <-  paste0(
                        "coef_", 
                        names(fit$fit$coefficients))
                      # Get statistic names
                      stat_name <- names(fit)[grep("stat", names(fit))]
                      c(psi            = as.numeric(fit$coef[2]),
                        lcl_psi        = as.numeric(fit$ci[1]),
                        ucl_psi        = as.numeric(fit$ci[2]),
                        test_statistic = as.numeric(fit[[stat_name]][2]),
                        p_value        = as.numeric(fit$pval[2]), 
                        fit$fit$coefficients[-1])
                      
                    })
    } else {
      res <- lapply(dt_list, 
                    FUN = function(x){ 
                      dt_l <- x[,c(colnames(x) %in% 
                                     c(expnms, "feature_value", covars)),
                                with = FALSE]
                      fit <- qgcomp::qgcomp.noboot(feature_value~.,
                                            expnms = expnms, 
                                            q = q,
                                            alpha = 0.05,
                                            data = dt_l, 
                                            family = family_type)
                      # Get names of coefficients
                      names(fit$fit$coefficients) <-  paste0(
                        "coef_", 
                        names(fit$fit$coefficients))
                      # Get statistic names
                      stat_name <- names(fit)[grep("stat", names(fit))]
                      c(psi            = as.numeric(fit$coef[2]),
                        lcl_psi        = as.numeric(fit$ci[1]),
                        ucl_psi        = as.numeric(fit$ci[2]),
                        test_statistic = as.numeric(fit[[stat_name]][2]),
                        p_value        = as.numeric(fit$pval[2]), 
                        fit$fit$coefficients[-1])
                    })
      
    }
    
    # Add column for estimate
    res_2 <- t(as.data.frame(res))
    res_2 <- as.data.frame(res_2)
    res_2$feature <- rownames(res_2)
    rownames(res_2) <- NULL
    
    # Calculate adjusted p value
    res_2$adjusted_pval <- p.adjust(res_2$p_value, method = "fdr")
    
    res_2$threshold <- ifelse(res_2$adjusted_pval < alpha,
                              "Significant",
                              "Non-significant")
    
    # Reorder (select feature then all other cols)
    final_results <- res_2[
      c("feature", "psi", 
        "lcl_psi", "ucl_psi", 
        "p_value", "test_statistic",
        "adjusted_pval", "threshold", 
        paste0("coef_", expnms))
    ]
    
    # Add column for covariates
    if(is.null(covars)){
      final_results$covariates <- "None"
    } else {
      final_results$covariates <- paste0(covars, collapse = ", ")
    }
    
    return(final_results)
    
  }
  
)
