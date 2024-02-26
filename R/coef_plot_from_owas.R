#' Create volcano plot using results from owas
#' @description 
#' Creates a coefficient plot based on ggplot using the results from the 
#' \code{owas} function.  
#' @import ggplot2 ggrepel
#' @importFrom stats reorder
#' @export 
#' @param df output from \code{owas} function call, using conf_int = TRUE. 
#' @param main_cat_var Which variable should be the primary categorical 
#' variable? Should be either var_name or feature_name. Only relevant if 
#' both var_name and feature_name have more than one level. Default is NULL,
#' and the y-axis is chosen as the variable that has more levels. 
#' @param order_effects Should features be ordered by the mean effect estimate? 
#' Default is TRUE.
#' @param highlight_adj_p Should features which meet a specific adjusted p-value 
#' threshold be highlighted? Default is TRUE. 
#' @param highlight_adj_p_threshold If \code{highlight_adj_p} = TRUE, can set 
#' annotation_adj_p_threshold to change the adjusted p-value threshold for 
#' which features will be highlighted. Defaults to 0.05.
#' @param effect_ratio Are the effect estimates on the ratio scale (ie, should
#' the null effect line be centered at 1)? Defaults to FALSE.
#' @param flip_axis Flip the x and y axis? Default is FALSE, and the y-axis is 
#' plotted with the features or variable names.  
#' @param filter_p_less_than P-value threshold for which features/variables will 
#' be included in the plot. Default is 1, and all features will be included.
#' @return
#' A ggplot figure
#' 
#' @examples 
#' data("example_data")
#' 
#' # Get names of omics
#' colnames_omic_fts <- colnames(example_data)[
#'   grep("feature_",
#'        colnames(example_data))][1:5]
#' 
#' # Run function with continuous exposure as the variable of interest
#' owas_out <- owas(df = example_data,
#'                  var = "exposure1",
#'                  omics = colnames_omic_fts,
#'                  covars = c("age", "sex"),
#'                  var_exposure_or_outcome = "exposure",
#'                  family = "gaussian", 
#'                  conf_int = TRUE)
#' 
#' coef_plot_from_owas(owas_out)
#' 
coef_plot_from_owas <- compiler::cmpfun(
  function(df, 
           main_cat_var = NULL, 
           order_effects = TRUE, 
           highlight_adj_p = TRUE, 
           highlight_adj_p_threshold = 0.05, 
           effect_ratio = FALSE, 
           flip_axis = FALSE, 
           filter_p_less_than = 1){
    
    estimate <- p_value <- adjusted_pval <- facet_var <- 
      conf_high <- conf_low <- fdr_sig <- NULL
    
    # Determine number of unique variables and unique features 
    n_var <- length(unique(df$var_name))
    n_ftr <- length(unique(df$feature_name))
    
    # Error if n_var and n_ftr are both > 1 and main_cat_var isn't specified
    if(n_var > 1 & n_ftr > 1 & is.null(main_cat_var)){
      stop("main_cat_var must be specified if n_var and n_ftr are both > 1")
    }
    
    # Error if conf intervals not found
    if(!all(c("conf_low", "conf_high") %in% colnames(df))){
      stop("Confidence intervals named conf_low and conf_high must be in data frame. 
           Hint: In owas function, check to make sure conf_int is set to TRUE.")
    }
    # Give error if y axis not found 
    if(!is.null(main_cat_var)){
      if(!(main_cat_var %in% c("var_name", "feature_name", "estimate"))){
        stop("Y axis column name not found in column names")
      }
    }
    
    # Give error if highlight_adj_p_threshold is not between 
    # 0<highlight_adj_p_threshold<=1
    if(highlight_adj_p_threshold > 1 | highlight_adj_p_threshold <= 0){
      stop("highlight_adj_p_threshold must be between 0 and 1")
    }
    
    
    # Set main_cat_var if not specified in function call
    if(is.null(main_cat_var)){
      if(n_var == 1 & n_ftr > 1){
        main_cat_var <- "feature_name"
      } else if(n_var > 1 & n_ftr == 1){
        main_cat_var <- "var_name"
      } else if(n_var == 1 & n_ftr == 1){
        main_cat_var <- "feature_name"
      } 
    }
    
    # Set facet_var, if needed
    if(n_var > 1 & n_ftr > 1){
      facet_var <- setdiff(c("feature_name", "var_name"), main_cat_var)
    }  
    
    # Filter out values based on p-value threshold
    df <- df[df$p_value < filter_p_less_than, ]
    
    # Set main_cat_var
    df$main_cat_var <- df[[main_cat_var]]
    # Set facet_var 
    if(!is.null(facet_var)) {
      df$facet_var <- df[[facet_var]]
    }
    
    # Set x line annotation at 0 or 1
    if(effect_ratio){x_int = 1} else {x_int = 0}
    
    # reaorder main_cat_var by the mean of the effect estimates per facet_var
    if(order_effects){
      df$main_cat_var <- reorder(df$main_cat_var, df$estimate, FUN = mean)
    }
    

    # Color by fdr sig
    df$fdr_sig <- df$adjusted_pval < highlight_adj_p_threshold
    
    # Initialize plot ----
    main_plot <- ggplot(df, 
                         aes(x = estimate,
                             y = main_cat_var)) + 
       geom_vline(xintercept = x_int, linetype = 2) +
       xlab("Estimate (95% CI)") + 
       ylab(NULL) + 
      theme_classic()
    
    
    # Color by significance
    if(highlight_adj_p){
      main_plot <- main_plot + 
        geom_point(aes(color = fdr_sig)) + 
        geom_errorbar(aes(xmin = conf_low, 
                          xmax = conf_high, 
                          color = fdr_sig), 
                      width = 0) +
        scale_color_manual(values = c("black", "red"), 
                           labels=c(paste0("Adjusted p-value > ", 
                                           highlight_adj_p_threshold), 
                                    paste0("Adjusted p-value < ", 
                                           highlight_adj_p_threshold)), 
                           name = "Adjusted p-value")
    } else {
      main_plot <- main_plot + 
        geom_point() + 
        geom_errorbar(aes(xmin = conf_low, 
                          xmax = conf_high), 
                      width = 0) 
    }
    
    # Flip axis if needed
    if(flip_axis){
      main_plot <- main_plot + 
        coord_flip() + 
        theme(axis.text.x = element_text(angle = 90))
    }
    
    # Facet wrap if multiple exposures
    if(!is.null(facet_var) & flip_axis){
      main_plot <- main_plot + 
         facet_grid(facet_var~.)
    } else if(!is.null(facet_var)) {
      main_plot <- main_plot + 
        facet_grid(~facet_var)
    }
    
    # return plot
    return(main_plot) 
  }
)
