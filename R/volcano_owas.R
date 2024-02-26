#' Create volcano plot using results from owas
#' @description 
#' Creates a volcano plot based on ggplot using the results from the 
#' \code{owas} function.  
#' @import ggplot2 ggrepel
#' @export 
#' @param df output from \code{owas} function call
#' @param annotate_ftrs Should features be annotated with the feature name? 
#' Default is TRUE. If necessary can change the p_value_threshold as well. 
#' @param annotation_p_threshold If \code{annotate_ftrs} = TRUE, can set 
#' annotation_p_threshold to change the p-value threshold for which features 
#' will be annotated. Defaults to 0.05. 
#' @param highlight_adj_p Should features which meet a specific adjusted p-value 
#' threshold be highlighted? Default is TRUE. 
#' @param highlight_adj_p_threshold If \code{highlight_adj_p} = TRUE, can set 
#' annotation_adj_p_threshold to change the adjusted p-value threshold for 
#' which features will be highlighted. Defaults to 0.05.
#' @param horizontal_line_p_value Set the p-value for the horizontal line 
#' for the threshold of significance. 
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
#'                  family = "gaussian")
#' 
#' vp <- volcano_owas(owas_out)
#' 
volcano_owas <- compiler::cmpfun(
  function(df, 
           annotate_ftrs = TRUE, 
           annotation_p_threshold = 0.05, 
           highlight_adj_p = TRUE, 
           highlight_adj_p_threshold = 0.05, 
           horizontal_line_p_value = 0.05){

    estimate <- p_value <- ftr_names_label <- adjusted_pval <- NULL
    # Annotated feature names
    df$ftr_names_label <- ifelse(df$p_value < annotation_p_threshold, 
                                 as.character(df$feature_name), "")
    
    # Make plot
    main_plot <- ggplot(df, 
                        aes(x = estimate,
                            y = -log10(p_value))) + 
      geom_hline(yintercept = -log10(horizontal_line_p_value), 
                 linetype = 2) + 
      ylab("-log10 p-value") 
    
    # Feature Annotation
    if(annotate_ftrs){
      main_plot <- main_plot + 
        geom_text_repel(aes(label = ftr_names_label))
    }
    
    # Feature Annotation
    if(highlight_adj_p){
      main_plot <- main_plot + 
        geom_point(aes(color = adjusted_pval<highlight_adj_p_threshold)) + 
        scale_color_manual(values = c("black", "red"), 
                           labels=c(paste0("Adjusted p-value > ", 
                                           highlight_adj_p_threshold), 
                                    paste0("Adjusted p-value < ", 
                                           highlight_adj_p_threshold)), 
                           name = "Adjusted p-value")
    } else {
      main_plot <- main_plot + 
        geom_point() 
    }
    
    # Facet wrap if multiple exposures
    if(length(unique(df$var_name))>1){
      main_plot <- main_plot + 
        facet_wrap(~var_name)
    }
  # return plot
  main_plot + theme_classic()
  }
)
