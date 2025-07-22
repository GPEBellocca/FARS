#' @title Plot Loadings from \code{mldfm} Object
#'
#' @description Displays bar plots of the estimated factor loadings with 95% confidence intervals.
#'
#' @param object An object of class \code{mldfm}.
#' @param var_names Optional vector of variable names. If NULL, default names are used.
#' @param ... Additional arguments (ignored).
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate filter
#' @importFrom ggplot2 ggplot geom_bar geom_errorbar geom_hline theme_bw theme element_blank element_text scale_y_continuous ggtitle coord_flip
#' @importFrom forcats fct_rev
#'
#' @keywords internal
plot_loadings.mldfm <- function(object, var_names = NULL, ...) {
  stopifnot(inherits(object, "mldfm"))
  
  factors   <- get_factors(object)
  loadings  <- get_loadings(object)
  residuals <- get_residuals(object)
  
  t <- nrow(residuals)
  N <- ncol(residuals)
  
  loadings_df <- as.data.frame(loadings)
  
  # Factor names
  keys   <- names(object$factors_list)
  values <- unlist(object$factors_list)
  factor_names <- unlist(
    mapply(function(key, val) {
      clean <- paste0("F", gsub("-", "", key))
      if (val > 1) paste0(clean, "n", seq_len(val)) else clean
    }, keys, values, SIMPLIFY = FALSE)
  )
  colnames(loadings_df) <- factor_names
  
  # Variable names
  loadings_df$Variables <- if (is.null(var_names)) {
    paste0("Var", seq_len(nrow(loadings_df)))
  } else {
    var_names
  }
  
  loadings_df <- loadings_df[, c("Variables", setdiff(names(loadings_df), "Variables"))]
  
  # Long format
  loadings_long <- loadings_df %>%
    pivot_longer(cols = -Variables, names_to = "Factor", values_to = "Loading")
  
  # Standard errors and CIs
  se_vector <- apply(residuals, 2, sd) / sqrt(t)
  
  loadings_long <- loadings_long %>%
    mutate(SE = rep(se_vector, times = length(unique(Factor))),
           Loading_lower = Loading - 1.96 * SE,
           Loading_upper = Loading + 1.96 * SE)
  
  # Plot
  unique_factors <- unique(loadings_long$Factor)
  y_min <- -1
  y_max <- 1
  
  for (factor_name in unique_factors) {
    df_i <- loadings_long %>%
      filter(Factor == factor_name & Loading != 0) %>%
      mutate(Variables = factor(Variables, levels = unique(Variables)))
    
    p <- ggplot(df_i, aes(x = fct_rev(Variables), y = Loading)) +
      geom_bar(stat = "identity", fill = "grey", alpha = 0.7) +
      geom_hline(yintercept = 0, color = "red") +
      geom_errorbar(aes(ymin = Loading_lower, ymax = Loading_upper),
                    width = 0.5, color = "black", alpha = 1, size = 0.2) +
      coord_flip() +
      theme_bw() +
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5)
      ) +
      scale_y_continuous(limits = c(y_min, y_max)) +
      ggtitle(factor_name)
    
    print(p)
  }
  
  invisible(NULL)
}
