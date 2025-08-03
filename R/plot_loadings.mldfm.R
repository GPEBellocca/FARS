#' @title Plot Loadings from \code{mldfm} Object
#'
#' @description Displays bar plots of the estimated factor loadings with 95% confidence intervals.
#'
#' @param x An object of class \code{mldfm}.
#' @param var_names Optional vector of variable names. If NULL, default names are used.
#' @param ... Additional arguments (ignored).
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate filter
#' @importFrom ggplot2 ggplot geom_bar geom_errorbar geom_hline theme_bw theme element_blank element_text scale_y_continuous ggtitle coord_flip
#' @importFrom forcats fct_rev
#' @importFrom rlang .data
#'
#' @keywords internal
plot_loadings.mldfm <- function(x, var_names = NULL, ...) {
  stopifnot(inherits(x, "mldfm"))
  
  factors   <- get_factors(x)
  loadings  <- get_loadings(x)
  residuals <- get_residuals(x)
  
  t <- nrow(residuals)
  N <- ncol(residuals)
  
  loadings_df <- as.data.frame(loadings)
  
  # Factor names
  keys   <- names(x$factors_list)
  values <- unlist(x$factors_list)
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
    pivot_longer(cols = - .data$Variables, names_to = "Factor", values_to = "Loading")
  
  # Standard errors and CIs
  se_vector <- apply(residuals, 2, sd) / sqrt(t)
  
  loadings_long <- loadings_long %>%
    mutate(SE = rep(se_vector, times = length(unique(.data$Factor))),
           Loading_lower = .data$Loading - 1.96 * .data$SE,
           Loading_upper = .data$Loading + 1.96 * .data$SE)
  
  
  # Plot
  unique_factors <- unique(loadings_long$Factor)
  y_min <- -1
  y_max <- 1
  
  for (factor_name in unique_factors) {
    df_i <- loadings_long %>%
      filter(.data$Factor == factor_name & .data$Loading != 0) %>%
      mutate(Variables = factor(.data$Variables, levels = unique(.data$Variables)))
    
    p <- ggplot(df_i, aes(x = forcats::fct_rev(.data$Variables), y = .data$Loading)) +
      geom_bar(stat = "identity", fill = "grey", alpha = 0.7) +
      geom_hline(yintercept = 0, color = "red") +
      geom_errorbar(aes(ymin = .data$Loading_lower, ymax = .data$Loading_upper),
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
