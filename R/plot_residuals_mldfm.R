#' @title Plot Residuals from \code{mldfm} Object
#'
#' @description Displays a correlation heatmap of the residuals.
#' 
#' @param x An object of class \code{mldfm}.
#' @param var_names Optional vector of variable names. If NULL, default names are used.
#' @param ... Additional arguments (ignored).
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_distiller labs theme_minimal theme element_text element_blank
#' @importFrom stats cor
#'
#' @keywords internal
plot_residuals.mldfm <- function(x, var_names = NULL, ...) {
  stopifnot(inherits(x, "mldfm"))
  
  residuals <- get_residuals(x)
  n_vars <- ncol(residuals)
  
  country_names <- if (is.null(var_names)) {
    paste0("Var", seq_len(n_vars))
  } else {
    var_names
  }
  
  corr_matrix <- cor(residuals)
  rownames(corr_matrix) <- country_names
  colnames(corr_matrix) <- country_names
  
  corr_df <- as.data.frame(as.table(corr_matrix))
  colnames(corr_df) <- c("Country1", "Country2", "Correlation")
  
  Country1 <- Country2 <- Correlation <- NULL
  
  g <- ggplot(corr_df, aes(x = Country1, y = Country2, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_distiller(palette = "RdYlBu", limits = c(-1, 1), name = "Correlation") +
    labs(title = "Residual Correlation Matrix",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank(),
      axis.ticks = element_blank()
    )
  
  print(g)
  invisible(NULL)
}

