#' @title Plot Factors from \code{mldfm} Object
#'
#' @description Displays time series plots of the estimated factors with 95% confidence bands.
#'
#' @param object An object of class \code{mldfm}.
#' @param dates Optional vector of dates. If NULL, uses 1:n as default.
#' @param ... Additional arguments (ignored).
#' 
#' @importFrom dplyr mutate filter
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon ggtitle coord_cartesian theme_bw theme element_blank element_text scale_y_continuous
#' @importFrom MASS ginv
#' @importFrom magrittr %>%
#'
#' @keywords internal
plot_factors.mldfm <- function(object, dates = NULL, ...) {
  stopifnot(inherits(object, "mldfm"))
  
  factors   <- get_factors(object)
  loadings  <- get_loadings(object)
  residuals <- get_residuals(object)
  
  T_obs  <- nrow(residuals)
  N_vars <- ncol(residuals)
  
  # Compute standard deviation for confidence bands
  PP      <- MASS::ginv((t(loadings) %*% loadings) / N_vars)
  sigma_e <- sum(diag(t(residuals) %*% residuals)) / (N_vars * T_obs)
  gamma   <- sigma_e * (t(loadings) %*% loadings) / N_vars
  SD      <- sqrt(diag(PP %*% gamma %*% PP) / N_vars)
  
  # Factor names
  keys         <- names(object$factors_list)
  values       <- unlist(object$factors_list)
  factor_names <- unlist(
    mapply(function(key, val) {
      clean <- paste0("F", gsub("-", "", key))
      if (val > 1) paste0(clean, "n", seq_len(val)) else clean
    }, keys, values, SIMPLIFY = FALSE)
  )
  if (is.null(factor_names) || length(factor_names) != ncol(factors)) {
    factor_names <- paste0("F", seq_len(ncol(factors)))
  }
  colnames(factors) <- factor_names
  
  if (is.null(dates)) {
    dates <- seq_len(nrow(factors))
  }
  
  df_long <- as.data.frame(factors) %>%
    mutate(Date = as.Date(dates)) %>%
    pivot_longer(cols = -Date, names_to = "Factor", values_to = "Value") %>%
    mutate(index = as.numeric(factor(Factor, levels = factor_names)),
           LB = Value - 2 * SD[index],
           UB = Value + 2 * SD[index])
  
  y_min <- min(df_long$LB, na.rm = TRUE)
  y_max <- max(df_long$UB, na.rm = TRUE)
  
  for (factor_name in factor_names) {
    df_i <- df_long %>% filter(Factor == factor_name)
    
    p <- ggplot(df_i, aes(x = Date, y = Value)) +
      geom_line(color = "blue", alpha = 0.6) +
      geom_ribbon(aes(ymin = LB, ymax = UB), fill = "grey70", alpha = 0.3) +
      ggtitle(factor_name) +
      coord_cartesian(ylim = c(y_min, y_max)) +
      theme_bw() +
      theme(
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none"
      )
    
    print(p)
  }
  
  invisible(NULL)
}