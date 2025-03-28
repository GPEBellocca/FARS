#' Generate Plot for a Single Factor
#'
#' @import dplyr
#' @import ggplot2
#' @keywords internal
generate_factor_plot <- function(factor_name, data_plot, ymin, ymax) {
  data_plot %>%
    filter(Factors == factor_name) %>%
    ggplot(aes(x = Date, y = value)) +
    geom_line(color = "blue", alpha = 0.5) +
    geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.3) +
    facet_wrap(vars(Factors), nrow = 3) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_blank()
    )
}

#' Plot Extracted Factors from MLDFM
#'
#' @export
plot_factors.mldfm <- function(x, dates = NULL) {
  Factors <- x$Factors
  Lambda <- x$Lambda
  Residuals <- x$Residuals
  
  T_obs <- nrow(Residuals)
  N_vars <- ncol(Residuals)
  
  # Confidence intervals (assuming uncorrelated idiosyncratic components)
  PP <- MASS::ginv((t(Lambda) %*% Lambda) / N_vars)
  
  sigma_e <- sum(diag(t(Residuals) %*% Residuals)) / (N_vars * T_obs)
  Gamma <- sigma_e * (t(Lambda) %*% Lambda) / N_vars
  SD <- sqrt(diag(PP %*% Gamma %*% PP) / N_vars)
  
  # Factor names
  factor_df <- as.data.frame(Factors)
  keys <- names(x$Factors_list)
  values <- unlist(x$Factors_list)
  
  factor_names <- unlist(
    mapply(function(key, val) {
      clean <- paste0("F", gsub("-", "", key))
      if (val > 1) paste0(clean, "n", seq_len(val)) else clean
    }, keys, values, SIMPLIFY = FALSE)
  )
  
  # Fallback: assign default column names if needed (1 factor case)
  if (is.null(factor_names) || length(factor_names) != ncol(factor_df)) {
    factor_names <- paste0("F", seq_len(ncol(factor_df)))
  }
  colnames(factor_df) <- factor_names
  
  # Add date
  if (is.null(dates)) dates <- 1:nrow(factor_df)
  
  df_long <- factor_df %>%
    mutate(Date = as.Date(dates)) %>%
    pivot_longer(cols = -Date, names_to = "Factors", values_to = "value") %>%
    mutate(index = as.numeric(factor(Factors, levels = factor_names)),
           LB = value - 2 * SD[index],
           UB = value + 2 * SD[index])
  
  
  # Compute global scale for ribbon
  y_min <- min(df_long$LB, na.rm = TRUE)
  y_max <- max(df_long$UB, na.rm = TRUE)
  
  # Create plots
  plot_list <- lapply(factor_names, function(fname) {
    generate_factor_plot(fname, df_long, y_min, y_max)
  })
  
  names(plot_list) <- factor_names
  return(plot_list)
}