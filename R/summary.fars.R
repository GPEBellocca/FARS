#' @title Summary Method for \code{fars} Object
#'
#' @description Prints a complete summary of the fars object, including information on estimated quantiles, stressed quantiles,
#' regression coefficients, standard errors, and p-values.
#'
#' @param object An object of class \code{fars}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input \code{fars} object, returned invisibly.
#'
#' @method summary fars
#' @export
summary.fars <- function(object, ...) {
  cat("Factor-Augmented Quantile Regressions (FARS)\n")
  cat("===========================================\n")
  cat("Summary of Quantile Regressions\n\n")
  
  levels <- object$levels
  coeff <- object$coeff
  stderr <- object$std_error
  pval <- object$pvalue
  variables <- rownames(coeff)
  
  for (i in seq_along(levels)) {
    cat(sprintf("-------------------------\n Quantile: %.2f \n-------------------------\n", levels[i]))
    est <- coeff[, i]
    se <- stderr[, i]
    p <- pval[, i]
    
    summary_df <- data.frame(
      Estimate = formatC(est, digits = 3, format = "f"),
      `Std. Error` = formatC(se, digits = 3, format = "f"),
      `P-value` = formatC(p, digits = 3, format = "f"),
      row.names = variables,
      check.names = FALSE
    )
    
    print(summary_df)
    cat("\n")
  }
  
  invisible(object)
}
