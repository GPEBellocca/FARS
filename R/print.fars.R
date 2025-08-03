#' @title Print Method for \code{fars} Object
#'
#' @description Prints a short summary of the fars object.
#'
#' @param x An object of class \code{fars}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input \code{fars} object, returned invisibly. 
#'
#' @method print fars
#' @export
print.fars <- function(x, ...) {
  stopifnot(inherits(x, "fars"))
  
  cat("Factor-Augmented Quantile Regressions (FARS)\n")
  cat("===========================================\n")
  
  # Summary of forecasted quantiles
  cat("Forecasted quantiles\n")
  cat("Number of periods: ", nrow(x$quantiles), "\n")
  cat("Quantile levels: ", formatC(x$levels, format = "f", digits = 2), "\n")
  
  # Check if stressed quantiles are available
  if (!is.null(x$stressed_quantiles)) {
    cat("Stressed quantiles: YES\n")
  } else {
    cat("Stressed quantiles: NO\n")
  }
  
  invisible(x)
}
