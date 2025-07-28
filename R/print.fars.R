#' @title Print Method for \code{fars} Object
#'
#' @description Prints a short summary of the fars object.
#'
#' @param object An object of class \code{fars}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input \code{object}, returned invisibly. 
#'
#' @method print fars
#' @export
print.fars <- function(object, ...) {
  stopifnot(inherits(object, "fars"))
  
  cat("Factor-Augmented Quantile Regressions (FARS)\n")
  cat("===========================================\n")
  
  # Summary of forecasted quantiles
  cat("Forecasted quantiles\n")
  cat("Number of periods: ", nrow(object$quantiles), "\n")
  cat("Quantile levels: ", formatC(object$levels, format = "f", digits = 2), "\n")
  
  # Check if stressed quantiles are available
  if (!is.null(object$stressed_quantiles)) {
    cat("Stressed quantiles: YES\n")
  } else {
    cat("Stressed quantiles: NO\n")
  }
  
  invisible(object)
}
