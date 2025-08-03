#' @title Print Method for \code{fars_scenario} Object
#'
#' @description Prints a short summary of the FARS scenario object.
#' 
#' @param x An object of class \code{fars_scenario}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input \code{fars_scenario} object, invisibly.
#'
#' @method print fars_scenario
#' @export
print.fars_scenario <- function(x, ...) {
  stopifnot(inherits(x, "fars_scenario"))
  
  cat("FARS Scenario\n")
  cat("=====================\n")
  cat("Number of periods    :", x$periods, "\n")
  cat("Ellipsoid dimensions :", ncol(x$center), "\n")
  cat("Points per ellipsoid :", x$n_points, "\n")
  cat("Confidence level     :", round(x$alpha * 100), "%\n")
  cat("FPR Gamma            :", ifelse(isTRUE(x$call$fpr), "TRUE", "FALSE"), "\n")
  
  
  invisible(x)
}
