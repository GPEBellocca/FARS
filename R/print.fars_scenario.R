#' @title Print Method for \code{fars_scenario} Object
#'
#' @description Prints a short summary of the FARS scenario object.
#' 
#' @param object An object of class \code{fars_scenario}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input \code{fars_scenario} object, invisibly.
#'
#' @method print fars_scenario
#' @export
print.fars_scenario <- function(object, ...) {
  stopifnot(inherits(object, "fars_scenario"))
  
  cat("FARS Scenario\n")
  cat("=====================\n")
  cat("Number of periods    :", object$periods, "\n")
  cat("Ellipsoid dimensions :", ncol(object$center), "\n")
  cat("Points per ellipsoid :", object$n_points, "\n")
  cat("Confidence level     :", round(object$alpha * 100), "%\n")
  cat("FPR Gamma            :", ifelse(isTRUE(object$call$fpr), "TRUE", "FALSE"), "\n")
  
  
  invisible(object)
}
