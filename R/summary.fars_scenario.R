#' Summary Method for \code{fars_scenario} Object
#'
#' @description Provides a detailed summary of the FARS scenario object.
#'
#' @param object An object of class \code{fars_scenario}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input \code{fars_scenario} object, invisibly.
#'
#' @method summary fars_scenario
#' @export
summary.fars_scenario <- function(object, ...) {
  stopifnot(inherits(object, "fars_scenario"))
  
  cat("FARS Scenario Summary\n")
  cat("======================\n")
  cat("Number of periods    :", object$periods, "\n")
  cat("Ellipsoid dimensions :", ncol(object$center), "\n")
  cat("Points per ellipsoid :", object$n_points, "\n")
  cat("Confidence level     :", round(object$alpha * 100), "%\n")
  cat("FPR Gamma            :", ifelse(isTRUE(object$call$fpr), "TRUE", "FALSE"), "\n")
  
 
  center_vals <- as.vector(object$center)
  cat("\nCenter (factor estimates):\n")
  cat("  Mean     :", round(mean(center_vals), 4), "\n")
  cat("  Std. Dev :", round(sd(center_vals), 4), "\n")
  cat("  Min      :", round(min(center_vals), 4), "\n")
  cat("  Max      :", round(max(center_vals), 4), "\n")
  
  
  sigma_vals <- unlist(lapply(object$sigma, function(S) diag(S)))
  cat("\nEllipsoid variability (diagonal of Sigma):\n")
  cat("  Mean     :", round(mean(sigma_vals), 4), "\n")
  cat("  Std. Dev :", round(sd(sigma_vals), 4), "\n")
  cat("  Min      :", round(min(sigma_vals), 4), "\n")
  cat("  Max      :", round(max(sigma_vals), 4), "\n")
  
  invisible(object)
}
