#' @title Summary Method for \code{mldfm} Object
#'
#' @description Provides a complete summary of the MLDFM object.
#'
#' @param object An object of class \code{mldfm}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input \code{mldfm} object, invisibly. 
#'
#' @method summary mldfm
#' @export
summary.mldfm <- function(object, ...) {
  stopifnot(inherits(object, "mldfm"))
  factors <- get_factors(object)
  residuals <- get_residuals(object)
  
  cat("Summary of Multilevel Dynamic Factor Model (MLDFM)\n")
  cat("===================================================\n")
  cat("Number of periods               :", nrow(factors), "\n")
  cat("Number of factors               :", ncol(factors), "\n")
  cat("Number of nodes                 :", length(object$factors_list), "\n")
  
  if (!is.null(object$method)) {
    cat("Initialization method           :", object$method, "\n")
  }
  
  if (!is.null(object$iterations)) {
    cat("Number of iterations to converge:", object$iterations, "\n")
  }
  
  cat("\nFactor structure:\n")
  for (key in names(object$factors_list)) {
    cat(" -", key, ": ", object$factors_list[[key]], "factor(s)\n")
  }
  
  if (!is.null(residuals)) {
    rss <- sum(residuals^2)
    avg_rss <- mean(rowSums(residuals^2))
    cat("\nResidual diagnostics:\n")
    cat(" - Total residual sum of squares (RSS): ", formatC(rss, format = "f", digits = 2), "\n")
    cat(" - Average RSS per time period        : ", formatC(avg_rss, format = "f", digits = 2), "\n")
  }
  
  invisible(object)
}