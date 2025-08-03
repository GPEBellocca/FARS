#' Print Method for \code{mldfm} Object
#'
#' @description Prints a short summary of the MLDFM object.
#' 
#' @param x An object of class \code{mldfm}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input \code{mldfm} object, invisibly.
#'
#' @method print mldfm
#' @export
print.mldfm <- function(x, ...) {
  stopifnot(inherits(x, "mldfm"))
  factors <- get_factors(x)
  
  cat("Multilevel Dynamic Factor Model (MLDFM)\n")
  cat("=======================================\n")
  cat("Number of periods:", nrow(factors), "\n")
  cat("Number of factors:", ncol(factors), "\n")
  cat("Number of nodes  :", length(x$factors_list), "\n")
  
  
  
  cat("\nFactor structure:\n")
  for (key in names(x$factors_list)) {
    cat(" -", key, ": ", x$factors_list[[key]], "factor(s)\n")
  }
  
  
  invisible(x)
}
