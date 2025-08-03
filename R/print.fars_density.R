#' @title Print Method for \code{fars_density} Object
#'
#' @description Displays a brief summary of the \code{fars_density} object.
#'
#' @param x An object of class \code{fars_density}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input \code{fars_density} object, invisibly.
#' 
#' @method print fars_density
#' @export
print.fars_density <- function(x, ...) {
  stopifnot(inherits(x, "fars_density"))
  cat("FARS Density\n")
  cat("====================\n")
  cat("Time observations  :", nrow(x$density), "\n")
  cat("Estimation points  :", ncol(x$density), "\n")
  cat("Random samples     :", ncol(get_distribution(x)), "\n")
  cat("Support range      : [", min(x$eval_points), ",", max(x$eval_points), "]\n")
  cat("Optimization       :", x$optimization,"\n")
  invisible(x)
}

