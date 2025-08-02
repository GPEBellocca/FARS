#' @title Print Method for \code{fars_density} Object
#'
#' @description Displays a brief summary of the \code{fars_density} object.
#'
#' @param object An object of class \code{fars_density}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input \code{fars_density} object, invisibly.
#' 
#' @method print fars_density
#' @export
print.fars_density <- function(object, ...) {
  stopifnot(inherits(object, "fars_density"))
  cat("FARS Density\n")
  cat("====================\n")
  cat("Time observations  :", nrow(object$density), "\n")
  cat("Estimation points  :", ncol(object$density), "\n")
  cat("Random samples     :", ncol(get_distribution(object)), "\n")
  cat("Support range      : [", min(object$eval_points), ",", max(object$eval_points), "]\n")
  cat("Optimization       :", object$optimization,"\n")
  invisible(object)
}

