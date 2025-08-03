#' @title Get Ellipsoids from \code{fars_scenario}
#'
#' @description Returns the list of ellipsoids from a \code{fars_scenario} object.
#'
#' @param x An object of class \code{fars_scenario}.
#'
#' @return A list of matrices defining the ellipsoids at each time.
#' @export
get_ellipsoids <- function(x) {
  stopifnot(inherits(x, "fars_scenario"))
  x$ellipsoids
}
