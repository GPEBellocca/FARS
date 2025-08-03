#' @title Get Sigma List from \code{fars_scenario}
#'
#' @description Returns the list of covariance matrices used to construct the ellipsoids.
#'
#' @param x An object of class \code{fars_scenario}.
#'
#' @return A list of covariance matrices (one per period).
#' @export
get_sigma_list <- function(x) {
  stopifnot(inherits(x, "fars_scenario"))
  x$sigma
}
