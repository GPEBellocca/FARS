#' @title Get Sigma List from \code{fars_scenario}
#'
#' @description Returns the list of covariance matrices used to construct the ellipsoids.
#'
#' @param object An object of class \code{fars_scenario}.
#'
#' @return A list of covariance matrices (one per period).
#' @export
get_sigma_list <- function(object) {
  stopifnot(inherits(object, "fars_scenario"))
  object$sigma
}
