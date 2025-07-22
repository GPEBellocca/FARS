#' @title Generic function to extract residuals
#'
#' @param object An object from which to extract the residuals.
#' @param ... Additional arguments.
#'
#' @return A matrix of residuals.
#' @export
get_residuals <- function(object, ...) {
  UseMethod("get_residuals")
}

#' @title Extract Residuals from a \code{mldfm} Object
#'
#' @param object An object of class \code{mldfm}.
#' @param ... Further arguments (ignored).
#'
#' @return A matrix containing the residuals.
#'
#' @examples
#' mldfm_result <- mldfm(data = matrix(rnorm(100 * 5), 100, 5), blocks = 1, global = 2)
#' get_residuals(mldfm_result)
#'
#' @export
get_residuals.mldfm <- function(object, ...) {
  stopifnot(inherits(object, "mldfm"))
  object$residuals
}
