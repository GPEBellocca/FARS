#' @title Get List of MLDFMs from a \code{mldfm_subsample} Object
#'
#' @description Returns the list of all \code{mldfm} stored in a \code{mldfm_subsample} object.
#'
#' @param object An object of class \code{mldfm_subsample}.
#' @param ... Additional arguments (ignored).
#'
#' @return A list of \code{mldfm} objects.
#' @export
get_mldfm_list <- function(object, ...) {
  UseMethod("get_mldfm_list")
}

#' @export
get_mldfm_list.mldfm_subsample <- function(object, ...) {
  stopifnot(inherits(object, "mldfm_subsample"))
  object$models
}