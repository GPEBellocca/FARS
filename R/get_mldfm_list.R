#' @title Get List of MLDFMs from a \code{mldfm_subsample} Object
#'
#' @description Returns the list of all \code{mldfm} stored in a \code{mldfm_subsample} object.
#'
#' @param x An object of class \code{mldfm_subsample}.
#' @param ... Additional arguments (ignored).
#'
#' @return A list of \code{mldfm} objects.
#' @export
get_mldfm_list <- function(x, ...) {
  UseMethod("get_mldfm_list")
}

#' @export
get_mldfm_list.mldfm_subsample <- function(x, ...) {
  stopifnot(inherits(x, "mldfm_subsample"))
  x$models
}