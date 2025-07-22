#' @title Get a Single \code{mldfm} Object from a \code{mldfm_subsample} Object
#'
#' @description Returns the \code{mldfm} object at the specified position in a \code{mldfm_subsample} object.
#'
#' @param object An object of class \code{mldfm_subsample}.
#' @param index Integer. The position of the desired model (between 1 and \code{n_samples}).
#' @param ... Additional arguments (ignored).
#'
#' @return A single \code{mldfm} object.
#' @export
get_mldfm_model <- function(object, index, ...) {
  UseMethod("get_mldfm_model")
}

#' @export
get_mldfm_model.mldfm_subsample <- function(object, index, ...) {
  stopifnot(inherits(object, "mldfm_subsample"))
  
  n_models <- length(object$models)
  if (!is.numeric(index) || length(index) != 1 || index < 1 || index > n_models) {
    stop(sprintf("`index` must be a single integer between 1 and %d", n_models))
  }
  
  object$models[[index]]
}
