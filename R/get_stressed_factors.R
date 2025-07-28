#' @title Generic function to extract stressed factors
#'
#' @param object An object from which to extract the stressed factors.
#' @param ... Additional arguments.
#'
#' @return A matrix of stressed factors.
#' @export
get_stressed_factors <- function(object, ...) {
  UseMethod("get_stressed_factors")
}
#' @title Extract Stressed Factors from a \code{fars} Object
#'
#' @description Extracts the stressed factors from a \code{fars} object. If stressed factors are not available,
#'             it returns NULL.
#'
#' @param object An object of class \code{fars}.
#' @param ... Further arguments (ignored).
#'
#' @return A matrix containing the stressed factors if available, otherwise NULL.
#'
#' @examples
#' fars_result <- compute_fars(dep_variable = rnorm(100), factors = matrix(rnorm(100 * 3), ncol = 3))
#' get_stressed_factors(fars_result)  
#'
#' @export
get_stressed_factors.fars <- function(object, ...) {
  stopifnot(inherits(object, "fars"))
  
  # If stressed factors are available
  if (!is.null(object$stressed_factors)) {
    return(object$stressed_factors)
  }
  
  # If stressed factors are not available
  return(NULL)
}
