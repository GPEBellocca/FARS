#' @title Generic function to extract distribution from fars_density object
#'
#' @param object An object from which to extract the distribution.
#' @param ... Additional arguments.
#'
#' @return A matrix of random draws from the fitted skew-t distribution.
#' @export
get_distribution <- function(object, ...) {
  UseMethod("get_distribution")
}

#' @title Extract Distribution from a \code{fars_density} Object
#'
#' @description Extracts the distribution from a \code{fars_density} object. 
#'
#' @param object An object of class \code{fars_density}.
#' @param ... Further arguments (ignored).
#'
#' @return A matrix containing the random draws from the fitted skew-t distribution if available, otherwise NULL.
#'
#' @examples
#' \donttest{
#' fars_density_result <- compute_density(quantiles = matrix(rnorm(100 * 5), nrow = 100, ncol = 5))
#' get_distribution(fars_density_result)
#'}
#' @export
get_distribution.fars_density <- function(object, ...) {
  stopifnot(inherits(object, "fars_density"))
  
  if (!is.null(object$distribution)) {
    return(object$distribution)
  }
  return(NULL)
}
