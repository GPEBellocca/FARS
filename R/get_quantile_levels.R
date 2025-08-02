#' @title Extract Quantile Levels from a \code{fars} Object
#'
#' @description Returns the quantile levels from an object of class \code{fars}.
#'
#' @param object An object of class \code{fars}, typically the result of a computation such as \code{compute_fars}.
#' @param ... Further arguments (ignored).
#'
#' @return 
#' A vector of quantile levels stored within the \code{fars} object.
#' 
#' @examples
#' fars_result <- compute_fars(dep_variable = rnorm(100), factors = matrix(rnorm(100 * 3), ncol = 3))
#' get_quantile_levels(fars_result)  
#'
#' @export
get_quantile_levels.fars <- function(object, ...) {
  stopifnot(inherits(object, "fars"))
  return(object$levels)
}
