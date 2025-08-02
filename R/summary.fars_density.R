#' @title Summary Method for \code{fars_density} Object
#'
#' @description Displays a complete summary of the \code{fars_density} object.
#'
#' @param object An object of class \code{fars_density}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input \code{fars_density} object, invisibly.
#' 
#' @importFrom stats median
#' 
#' @method summary fars_density
#' @export
summary.fars_density <- function(object, ...) {
  stopifnot(inherits(object, "fars_density"))
  
  cat("FARS Density Summary\n")
  cat("=========================\n")

  means <- rowMeans(get_distribution(object))
  medians <- apply(get_distribution(object), 1, median)
  sds <- apply(get_distribution(object), 1, sd)
  
  summary_df <- data.frame(
    Mean = round(means, 4),
    Median = round(medians, 4),
    StdDev = round(sds, 4)
  )
  
  print(summary_df)
  invisible(summary_df)
}
