#' @title Plot Method for \code{mldfm_subsample} Object
#'
#' @description Plots a histogram of the number of iterations used in each subsample estimation.
#'
#' @param object An object of class \code{mldfm_subsample}.
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object (invisibly).
#'
#' @importFrom ggplot2 ggplot aes geom_histogram labs theme_minimal
#' @importFrom ggplot2 element_text
#' @method plot mldfm_subsample
#' @export
plot.mldfm_subsample <- function(object, ...) {
  stopifnot(inherits(object, "mldfm_subsample"))
  
  iterations <- sapply(object$models, function(m) m$iterations)
  df <- data.frame(Iterations = iterations)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Iterations)) +
    ggplot2::geom_histogram(binwidth = 1, fill = "steelblue", color = "white", ...) +
    ggplot2::labs(
      title = "Sequential Least Squares Iterations",
      x = "Number of Iterations",
      y = "Frequency"
    ) +
    ggplot2::theme_minimal()
  
  print(p)
  invisible(p)
}
