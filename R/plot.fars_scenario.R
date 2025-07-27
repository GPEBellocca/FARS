#' Plot Method for \code{fars_scenario} Object
#'
#' @description Plots the hyperellipsoid for a given time observation (only for 1D or 2D cases). 
#'
#' @param x An object of class \code{fars_scenario}.
#' @param obs Integer. Time index to plot (default = 1).
#' @param ... Additional arguments passed to plotting functions.
#'
#' @method plot fars_scenario
#' @export
plot.fars_scenario <- function(x, obs = 1, ...) {
  stopifnot(inherits(x, "fars_scenario"))
  
  K <- ncol(x$center)
  T <- x$periods
  if (obs < 1 || obs > T) stop("Invalid observation index: out of bounds.")
  
  center <- x$center[obs, ]
  shape <- x$ellipsoids[[obs]]
  
  if (K == 1) {
    lower <- shape[1]
    upper <- shape[2]
    
    plot(NA,
         xlim = c(lower, upper),
         ylim = c(0.9, 1.1),
         xlab = "Factor",
         ylab = "",
         yaxt = "n",
         xaxt = "n",
         main = paste("1D Confidence Interval (t =", obs, ")"))
    
    segments(lower, 1, upper, 1, col = "lightblue", lwd = 3)
    points(center, 1, pch = 19, col = "darkblue")
    
    axis(1,
         at = round(c(lower, center, upper), 3),
         labels = round(c(lower, center, upper), 3))
  } else if (K == 2) {

    plot(shape, type = "l", col = "lightblue", lwd = 2,
         xlab = "Factor 1", ylab = "Factor 2", asp = 1,
         main = paste("2D Ellipsoid (t =", obs, ")"), ...)
    points(center[1], center[2], pch = 19, col = "darkblue")
  } else {
    warning("Plotting is only supported for 1D and 2D cases")
  }
}
