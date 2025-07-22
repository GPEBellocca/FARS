#' @title Summary Method for \code{mldfm_subsample} Object
#'
#' @description Provides a structured summary of a \code{mldfm_subsample} object.
#'
#' @param object An object of class \code{mldfm_subsample}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object \code{object}, invisibly.
#'
#' @method summary mldfm_subsample
#' @export
summary.mldfm_subsample <- function(object, ...) {
  stopifnot(inherits(object, "mldfm_subsample"))
  
  cat("Summary of MLDFM Subsampling\n")
  cat("=============================\n")
  cat("Number of subsamples :", object$n_samples, "\n")
  cat("Sample size fraction :", object$sample_size, "\n")
  if (!is.null(object$seed)) cat("Seed used            :", object$seed, "\n")
  
  if (length(object$models) > 0 && inherits(object$models[[1]], "mldfm")) {
    T_obs <- nrow(get_residuals(object$models[[1]]))
    N_vars <- ncol(get_residuals(object$models[[1]]))
    cat("Data dimensions      :", T_obs, "periods,", N_vars, "variables\n")
    
    # Factor structure from first model
    cat("\nFactor structure:\n")
    f_list <- object$models[[1]]$factors_list
    for (key in names(f_list)) {
      cat(" •", key, ": ", f_list[[key]], "factor(s)\n")
    }
    
    # Estimation method 
    cat("\nEstimation method:", object$models[[1]]$method, "\n")
    
    # Iterations 
    iterations <- sapply(object$models, function(m) m$iterations)
    cat("Iterations       : mean =", round(mean(iterations), 2),
        "| min =", min(iterations),
        "| max =", max(iterations), "\n")
    
    # Compute RSS from residuals
    rss_vals <- sapply(object$models, function(m) {
      res <- get_residuals(m)
      sum(res^2, na.rm = TRUE)
    })
    cat("Final RSS        : mean =", round(mean(rss_vals), 2),
        "| min =", round(min(rss_vals), 2),
        "| max =", round(max(rss_vals), 2), "\n")
  }
  
  invisible(object)
}
