#' @title Print Method for \code{mldfm_subsample} Object
#'
#' @description Prints a brief summary of the \code{mldfm_subsample} object.
#'
#' @param x An object of class \code{mldfm_subsample}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object \code{object}, invisibly.
#'
#' @method print mldfm_subsample
#' @export
print.mldfm_subsample <- function(x, ...) {
  stopifnot(inherits(x, "mldfm_subsample"))
  
  cat("MLDFM Subsampling\n")
  cat("==========================\n")
  cat("Number of subsamples :", x$n_samples, "\n")
  cat("Sample size fraction :", x$sample_size, "\n")
  if (!is.null(x$seed)) cat("Seed used            :", x$seed, "\n")
  
  models <- get_mldfm_list(x)
  if (length(models) > 0 && inherits(models[[1]], "mldfm")) {
    model1 <- get_mldfm_model(x, 1)
    T_obs <- nrow(model1$residuals)
    N_vars <- ncol(model1$residuals)
    cat("Data dimensions      :", T_obs, "periods,", N_vars, "variables\n")
    
    cat("\nFactor structure:\n")
    f_list <- model1$factors_list
    for (key in names(f_list)) {
      cat(" -", key, ": ", f_list[[key]], "factor(s)\n")
    }
  }
  
  invisible(x)
}
