#' @title Print Method for \code{mldfm_subsample} Object
#'
#' @description Prints a brief summary of the \code{mldfm_subsample} object.
#'
#' @param object An object of class \code{mldfm_subsample}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object \code{object}, invisibly.
#'
#' @method print mldfm_subsample
#' @export
print.mldfm_subsample <- function(object, ...) {
  stopifnot(inherits(object, "mldfm_subsample"))
  
  cat("MLDFM Subsampling\n")
  cat("==========================\n")
  cat("Number of subsamples :", object$n_samples, "\n")
  cat("Sample size fraction :", object$sample_size, "\n")
  if (!is.null(object$seed)) cat("Seed used            :", object$seed, "\n")
  
  models <- get_mldfm_list(object)
  if (length(models) > 0 && inherits(models[[1]], "mldfm")) {
    model1 <- get_mldfm_model(object, 1)
    T_obs <- nrow(model1$residuals)
    N_vars <- ncol(model1$residuals)
    cat("Data dimensions      :", T_obs, "periods,", N_vars, "variables\n")
    
    cat("\nFactor structure:\n")
    f_list <- model1$factors_list
    for (key in names(f_list)) {
      cat(" •", key, ": ", f_list[[key]], "factor(s)\n")
    }
  }
  
  invisible(object)
}
