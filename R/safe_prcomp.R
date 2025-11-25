#' @title Safe principal component calculation
#'
#' @description
#' Wrapper around [stats::prcomp()] that first tries the standard
#' LAPACK-based PCA and, if it fails (for example due to a
#' `La.svd` / `dgesdd` error), falls back to more robust
#' alternatives. It first tries [svd()] on the data matrix,
#' and if that also fails, it computes the principal components
#' via an eigendecomposition of \eqn{X'X}. This is useful in
#' simulation or subsampling contexts where some matrices of
#' residuals may be nearly singular or numerically problematic.
#'
#' @param x A numeric matrix for which principal components should be
#'   computed. Must not contain `NA`, `NaN` or `Inf` values.
#' @param k Optional integer specifying the maximum number of right
#'   singular vectors (principal component loadings) to retain.
#'   If NULL, it defaults to min(nrow(x), ncol(x)). The value is
#'   always truncated to min(nrow(x), ncol(x)).
#'
#' @return
#' An object of class \code{prcomp} with at least the components:
#' \itemize{
#'   \item{sdev} The singular values (standard deviations of the components).
#'   \item{rotation} The matrix of variable loadings (right singular vectors).
#' }
#'
#' @importFrom stats prcomp
#'
#' @keywords internal
safe_prcomp <- function(x, k = NULL) {
  if (any(!is.finite(x))) {
    stop("safe_prcomp: x contains NA/NaN/Inf")
  }
  
  x  <- as.matrix(x)
  nr <- nrow(x); nc <- ncol(x)
  
  if (is.null(k)) {
    k <- min(nr, nc)
  }
  k <- min(k, nr, nc)
  
  ## 1)  prcomp
  pc <- try(prcomp(x, center = FALSE, scale. = FALSE), silent = TRUE)
  if (!inherits(pc, "try-error")) {
    return(pc)
  }
  
  ## 2) svd on X
  sv <- try(svd(x, nu = 0, nv = k), silent = TRUE)
  if (!inherits(sv, "try-error")) {
    k_eff <- min(k, length(sv$d), ncol(sv$v))
    return(
      structure(
        list(
          sdev     = sv$d[seq_len(k_eff)],
          rotation = sv$v[, seq_len(k_eff), drop = FALSE]
        ),
        class = "prcomp"
      )
    )
  }
  
  ## 3) PCA via eigen(X'X)
  S   <- crossprod(x)  # p x p
  eig <- try(eigen(S, symmetric = TRUE), silent = TRUE)
  if (!inherits(eig, "try-error")) {
    vals <- eig$values
    vals[vals < 0] <- 0
    k_eff <- min(k, length(vals))
    
    rotation <- eig$vectors[, seq_len(k_eff), drop = FALSE]
    sdev     <- sqrt(vals[seq_len(k_eff)])
    
    return(
      structure(
        list(
          sdev     = sdev,
          rotation = rotation
        ),
        class = "prcomp"
      )
    )
  }
  
  ## 4) All method failed. Diagnostic.
  rng <- try(range(x[is.finite(x)]), silent = TRUE)
  msg <- paste0(
    "safe_prcomp: prcomp(), svd() and eigen(crossprod(x)) all failed.\n",
    "  dim(x): ", paste(dim(x), collapse = " x "), "\n",
    "  any(!is.finite(x)): ", any(!is.finite(x)), "\n",
    "  range(finite x): ",
    if (inherits(rng, "try-error")) "unavailable" else paste(rng, collapse = " ")
  )
  stop(msg, call. = FALSE)
}
