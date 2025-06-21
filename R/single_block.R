#' Multi-Level Dynamic Factor Model - Single block
#'
#'
#' @importFrom stats prcomp
#'
#' @keywords internal

# Internal function: Single-block MLDFM computation
single_block <- function(data, r) {
  
  # Standardize the data 
  X <- scale(data, center = TRUE, scale = TRUE)

 
  
  # Dimensions
  T <- nrow(X)
  N <- ncol(X)
  
  # Eigen decomposition
  eig_res <- eigen(X %*% t(X))
  values <- eig_res$values[1:r]
  vectors <- eig_res$vectors[, 1:r, drop = FALSE]

  # Extract factors (F_tilde) and loadings (P_tilde)
  Factors <- sqrt(T) * vectors
  Lambda <- (1 / T) * t(Factors) %*% X

  
  # pca_res <- prcomp(X, center = FALSE, scale. = FALSE)
  # pca_Factors <- pca_res$x[, 1:r]      # T × r
  # pca_Lambda <- t(pca_res$rotation[, 1:r])  # r × N
  # 
  # Factors <- as.matrix(Factors)
  # pca_Factors <- as.matrix(pca_Factors)
  # for (i in 1:ncol(Factors)) {
  #   corr <- cor(Factors[, i], pca_Factors[, i])
  #   cat(sprintf("Fattore %d: correlazione = %.4f\n", i, corr))
  # }
  # 
  
  # Residuals 
  Residuals <- X - Factors %*% Lambda

  
  # Save factor structure 
  Factors_list <- list()
  Factors_list[["1"]] <- r  # list index must be numeric for proper iteration later
  
  
  iteration <- 0
  
  return(list(
    Factors = Factors,
    Lambda = t(Lambda),
    Residuals = Residuals,
    Method = "PCA",
    Iterations = iteration,
    Factors_list = Factors_list
  ))
}



