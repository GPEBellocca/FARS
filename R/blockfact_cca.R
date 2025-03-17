blockfact_cca <- function(y, Nregio, r_glob, r_reg) {
  library(CCA)  # Load CCA package
  
  g <- length(Nregio)  # Number of groups
  RegInd <- c(0, cumsum(Nregio))  # Group indices for slicing
  r <- r_glob + r_reg  # Total number of factors per group (global + local)
  rsum <- cumsum(r)  # Cumulative sum of factors
  
  # Step 1: Extract Regional Factors using CCA
  fi <- c()
  for (i in 1:g) {
    # Extract variables for this group
    y_region <- y[, (RegInd[i] + 1):RegInd[i+1]]
    
    # Perform PCA to reduce dimension (ensuring we don't have more factors than variables)
    pca_result <- prcomp(y_region, center = TRUE, scale. = TRUE)
    y_region_reduced <- pca_result$x[, 1:r[i]]  # Keep r[i] principal components
    
    # If more than one group, apply CCA between this group and the first group
    if (i == 1) {
      fi <- y_region_reduced
    } else {
      cca_result <- cancor(fi, y_region_reduced)  # Canonical Correlation Analysis
      fi <- cbind(fi, y_region_reduced %*% cca_result$ycoef[, 1:r[i]])  # Extract regional factors
    }
  }
  
  # Normalize regional factors to ensure orthonormality
  fi <- fi / kronecker(matrix(1, nrow = nrow(fi), ncol = 1), t(matrix(sqrt(diag(t(fi) %*% fi)))))
  
  # Step 2: Extract Global Factors using CCA Across All Regional Factors
  fhat <- c()
  for (i in 1:(g-1)) {
    cca_result <- cancor(fi[, 1:rsum[1]], fi[, (rsum[i] + 1):rsum[i+1]])  # CCA between first block and others
    X_factors <- fi[, 1:rsum[1]] %*% cca_result$xcoef[, 1:r_glob]
    Y_factors <- fi[, (rsum[i] + 1):rsum[i+1]] %*% cca_result$ycoef[, 1:r_glob]
    fhat <- cbind(fhat, (X_factors + Y_factors) / 2)  # Average canonical factors
  }
  
  # Step 3: Final Processing of Global Factors (Optional PCA for Final Dimensionality Reduction)
  C <- t(fhat) %*% fhat
  evec <- eigen(C)$vectors
  fhatblock <- fhat %*% evec[, 1:r_glob]  # Extract the top r_glob components
  
  return(fhatblock)
}
