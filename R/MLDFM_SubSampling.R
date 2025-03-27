#' Subsampling Procedure for MLDFM Estimation
#'
#' Repeatedly applies the MLDFM estimation to randomly drawn subsamples of the input data.
#'
#' @param data A numeric matrix or data frame containing the time series data. Rows represent time points; columns represent observed variables.
#' @param blocks Integer. The number of blocks into which the data is divided.
#' @param block_ind A vector of integers indicating the end index of each block. Must be of length \code{blocks} and in increasing order. Required if \code{blocks > 1}.
#' @param r A vector of integers specifying the number of factors to extract for each node in the block hierarchy. Its length must equal \code{2^blocks - 1}, corresponding to all nodes in the hierarchical tree.
#' @param method Integer. The method used to initialize the factors: \code{0} for Canonical Correlation Analysis (CCA), \code{1} for Principal Component Analysis (PCA).
#' @param tol Numeric. The tolerance level for the residual sum of squares (RSS) minimization process. Used as a convergence criterion.
#' @param max_iter Integer. The maximum number of iterations allowed for the RSS minimization process.
#' @param n_sample Number of subsamples to generate.
#' @param sample_size Proportion of the original sample to retain (e.g., 0.9 for 90%).
#' @param seed Optional integer. Seed for reproducibility of the subsampling process. If \code{NULL}, random draws will differ each run.
#'
#' @return A list of \code{mldfm} objects, one for each subsample.
#' @export

mldfm_subsampling <- function(data, blocks = 1, block_ind = NULL, r = c(1), 
                              method = 0, tol = 1e-6, max_iter = 1000, 
                              n_sample = 10, sample_size = 1, seed = NULL) {
  
  
  
  # Argument checks
  if (!is.matrix(data) && !is.data.frame(data)) stop("data must be a matrix or data frame.")
  if (!is.numeric(blocks) || length(blocks) != 1) stop("blocks must be a single numeric value.")
  if (!is.numeric(r) || length(r) != (2^blocks - 1)) stop("r must be a numeric vector of length 2^blocks - 1.")
  if (!is.numeric(tol) || tol <= 0) stop("tol must be a positive numeric value.")
  if (!is.numeric(max_iter) || max_iter < 1) stop("max_iter must be a positive integer.")
  if (!method %in% c(0, 1)) stop("method must be 0 (CCA) or 1 (PCA).")
  if (blocks > 1 && is.null(block_ind)) stop("block_ind must be provided when blocks > 1.")
  if (!is.numeric(n_sample) || n_sample < 1) stop("n_sample must be a positive integer.")
  if (!is.numeric(sample_size) || sample_size <= 0 || sample_size > 1) stop("sample_size must be a number in (0, 1].")
  
  
  n_obs <- nrow(data)
  result <- vector("list", n_sample)
  
 
  for (i in 1:n_sample) {
    
    # Compute subsample 
    sub_sample_result <- compute_subsample(data, 
                                   block_ind = block_ind, 
                                   n = n_blocks,
                                   r = r,
                                   sample_size,
                                   seed = if (!is.null(seed)) seed + i - 1 else NULL )
    
    # Call MLDFM on subsample
    mldfm_result <- mldfm(sub_sample_result$sample_data, 
                          blocks = blocks, 
                          block_ind = sub_sample_result$sample_block_ind,
                          r = r, 
                          method = method, 
                          tol = tol, 
                          max_iter = max_iter)
    
    # Store results
    result[[i]] <- mldfm_result
  }
  
  
  return(result)
}
  
