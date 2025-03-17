orthonormalize_factors <- function(factors) {
  factors <- as.matrix(factors)
  qr_decomp <- qr(factors)
  Q <- qr.Q(qr_decomp)
  #Q <- Q * sqrt(nrow(factors))  # Scale to ensure F'F/T = I
  return(Q)
}