#' Function to generate an AR(1) variance-covariance matrix with parameter rho s.t. \eqn{\lvert \rho\rvert  < 1}. 
#'  
#' @param rho Numeric value representing the autocorrelation parameter. Must satisfy |rho| < 1.
#' @param n Integer representing the size of the matrix (n x n).
#' @keywords internal 
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm
#' @return A variance-covariance matrix of size n x n based on the AR(1) process.

ARcov <- function(rho, n) {
  
  # Check if rho is within the valid range
  if (abs(rho) >= 1) {
    stop("Parameter rho must be between -1 and 1.")
  }
  
  # Sequence of powers of rho
  sqnc <- rho^seq(0, n, by = 1)
  
  # Construct Toeplitz matrix
  Omega <- toeplitz(sqnc[1: n])
  
  # Compute variance-covariance matrix
  sig <- (1/(1 - rho^2)) * Omega
  
  return(sig)
}
