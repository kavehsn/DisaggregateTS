#' Function to generate an ARIMA(1,1,0) variance-covariance matrix for the Litterman method with parameter \eqn{\rho} such that \eqn{\lvert \rho \rvert < 1}.
#'  
#' @param rho Numeric value representing the autocorrelation parameter. Must satisfy \eqn{\lvert \rho \rvert < 1}.
#' @param n Integer representing the size of the matrix \eqn{(n x n)}.
#' @keywords internal 
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm
#' @return A variance-covariance matrix of size \eqn{(n x n)} for the ARIMA(1,1,0) process, used in the Litterman method.

ARcov_lit <- function(rho, n) {

  # Check if rho is within the valid range
  if (abs(rho) >= 1) {
    stop("Parameter rho must be between -1 and 1.")
  }

  # Construct the Delta_t matrix
  diags <- list(rep(1, times = n), rep(-1, times = n-1))
  Delta_t <- bandSparse(n, k = 0:1, diagonals = diags, symmetric = FALSE)
  
  # Transpose to get Delta matrix
  Delta <- t(Delta_t)
  
  # Construct the H_r_t matrix
  diags2 <- list(rep(1, times = n), rep(-rho, times = n-1))
  H_r_t <- bandSparse(n, k = 0:1, diagonals = diags2, symmetric = FALSE)
  
  # Transpose to get H_r matrix
  H_r <- t(H_r_t)
  
  # Compute the variance-covariance matrix
  Sigma <- solve(Delta_t %*% H_r_t %*% H_r %*% Delta)
  
  # Force the result to be symmetric
  Sigma <- forceSymmetric(Sigma)
  
  return(Sigma)
}
