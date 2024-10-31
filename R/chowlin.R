#' Function to perform Chow-Lin temporal disaggregation from \insertCite{chow1971best;textual}{DisaggregateTS}
#' and its special case counterpart, Litterman \insertCite{litterman1983random;textual}{DisaggregateTS}.
#'
#' Used in \code{\link{disaggregate}} to find estimates given the optimal \eqn{rho} parameter.
#'
#' @param Y         The low-frequency response series (a \eqn{n_l \times 1} matrix).
#' @param X         The high-frequency indicator series (a \eqn{n \times p} matrix).
#' @param rho       The AR(1) residual parameter. Must be strictly between \eqn{-1} and \eqn{1}.
#' @param aggMat    Aggregation matrix method: 'first', 'sum', 'average', 'last'. Default is 'sum'.
#' @param aggRatio  Aggregation ratio, e.g. 4 for annual-to-quarterly, 3 for quarterly-to-monthly. Default is 4.
#' @param litterman Boolean. If TRUE, use Litterman variance-covariance method, otherwise use Chow-Lin. Default is FALSE.
#' @return A list containing the following elements:
#'   \itemize{
#'     \item \code{y}: Estimated high-frequency response series (an \eqn{n \times 1} matrix).
#'     \item \code{betaHat}: Estimated coefficient vector (a \eqn{p \times 1} matrix).
#'     \item \code{u_l}: Estimated aggregate residual series (an \eqn{n_l \times 1} matrix).
#'   }
#' @keywords internal
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm

chowlin <- function(Y, X, rho, aggMat = 'sum', aggRatio = 4, litterman = FALSE) {

  # Input validation
  if (abs(rho) >= 1) stop("'rho' must be strictly between -1 and 1.")
  if (aggRatio <= 0) stop("'aggRatio' must be a positive integer.")

  n_l <- dim(Y)[1]
  n <- dim(X)[1]
  p <- dim(X)[2]
  nfull <- aggRatio * n_l
  extr <- n - nfull  # number of extrapolations

  # Check that dimensions match for disaggregation
  if (n_l * aggRatio > n) {
    stop("Insufficient observations: X does not have enough rows for the given aggregation ratio.")
  }

  # Generate the aggregation matrix C based on 'aggMat' parameter
  if (aggMat == 'sum') {
    C <- kronecker(diag(n_l), matrix(data = 1, nrow = 1, ncol = aggRatio))
  } else if (aggMat == 'average') {
    C <- kronecker(diag(n_l), matrix(data = 1 / aggRatio, nrow = 1, ncol = aggRatio))
  } else if (aggMat == 'first') {
    C <- kronecker(diag(n_l), matrix(data = c(1, rep(0, times = aggRatio - 1)), nrow = 1, ncol = aggRatio))
  } else if (aggMat == 'last') {
    C <- kronecker(diag(n_l), matrix(data = c(rep(0, times = aggRatio - 1), 1), nrow = 1, ncol = aggRatio))
  } else {
    stop("Invalid 'aggMat' parameter. Choose from 'first', 'sum', 'average', 'last'.")
  }
  C <- cbind(C, matrix(0L, n_l, extr))  # Add extrapolation columns if necessary

  # Perform variance-covariance matrix calculation using either Chow-Lin or Litterman method
  vcov <- if (litterman) {
    ARcov_lit(rho, n)
  } else {
    ARcov(rho, n)
  }

  # Compute aggregated variance-covariance matrix and its Cholesky decomposition
  vcov_agg <- forceSymmetric(C %*% vcov %*% t(C))  # Ensure symmetry for decomposition
  Uchol <- chol(vcov_agg)
  Lchol <- t(Uchol)

  # Preconditioning the variables for GLS estimation
  X_l <- C %*% X
  X_F <- solve(Lchol) %*% X_l
  Y_F <- solve(Lchol) %*% Y

  # Estimate betaHat using GLS with rho-dependent variance-covariance structure
  betaHat <- solve(t(X_F) %*% X_F) %*% t(X_F) %*% Y_F

  # Compute distribution matrix D
  D <- vcov %*% t(C) %*% solve(vcov_agg)

  # Calculate aggregate residuals
  u_l <- Y - X_l %*% betaHat

  # Generate the high-frequency series
  y <- X %*% betaHat + (D %*% u_l)

  # Output results
  output <- list('y' = y, 'betaHat' = betaHat, 'u_l' = u_l)

  return(output)
}
