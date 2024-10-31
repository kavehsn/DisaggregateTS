#' Likelihood function for Chow-Lin or Litterman temporal disaggregation.
#'
#' This function computes the likelihood function used in temporal disaggregation to find the optimal \eqn{\rho} parameter.
#' It is used in conjunction with \code{\link{disaggregate}} to estimate the autocorrelation coefficient \eqn{\rho}.
#'
#' @param Y        The low-frequency response series (an \eqn{n_l \times 1} matrix).
#' @param X        The aggregated high-frequency indicator series (an \eqn{n_l \times p} matrix).
#' @param vcov     Aggregated variance-covariance matrix for the Chow-Lin or Litterman residuals.
#' @return         The log-likelihood value for the given parameters.
#' @keywords internal
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm
chowlin_likelihood <- function(Y, X, vcov) {

  # Error handling: Ensure that dimensions of Y, X, and vcov are compatible
  n_l <- dim(Y)[1]
  if (n_l != dim(X)[1]) stop("Dimensions of Y and X do not match.")
  if (n_l != dim(vcov)[1] || n_l != dim(vcov)[2]) stop("vcov must be a square matrix of dimensions matching Y.")

  # Cholesky factorization of the covariance matrix

  Uchol <- chol(vcov)
  Lchol <- t(Uchol)

  # Preconditioning the variables using the Cholesky factorization
  X_F <- solve(Lchol) %*% X
  Y_F <- solve(Lchol) %*% Y

  # Estimate betaHat using GLS assuming the covariance structure from vcov
  betaHat <- solve(t(X_F) %*% X_F) %*% t(X_F) %*% Y_F

  # Compute residuals based on the estimated betaHat
  u_l_sim <- Y - X %*% betaHat

  # Precondition the residuals
  u_l_sim_F <- solve(Lchol) %*% u_l_sim

  # Calculate the log-likelihood function
  log_det <- sum(log(diag(Lchol)))  # log-determinant of Lchol
  LF <- -(n_l / 2) * log(2 * pi) - log_det - n_l / 2 - (n_l / 2) * log((t(u_l_sim_F) %*% u_l_sim_F) / n_l)

  return(LF)
}
