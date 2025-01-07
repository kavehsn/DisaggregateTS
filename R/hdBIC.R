#' High-dimensional BIC score
#'
#' This function calculates a BIC score that performs better than the ordinary BIC in high-dimensional scenarios.
#' It uses the variance estimator given in \insertCite{yu2019estimating;textual}{DisaggregateTS}.
#'
#' @param X          Aggregated indicator series matrix that has been GLS rotated (an \eqn{n_l \times p} matrix).
#' @param Y          Low-frequency response vector that has been GLS rotated (an \eqn{n_l \times 1} vector).
#' @param covariance Aggregated AR covariance matrix (an \eqn{n_l \times n_l} matrix).
#' @param beta       Estimate of the regression coefficients (a \eqn{p \times 1} vector).
#' @return           The BIC score for model comparison.
#' @keywords internal
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm

hdBIC <- function(X, Y, covariance, beta) {

  # Ensure the dimensions of X, Y, and covariance are compatible
  n_l <- length(Y)
  if (n_l != dim(X)[1]) stop("Dimensions of X and Y do not match.")
  if (n_l != dim(covariance)[1] || n_l != dim(covariance)[2]) stop("Covariance matrix must be square with dimensions matching Y.")

  # Calculate the number of non-zero coefficients in beta
  support <- sum(beta != 0)

  # Calculate residuals
  residuals <- Y - X %*% beta

  # Calculate log-likelihood (breaking down for clarity)
  log_likelihood <- -(n_l - support) / 2 -
                    (n_l / 2) * log(2 * pi / (n_l - support) * crossprod(residuals)) -
                    (1 / 2) * log(det(covariance))

  # Calculate BIC score
  BIC <- -2 * log_likelihood + log(n_l) * support

  return(BIC)
}
