#' BIC Score for Sparse Temporal Disaggregation
#'
#' This function calculates the BIC score for sparse temporal disaggregation, as described in
#' \insertCite{10.1111/rssa.12952;textual}{DisaggregateTS}. It uses the LARS algorithm to find the
#' optimal beta coefficients and refits the models to compute BIC scores.
#'
#' @param Y        The low-frequency response series (\eqn{n_l \times 1} matrix).
#' @param X        The aggregated high-frequency indicator series (\eqn{n_l \times p} matrix).
#' @param vcov     Aggregated variance-covariance matrix of AR(\eqn{1}) residuals (\eqn{n_l \times n_l} matrix).
#' @return         The minimum BIC score from the refitted models.
#' @keywords internal
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm

sptd_BIC <- function(Y, X, vcov) {

  # Error handling: Ensure that Y, X, and vcov have compatible dimensions
  n_l <- dim(Y)[1]
  if (n_l != dim(X)[1]) stop("Y and X must have the same number of rows.")
  if (n_l != dim(vcov)[1] || n_l != dim(vcov)[2]) stop("vcov must be a square matrix with dimensions matching Y.")

  # Cholesky factorization of the covariance matrix
  Uchol <- chol(vcov)
  Lchol <- t(Uchol)

  # Precondition the variables
  X_F <- solve(Lchol) %*% X
  Y_F <- solve(Lchol) %*% Y

  # Fit LARS algorithm to the preconditioned data
  lars.fit <- lars(X_F, Y_F, intercept = FALSE, normalize = FALSE)
  betamat <- lars.fit$beta

  # Limit the support size to n_l/2 to prevent overfitting
  max_path <- k.index(betamat, n_l)

  # Initialize list for refitted beta values and BIC score array
  beta_refit <- list()
  BIC <- numeric(max_path)

  # Compute BIC for the first model
  BIC[1] <- hdBIC(X_F, Y_F, vcov, betamat[1, ])
  beta_refit[[1]] <- betamat[1, ]

  # Loop through each lambda path to compute BIC for refitted beta values
  for (lam in 2:max_path) {
    beta_refit[[lam]] <- refit(X_F, Y_F, betamat[lam, ])
    BIC[lam] <- hdBIC(X_F, Y_F, vcov, beta_refit[[lam]])
  }

  # Return the minimum BIC value
  return(min(BIC))
}
