#' Sparse Temporal Disaggregation
#'
#' This function performs sparse temporal disaggregation as described in \insertCite{10.1111/rssa.12952;textual}{DisaggregateTS}.
#' It estimates the high-frequency response series using LARS (Least Angle Regression) and applies either a LASSO or adaptive LASSO penalty.
#'
#' @param Y         The low-frequency response series (\eqn{n_l \times 1} matrix).
#' @param X         The high-frequency indicator series (\eqn{n \times p} matrix).
#' @param rho       The AR(\eqn{1}) residual parameter (must be strictly between \eqn{-1} and \eqn{1}).
#' @param aggMat    Aggregation matrix method ('first', 'sum', 'average', 'last'). Default is 'sum'.
#' @param aggRatio  Aggregation ratio (e.g., 4 for annual-to-quarterly, 3 for quarterly-to-monthly). Default is 4.
#' @param adaptive  Logical. If `TRUE`, use adaptive LASSO penalty. If `FALSE`, use standard LASSO penalty. Default is `FALSE`.
#' @return A list containing:
#' \itemize{
#'   \item \code{y}: Estimated high-frequency response series (\eqn{n \times 1} matrix).
#'   \item \code{betaHat}: Estimated coefficient vector (\eqn{p \times 1} matrix).
#'   \item \code{u_l}: Estimated aggregate residual series (\eqn{n_l \times 1} matrix).
#' }
#' @keywords internal
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm
sptd <- function(Y, X, rho, aggMat = 'sum', aggRatio = 4, adaptive = FALSE) {

  # Error handling for input parameters
  if (abs(rho) >= 1) stop("rho must be strictly between -1 and 1.")
  if (aggRatio <= 0) stop("aggRatio must be a positive integer.")
  if (!is.matrix(Y) || !is.matrix(X)) stop("Y and X must be matrices.")

  n_l <- dim(Y)[1]
  n <- dim(X)[1]
  p <- dim(X)[2]
  nfull <- aggRatio * n_l
  extr <- n - nfull  # number of extrapolations

  # Generate the aggregation matrix C
  if (aggMat == 'sum') {
    C <- kronecker(diag(n_l), matrix(data = 1, nrow = 1, ncol = aggRatio))
  } else if (aggMat == 'average') {
    C <- kronecker(diag(n_l), matrix(data = 1 / aggRatio, nrow = 1, ncol = aggRatio))
  } else if (aggMat == 'first') {
    C <- kronecker(diag(n_l), matrix(data = c(1, rep(0, times = aggRatio - 1)), nrow = 1, ncol = aggRatio))
  } else if (aggMat == 'last') {
    C <- kronecker(diag(n_l), matrix(data = c(rep(0, times = aggRatio - 1), 1), nrow = 1, ncol = aggRatio))
  } else {
    stop("Invalid aggMat parameter. Must be 'sum', 'average', 'first', or 'last'.")
  }

  C <- cbind(C, matrix(0L, n_l, extr))  # Add extrapolation columns if necessary
  X_l <- C %*% X

  # Cholesky factorization of the covariance matrix
  vcov <- ARcov(rho, n)
  vcov_agg <- C %*% vcov %*% t(C)
  Uchol <- chol(vcov_agg)
  Lchol <- t(Uchol)

  # Precondition the variables for LARS
  X_F <- solve(Lchol) %*% X_l
  Y_F <- solve(Lchol) %*% Y

  # Fit LARS algorithm to the data
  lars.fit <- lars(X_F, Y_F, intercept = FALSE, normalize = FALSE)
  betamat <- lars.fit$beta

  # Limit the support of beta estimates to n_l/2
  max_path <- k.index(betamat, n_l)

  # Find BIC for each re-fitted beta estimate
  beta_refit <- list()
  BIC <- numeric(max_path)
  BIC[1] <- hdBIC(X_F, Y_F, vcov_agg, betamat[1, ])
  beta_refit[[1]] <- betamat[1, ]

  for (lam in 2:max_path) {
    beta_refit[[lam]] <- refit(X_F, Y_F, betamat[lam, ])
    BIC[lam] <- hdBIC(X_F, Y_F, vcov_agg, beta_refit[[lam]])
  }

  # Find the index of the minimum BIC value
  min_bic_idx <- which.min(BIC)
  betaHat <- beta_refit[[min_bic_idx]]

  # Apply adaptive LASSO if specified
  if (adaptive) {
    ada_weight <- abs(betaHat)

    # Scale the data for adaptive LASSO
    X_F_scaled <- scale(X_F, center = FALSE, scale = 1 / ada_weight)

    # Reapply LARS to the scaled data
    lars.fit <- lars(X_F_scaled, Y_F, intercept = FALSE, normalize = FALSE)
    betamat <- lars.fit$beta

    # Rescale beta estimates
    betamat_scaled <- scale(betamat, center = FALSE, scale = 1 / ada_weight)

    # Tune using BIC for the adaptive LASSO
    max_path <- k.index(betamat_scaled, n_l)
    beta_refit <- list()
    BIC <- numeric(max_path)
    BIC[1] <- hdBIC(X_F, Y_F, vcov_agg, betamat_scaled[1, ])
    beta_refit[[1]] <- betamat_scaled[1, ]

    for (lam in 2:max_path) {
      beta_refit[[lam]] <- refit(X_F, Y_F, betamat_scaled[lam, ])
      BIC[lam] <- hdBIC(X_F, Y_F, vcov_agg, beta_refit[[lam]])
    }

    # Store the best BIC, betaHat, and lambdaHat
    min_bic_idx <- which.min(BIC)
    betaHat <- beta_refit[[min_bic_idx]]
  }

  # Calculate the distribution matrix
  D <- vcov %*% t(C) %*% solve(vcov_agg)

  # Calculate the residuals
  u_l <- Y - X_l %*% betaHat

  # Generate the high-frequency series
  y <- X %*% betaHat + (D %*% u_l)

  # Return results
  output <- list('y' = y, 'betaHat' = betaHat, 'u_l' = u_l)

  return(output)
}
