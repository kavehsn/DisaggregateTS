#' Refit LASSO estimate into GLS 
#' 
#' This function reduces the bias in LASSO estimates by re-fitting the 
#' active set of coefficients back into GLS (Generalized Least Squares). 
#' 
#' @param X     Aggregated indicator series matrix that has been GLS rotated (an \eqn{n_l \times p} matrix).
#' @param Y     Low-frequency response vector that has been GLS rotated (an \eqn{n_l \times 1} vector).
#' @param beta  Estimated beta coefficients from the LARS algorithm (a \eqn{p \times 1} vector).
#' @return      A debiased estimate of the beta coefficients (a \eqn{p \times 1} vector).
#' @keywords internal 
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm

refit <- function(X, Y, beta) {
  
  # Error handling: Ensure that X, Y, and beta have compatible dimensions
  p <- ncol(X)
  if (length(beta) != p) stop("'beta' must have the same number of elements as columns in 'X'.")
  if (nrow(X) != length(Y)) stop("'X' and 'Y' must have the same number of rows.")
  
  # Identify active (non-zero) coefficients in beta
  active <- which(beta != 0)
  
  # If no active coefficients, return a zero vector
  if (length(active) == 0) {
    return(rep(0, p))
  }
  
  # Subset X to the active set of coefficients
  X_active <- X[, active]
  
  # Fit OLS model using the active set
  lm_fit <- lm(Y ~ 0 + X_active)
  
  # Extract the coefficients from the OLS fit
  beta_refit <- rep(0, p)
  beta_refit[active] <- lm_fit$coefficients
  
  return(beta_refit)
}
