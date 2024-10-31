#' Index of support for LARS algorithm in high-dimensional settings
#' 
#' This function returns the index where the support of beta coefficients exceeds \eqn{n_l/2}, 
#' preventing the BIC from becoming erratic in high-dimensional scenarios.
#' 
#' @param coef_matrix A matrix of beta coefficients, where rows represent different models.
#' @param n_l         The length of the low-frequency response series.
#' @return            The index where the support of beta exceeds \eqn{n_l/2}, or the number of rows of the matrix if no such index is found.
#' @keywords internal 
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm

k.index <- function(coef_matrix, n_l) {
  
  # Error handling: Ensure that coef_matrix is a matrix and n_l is positive
  if (!is.matrix(coef_matrix)) stop("'coef_matrix' must be a matrix.")
  if (!is.numeric(n_l) || n_l <= 0) stop("'n_l' must be a positive numeric value.")
  
  # Calculate the count of non-zero coefficients in each row
  count <- apply(coef_matrix, 1, function(x) {sum(x != 0)})
  
  # Determine the index where the count exceeds n_l/2
  if (max(count) > n_l / 2) {
    kindex <- min(which(count > n_l / 2))
  } else {
    kindex <- nrow(coef_matrix)
  }
  
  return(kindex)
}
