#' Simulation Diagnostics
#'
#' This function provides diagnostics for evaluating the accuracy of simulated data. Specifically, it computes the Mean Squared Error (MSE) between the true and estimated response vectors, and optionally, the sign recovery percentage of the coefficient vector.
#'
#' The function takes in the generated high-frequency data (`data_True`) and the estimated high-frequency data (`data_Hat`), and returns the Mean Squared Error (MSE) between the true and estimated values of the response vector. If the `sgn` parameter is set to `TRUE`, the function additionally computes the percentage of correctly recovered signs of the coefficient vector.
#'
#' @param data_Hat List containing the estimated high-frequency data, with components `y_Est` (estimated response vector) and `beta_Est` (estimated coefficient vector).
#' @param data_True List containing the true high-frequency data, with components `y_Gen` (true response vector) and `Beta_Gen` (true coefficient vector).
#' @param sgn Logical value indicating whether to compute the sign recovery percentage. Default is `FALSE`.
#'
#' @return If `sgn` is `FALSE`, the function returns the Mean Squared Error (MSE) between the true and estimated response vectors. If `sgn` is `TRUE`, the function returns a list containing both the MSE and the sign recovery percentage.
#'
#' @examples
#' true_data <- list(y_Gen = c(1, 2, 3), Beta_Gen = c(1, -1, 0))
#' est_data <- list(y_Est = c(1.1, 1.9, 2.8), beta_Est = c(1, 1, 0))
#' mse <- simulDiagnosis(est_data, true_data)
#' results <- simulDiagnosis(est_data, true_data, sgn = TRUE)
#'
#' @export
#' @keywords simulation diagnostics MSE sign recovery
simulDiagnosis <- function(data_Hat, data_True, sgn = FALSE) {

  # Extract the true and estimated response vectors
  y_True <- data_True$y_Gen
  y_Hat  <- data_Hat$y_Est

  # Ensure that the lengths of the true and estimated response vectors are the same
  if (length(y_True) != length(y_Hat)) {
    stop("The lengths of 'y_Gen' (true response vector) and 'y_Est' (estimated response vector) differ!")
  }

  # Compute the Mean Squared Error (MSE)
  MSE <- mean((y_True - y_Hat)^2)

  if (sgn == FALSE) {
    return(MSE)
  } else {

    # Extract the true and estimated coefficient vectors
    beta_True <- data_True$Beta_Gen
    beta_Hat  <- data_Hat$beta_Est

    # Ensure that the lengths of the true and estimated coefficient vectors are the same
    if (length(beta_True) != length(beta_Hat)) {
      stop("The lengths of 'Beta_Gen' (true coefficient vector) and 'beta_Est' (estimated coefficient vector) differ!")
    }

    # Calculate the signs of the true and estimated coefficients
    beta_True_signs <- sign(beta_True)
    beta_Hat_signs  <- sign(beta_Hat)

    # Calculate the correct sign recovery percentage
    correct_signs <- sum(beta_True_signs == beta_Hat_signs)
    total_signs <- length(beta_True_signs)

    sign_recovery_percentage <- (correct_signs / total_signs) * 100

    # Return both MSE and sign recovery percentage
    return(list(MSE = MSE, Sign_Recovery_Percentage = sign_recovery_percentage))
  }
}
