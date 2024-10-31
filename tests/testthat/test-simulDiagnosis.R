# Test for basic MSE calculation without sign recovery
test_that("simulDiagnosis calculates MSE correctly", {
  true_data <- list(y_Gen = c(1, 2, 3), Beta_Gen = c(1, -1, 0))
  est_data <- list(y_Est = c(1.1, 1.9, 2.8), beta_Est = c(1, 1, 0))

  # Expected MSE
  expected_mse <- mean((c(1, 2, 3) - c(1.1, 1.9, 2.8))^2)

  # Run the function
  mse <- simulDiagnosis(est_data, true_data)

  expect_equal(mse, expected_mse)
})

# Test for sign recovery functionality
test_that("simulDiagnosis calculates MSE and sign recovery correctly", {
  true_data <- list(y_Gen = c(1, 2, 3), Beta_Gen = c(1, -1, 0))
  est_data <- list(y_Est = c(1.1, 1.9, 2.8), beta_Est = c(1, 1, 0))

  # Expected MSE
  expected_mse <- mean((c(1, 2, 3) - c(1.1, 1.9, 2.8))^2)

  # Expected sign recovery percentage
  true_signs <- sign(c(1, -1, 0))
  est_signs <- sign(c(1, 1, 0))
  expected_sign_recovery <- sum(true_signs == est_signs) / length(true_signs) * 100

  # Run the function
  result <- simulDiagnosis(est_data, true_data, sgn = TRUE)

  expect_equal(result$MSE, expected_mse)
  expect_equal(result$Sign_Recovery_Percentage, expected_sign_recovery)
})


# Test for mismatched lengths of y_Gen and y_Est without checking the error message
test_that("simulDiagnosis throws an error when y_Gen and y_Est lengths differ", {
  true_data <- list(y_Gen = c(1, 2, 3), Beta_Gen = c(1, -1, 0))
  est_data <- list(y_Est = c(1.1, 1.9), beta_Est = c(1, 1, 0))  # Different length for y_Est

  # Check if an error is thrown, regardless of the exact message
  expect_error(simulDiagnosis(est_data, true_data))
})


# Test for mismatched lengths of Beta_Gen and beta_Est without checking the error message
test_that("simulDiagnosis throws an error when Beta_Gen and beta_Est lengths differ", {
  true_data <- list(y_Gen = c(1, 2, 3), Beta_Gen = c(1, -1, 0))
  est_data <- list(y_Est = c(1.1, 1.9, 2.8), beta_Est = c(1, 1))  # Different length for beta_Est

  # Check if an error is thrown when Beta_Gen and beta_Est lengths differ
  expect_error(simulDiagnosis(est_data, true_data, sgn = TRUE))
})


# Test with correct data lengths and no sign recovery
test_that("simulDiagnosis works with correct data lengths and no sign recovery", {
  true_data <- list(y_Gen = c(1, 2, 3), Beta_Gen = c(1, -1, 0))
  est_data <- list(y_Est = c(1, 2, 3), beta_Est = c(1, -1, 0))

  result <- simulDiagnosis(est_data, true_data, sgn = FALSE)

  expect_equal(result, 0)  # MSE should be zero since y_Est and y_Gen are identical
})


