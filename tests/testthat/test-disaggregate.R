# Test for basic functionality with Chow-Lin method
test_that("disaggregate works correctly with Chow-Lin method", {
  Y <- matrix(rnorm(25), ncol = 1)  # Low-frequency response (25 observations)
  X <- matrix(rnorm(100 * 4), ncol = 4)  # High-frequency indicators (100 rows, 4 columns)

  result <- disaggregate(Y = Y, X = X, method = "Chow-Lin", aggRatio = 4)

  expect_true(is.list(result))  # Expect the result to be a list
  expect_true(all(c("y_Est", "beta_Est", "rho_Est", "ul_Est") %in% names(result)))
  expect_equal(length(result$y_Est), 100)  # The high-frequency result should have 100 rows
})

# Test for input validation for non-matrix inputs
test_that("disaggregate stops with non-matrix inputs", {
  Y <- rnorm(25)  # Not a matrix
  X <- rnorm(100) # Not a matrix

  expect_error(disaggregate(Y = Y, X = X), "X and Y must be a matrices!")
})

# Test for input validation with invalid method
test_that("disaggregate stops with wrong method input", {
  Y <- matrix(rnorm(25), ncol = 1)
  X <- matrix(rnorm(100 * 4), ncol = 4)

  expect_error(disaggregate(Y = Y, X = X, method = "InvalidMethod"),
               "Wrong method inputted")
})

# Test different aggMat options (sum, average, first, last)
test_that("disaggregate works with different aggMat options", {
  Y <- matrix(rnorm(25), ncol = 1)
  X <- matrix(rnorm(100 * 4), ncol = 4)

  # Test with 'sum'
  result_sum <- disaggregate(Y = Y, X = X, aggMat = "sum", method = "Chow-Lin")
  expect_equal(length(result_sum$y_Est), 100)

  # Test with 'average'
  result_avg <- disaggregate(Y = Y, X = X, aggMat = "average", method = "Chow-Lin")
  expect_equal(length(result_avg$y_Est), 100)

  # Test with 'first'
  result_first <- disaggregate(Y = Y, X = X, aggMat = "first", method = "Chow-Lin")
  expect_equal(length(result_first$y_Est), 100)

  # Test with 'last'
  result_last <- disaggregate(Y = Y, X = X, aggMat = "last", method = "Chow-Lin")
  expect_equal(length(result_last$y_Est), 100)
})

# Test for different aggRatio values
test_that("disaggregate works with different aggRatio values", {
  Y <- matrix(rnorm(25), ncol = 1)

  # Test with aggRatio = 3
  X_agg3 <- matrix(rnorm(75 * 4), ncol = 4)
  result_agg3 <- disaggregate(Y = Y, X = X_agg3, aggRatio = 3, method = "Chow-Lin")
  expect_equal(length(result_agg3$y_Est), 75)

  # Test with aggRatio = 4
  X_agg4 <- matrix(rnorm(100 * 4), ncol = 4)
  result_agg4 <- disaggregate(Y = Y, X = X_agg4, aggRatio = 4, method = "Chow-Lin")
  expect_equal(length(result_agg4$y_Est), 100)
})

# Test for Denton method with valid inputs
test_that("disaggregate works correctly with Denton method", {
  Y <- matrix(rnorm(25), ncol = 1)
  X <- matrix(rnorm(100), ncol = 1)  # High-frequency indicators with one column for Denton

  result <- disaggregate(Y = Y, X = X, method = "Denton", Denton = 'additive-first-diff')

  expect_true(is.list(result))
  expect_true(all(c("y_Est", "beta_Est", "rho_Est", "ul_Est") %in% names(result)))
  expect_equal(length(result$y_Est), 100)  # Should have 100 rows for the high-frequency output
})

# Test for Denton method with invalid number of indicators
test_that("disaggregate stops for Denton method with more than one indicator", {
  Y <- matrix(rnorm(25), ncol = 1)
  X <- matrix(rnorm(100 * 2), ncol = 2)  # Two indicators

  expect_error(disaggregate(Y = Y, X = X, method = "Denton"),
               "X has more than 1 indicator. The Denton/Denton-Cholette methods requires only one indicator.")
})

# Test for Fernandez method
test_that("disaggregate works correctly with Fernandez method", {
  Y <- matrix(rnorm(25), ncol = 1)
  X <- matrix(rnorm(100 * 4), ncol = 4)

  result <- disaggregate(Y = Y, X = X, method = "Fernandez")

  expect_true(is.list(result))
  expect_true(all(c("y_Est", "beta_Est", "rho_Est", "ul_Est") %in% names(result)))
  expect_equal(length(result$y_Est), 100)
})

# Test for rho optimization in Chow-Lin method
test_that("disaggregate optimizes rho correctly in Chow-Lin method", {
  Y <- matrix(rnorm(25), ncol = 1)
  X <- matrix(rnorm(100 * 4), ncol = 4)

  result <- disaggregate(Y = Y, X = X, method = "Chow-Lin")

  expect_true(!is.na(result$rho_Est))  # Rho should be estimated
})
