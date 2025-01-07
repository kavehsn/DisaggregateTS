# Test for basic functionality with default parameters
test_that("TempDisaggDGP works with default parameters", {
  data <- TempDisaggDGP(n_l = 25, n = 100, aggRatio = 4, p = 10, method = 'Chow-Lin')

  expect_true(is.list(data))  # Check if the output is a list
  expect_true(all(c("y_Gen", "Y_Gen", "X_Gen", "Beta_Gen", "e_Gen") %in% names(data)))  # Ensure the output contains all required components
  expect_equal(nrow(data$y_Gen), 100)  # High-frequency series should have 100 rows
  expect_equal(nrow(data$Y_Gen), 25)   # Low-frequency series should have 25 rows
  expect_equal(ncol(data$X_Gen), 10)   # The indicator matrix should have 10 columns
})

# Test for Chow-Lin method with valid rho
test_that("TempDisaggDGP works with Chow-Lin method and valid rho", {
  data <- TempDisaggDGP(n_l = 25, n = 100, aggRatio = 4, p = 5, rho = 0.5, method = 'Chow-Lin')

  expect_true(is.list(data))
  expect_equal(nrow(data$y_Gen), 100)
  expect_equal(nrow(data$Y_Gen), 25)
  expect_equal(ncol(data$X_Gen), 5)  # 5 indicator series
  expect_equal(length(data$e_Gen), 100)  # Residuals should have 100 rows
})

# Test for input validation: rho out of bounds
test_that("TempDisaggDGP stops with invalid rho for Chow-Lin method", {
  expect_error(TempDisaggDGP(n_l = 25, n = 100, rho = 1.5, method = 'Chow-Lin'),
               "For the Chow-Lin method 'rho' must be between -1 and 1.")
})

# Test for valid Denton-Cholette method
test_that("TempDisaggDGP works with Denton-Cholette method", {
  data <- TempDisaggDGP(n_l = 25, n = 100, aggRatio = 4, p = 1, method = 'Denton-Cholette')

  expect_true(is.list(data))
  expect_equal(nrow(data$y_Gen), 100)
  expect_equal(nrow(data$Y_Gen), 25)
  expect_equal(ncol(data$X_Gen), 1)  # Denton-Cholette requires p = 1
})

# Test for invalid p value with Denton-Cholette method
test_that("TempDisaggDGP stops with invalid p for Denton-Cholette method", {
  expect_error(TempDisaggDGP(n_l = 25, n = 100, aggRatio = 4, p = 2, method = 'Denton-Cholette'),
               "For the Denton-Cholette method, p must be 1.")
})

# Test for aggregation matrices (sum, average, first, last)
test_that("TempDisaggDGP works with different aggMat options", {
  # Test with 'sum'
  data_sum <- TempDisaggDGP(n_l = 25, n = 100, aggMat = 'sum', method = 'Chow-Lin')
  expect_equal(nrow(data_sum$Y_Gen), 25)

  # Test with 'average'
  data_avg <- TempDisaggDGP(n_l = 25, n = 100, aggMat = 'average', method = 'Chow-Lin')
  expect_equal(nrow(data_avg$Y_Gen), 25)

  # Test with 'first'
  data_first <- TempDisaggDGP(n_l = 25, n = 100, aggMat = 'first', method = 'Chow-Lin')
  expect_equal(nrow(data_first$Y_Gen), 25)

  # Test with 'last'
  data_last <- TempDisaggDGP(n_l = 25, n = 100, aggMat = 'last', method = 'Chow-Lin')
  expect_equal(nrow(data_last$Y_Gen), 25)
})

# Test for Fernandez method
test_that("TempDisaggDGP works with Fernandez method", {
  data <- TempDisaggDGP(n_l = 25, n = 100, aggRatio = 4, p = 5, method = 'Fernandez')

  expect_true(is.list(data))
  expect_equal(nrow(data$y_Gen), 100)
  expect_equal(nrow(data$Y_Gen), 25)
  expect_equal(ncol(data$X_Gen), 5)
})

# Test for Litterman method
test_that("TempDisaggDGP works with Litterman method", {
  data <- TempDisaggDGP(n_l = 25, n = 100, aggRatio = 4, p = 5, method = 'Litterman')

  expect_true(is.list(data))
  expect_equal(nrow(data$y_Gen), 100)
  expect_equal(nrow(data$Y_Gen), 25)
  expect_equal(ncol(data$X_Gen), 5)
})

# Test for input validation: invalid sparsity
test_that("TempDisaggDGP stops with invalid sparsity", {
  expect_error(TempDisaggDGP(n_l = 25, n = 100, p = 5, sparsity = 1.2),
               "The 'sparsity' input can only take values in (0,1].", fixed = TRUE)
})

# Test for sparse coefficient vector with random sparsity
test_that("TempDisaggDGP works with sparse coefficient vector", {
  data <- TempDisaggDGP(n_l = 25, n = 100, p = 5, sparsity = 0.6, sparse_option = "random")

  expect_true(is.list(data))
  expect_equal(nrow(data$X_Gen), 100)
  expect_equal(length(data$Beta_Gen), 5)  # Should have 5 coefficients
  expect_true(sum(data$Beta_Gen == 0) >= 1)  # Should contain some zero coefficients due to sparsity
})
