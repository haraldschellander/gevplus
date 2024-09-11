test_that("fgevplus handles method argument correctly", {
  data <- rnorm(100)

  # Test with default method (Lmoments)
  result <- fgevplus(data)
  expect_equal(result$method, "Lmoments")

  # Test with explicit method (MLE)
  result <- fgevplus(data, method = "MLE")
  expect_equal(result$method, "MLE")

  # Test for invalid method
  expect_error(fgevplus(data, method = "INVALID"), "'arg' should be one of \"Lmoments\", \"MLE\"")
})

test_that("fgevplus uses L-moments method when tau3 <= tau3_thresh", {
  # Mocking Lmoments to produce tau3 <= tau3_thresh
  mock_lmoments <- list(lambdas = c(10, 2), ratios = c(NA, NA, 0.05))
  mock_data <- rnorm(100)
  assign("Lmoments", function(.data, returnobject) mock_lmoments, envir = .GlobalEnv)

  result <- fgevplus(mock_data, method = "Lmoments")

  # Check that the shape parameter is 0 (Gumbel)
  expect_equal(as.numeric(result$params[3]), 0)
  expect_equal(result$type, "GEV+")
  expect_null(result$fit)
})

test_that("fgevplus uses MLE method when tau3 <= tau3_thresh", {
  # Mocking Lmoments to produce tau3 <= tau3_thresh
  mock_lmoments <- list(lambdas = c(10, 2), ratios = c(NA, NA, 0.05))
  mock_data <- rnorm(100)
  assign("Lmoments", function(.data, returnobject) mock_lmoments, envir = .GlobalEnv)

  # Check for MLE method
  result <- fgevplus(mock_data, method = "MLE")

  expect_equal(as.numeric(result$params[3]), 0)
  expect_equal(result$type, "GEV+")
  expect_s3_class(result$fit, "fevd") # Check if fit was generated using fevd
})

test_that("fgevplus uses GEV method when tau3 > tau3_thresh", {

  mock_data <- extRemes::revd(100, loc = 20, scale = 0.5, shape = 0.2)
  result <- fgevplus(mock_data, method = "Lmoments")
  expect_true(result$params[3] > 0)
  expect_s3_class(result$fit, "fevd")
})

test_that("fgevplus returns correct output structure", {
  mock_data <- rnorm(100)
  result <- fgevplus(mock_data)

  # Check if the output is a list with the correct components
  expect_named(result, c("type", "params", "method", "x", "fit"))
  expect_named(result$params, c("loc", "scale", "shape"))
})
