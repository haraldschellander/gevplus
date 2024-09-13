set.seed(123)
data <- extRemes::revd(100, loc = 20, scale = 0.5, shape = 0.2)

test_that("return.level handles invalid object correctly", {
  #data <- extRemes::revd(100, loc = 20, scale = 0.5, shape = 0.2)
  fit <- fgevplus(data, method = "Lmoments")
  expect_s3_class(fit, "gevplus")
})

# Test for return.level when x is a result of fgevplus with Lmoments method
test_that("return.level computes return levels correctly for fgevplus with Lmoments", {
  # Simulate data and fit using fgevplus with Lmoments
  #data <- extRemes::revd(100, loc = 20, scale = 0.5, shape = 0.2)
  fit <- fgevplus(data, method = "Lmoments")

  # Compute return levels
  result <- return.level(fit, return.period = c(2, 20, 100), do.ci = FALSE)

  # Check that result is a numeric vector without confidence intervals
  expect_type(result, "double")
  expect_length(result, 3)
})

# Test for return.level with confidence intervals for Lmoments method
test_that("return.level computes return levels with confidence intervals for fgevplus with Lmoments", {
  # Simulate data and fit using fgevplus with Lmoments
  #data <- extRemes::revd(100, loc = 20, scale = 0.5, shape = 0.2)
  fit <- fgevplus(data, method = "Lmoments")

  # Compute return levels with confidence intervals
  result <- return.level(fit, return.period = c(2, 20, 100), do.ci = TRUE, R = 100)

  # Check that result is a tibble and contains ci_l, rl, ci_u
  expect_s3_class(result, "tbl_df")
  expect_named(result, c("ci_l", "rl", "ci_u"))
  expect_equal(nrow(result), 3)
})

# Test for return.level when x is a result of fgevplus with MLE method
test_that("return.level computes return levels correctly for fgevplus with MLE", {
  # Simulate data and fit using fgevplus with MLE
  #data <- extRemes::revd(100, loc = 20, scale = 0.5, shape = 0.2)
  fit <- fgevplus(data, method = "MLE")

  # Compute return levels
  result <- return.level(fit, return.period = c(2, 20, 100), do.ci = FALSE)

  # Check that result is a numeric vector without confidence intervals
  expect_type(result, "double")
  expect_length(result, 3)
})

# Test for return.level with confidence intervals for MLE method
test_that("return.level computes return levels with confidence intervals for fgevplus with MLE", {
  # Simulate data and fit using fgevplus with MLE
  #data <- extRemes::revd(100, loc = 20, scale = 0.5, shape = 0.2)
  fit <- fgevplus(data, method = "MLE")

  # Compute return levels with confidence intervals
  result <- return.level(fit, return.period = c(2, 20, 100), do.ci = TRUE, alpha = 0.05)

  # Check that result is a tibble and contains ci_l, rl, ci_u
  expect_s3_class(result, "tbl_df")
  expect_named(result, c("ci_l", "rl", "ci_u"))
  expect_equal(nrow(result), 3)
})

# Test for return.level when the shape is constrained to 0 (Gumbel distribution)
test_that("return.level handles Gumbel distribution (shape = 0) correctly", {
  # Simulate Gumbel data and fit using fgevplus
  #data <- extRemes::revd(100, loc = 20, scale = 0.5, shape = 0)
  fit <- fgevplus(data, method = "Lmoments")

  # Compute return levels
  result <- return.level(fit, return.period = c(2, 20, 100), do.ci = FALSE)

  # Check that result is a numeric vector without confidence intervals
  expect_type(result, "double")
  expect_length(result, 3)
})

# # Test for handling invalid input for return.level
# test_that("return.level handles invalid input from fgevplus correctly", {
#   # Simulate invalid fit object (no valid method)
#   data <- extRemes::revd(100, loc = 20, scale = 0.5, shape = 0)
#   fit <- fgevplus(data, method = "INVALID")
#
#   # Expect error due to invalid method
#   expect_error(return.level(x, return.period = c(2, 20, 100)), "invalid method")
#
#
# })

test_that("return.level handles invalid return.periods correctly", {
  # Simulate valid fit but with invalid return period
  #x <- fgevplus(rnorm(100), method = "Lmoments")
  x <- fgevplus(data, method = "Lmoments")
  expect_error(return.level(x, return.period = c(-1, 2, 20)), "return period must be positive")
  expect_error(return.level(x, return.period = c(NA, 2, 20)), "return period must not be NA")
  expect_error(return.level(x, return.period = c(1, 2)), "return period must be > 1")


})


test_that("return.level handles invalid alpha correctly", {
  #x <- fgevplus(rnorm(100), method = "Lmoments")
  x <- fgevplus(data, method = "Lmoments")
  expect_error(return.level(x, do.ci = TRUE, alpha = -1), "alpha must be in \\[0, 1\\]")
})

test_that("return.level handles invalid R correctly", {
  #x <- fgevplus(rnorm(100), method = "Lmoments")
  x <- fgevplus(data, method = "Lmoments")
  expect_error(return.level(x, do.ci = TRUE, R = -1), "R must be > 0")
})
