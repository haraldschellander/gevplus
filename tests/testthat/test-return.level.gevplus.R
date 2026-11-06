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

test_that("return.level.gevplus rejects non-gevplus objects", {
  expect_error(return.level.gevplus(list(method = "Lmoments")),
               "x must be object of class 'gevplus'")
})

test_that("return.level rejects an unknown method", {
  x <- fgevplus(data, method = "Lmoments")
  x$method <- "INVALID"
  expect_error(return.level(x), "invalid method")
})

test_that("return.level matches the GEV quantile function", {
  x <- fgevplus(data, method = "Lmoments")
  rp <- c(2, 20, 100)
  expected <- extRemes::qevd(1 - 1 / rp, x$params[1], x$params[2], x$params[3])
  expect_equal(return.level(x, return.period = rp), expected)
  expect_true(all(diff(return.level(x, return.period = rp)) > 0))
})

test_that("bootstrap CIs for a Gumbel fit bracket the estimate", {
  set.seed(1)
  x <- fgevplus(rnorm(100, mean = 20, sd = 5), method = "Lmoments")
  expect_equal(unname(x$params[3]), 0)
  res <- return.level(x, return.period = c(2, 20, 100), do.ci = TRUE, R = 50)
  expect_named(res, c("ci_l", "rl", "ci_u"))
  expect_true(all(res$ci_l <= res$rl & res$rl <= res$ci_u))
})

test_that("bootstrap CIs get narrower with larger alpha", {
  x <- fgevplus(data, method = "Lmoments")
  set.seed(2)
  wide <- return.level(x, return.period = 100, do.ci = TRUE, R = 50, alpha = 0.05)
  set.seed(2)
  narrow <- return.level(x, return.period = 100, do.ci = TRUE, R = 50, alpha = 0.5)
  expect_lt(narrow$ci_u - narrow$ci_l, wide$ci_u - wide$ci_l)
})

test_that("return.level works for a Gumbel MLE fit", {
  set.seed(3)
  x <- fgevplus(rnorm(100, mean = 20, sd = 5), method = "MLE")
  expect_s3_class(x$fit, "fevd")
  expect_length(return.level(x, return.period = c(2, 20)), 2)
  res <- return.level(x, return.period = c(2, 20), do.ci = TRUE)
  expect_true(all(res$ci_l <= res$rl & res$rl <= res$ci_u))
})
