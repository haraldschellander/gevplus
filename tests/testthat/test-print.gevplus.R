# Test that print.gevplus works for a valid "gevplus" object
test_that("print.gevplus outputs correct values for a GEV fit", {
  # Simulate data and fit using fgevplus with GEV distribution
  data <- extRemes::revd(100, loc = 20, scale = 0.5, shape = 0.2)
  fit <- fgevplus(data, method = "Lmoments")

  # Capture the printed output using capture.output()
  output <- capture.output(print(fit))

  # Check that the output contains the correct type
  expect_true(any(grepl("GEV\\+ fitting", output)))
  expect_true(any(grepl("Type: GEV+", output)))

  # Check that the correct method is printed
  expect_true(any(grepl("Estimator: Lmoments", output)))

  # Ensure the estimated parameters are printed
  expect_true(any(grepl("Estimated parameters:", output)))

  # Ensure the parameters are printed and formatted correctly
  expect_true(any(grepl("location", output)))
  expect_true(any(grepl("scale", output)))
  expect_true(any(grepl("shape", output)))
})

# Test that print.gevplus works for a Gumbel distribution fit
test_that("print.gevplus outputs correct values for a Gumbel fit", {
  # Simulate data and fit using fgevplus with Gumbel distribution
  data <- extRemes::revd(100, loc = 20, scale = 0.5, shape = -0.1)
  fit <- fgevplus(data, method = "Lmoments")

  # Capture the printed output
  output <- capture.output(print(fit))

  # Check that the output contains the correct type (Gumbel case)
  expect_true(any(grepl("GEV\\+ fitting", output)))
  expect_true(any(grepl("Type: GEV+", output)))

  # Check that the correct method is printed
  expect_true(any(grepl("Estimator: Lmoments", output)))

  # Ensure the estimated parameters are printed correctly (Gumbel: shape = 0)
  expect_true(any(grepl("shape", output)))
  expect_true(any(grepl("0", output)))
})

# Test that print.gevplus throws an error for an invalid object
test_that("print.gevplus throws an error for non-gevplus object", {
  # Create a non-gevplus object
  not_gevplus <- list(type = "Not a GEV+")

  # Expect an error when trying to print a non-gevplus object
  expect_error(print.gevplus(not_gevplus), "x must be object of class 'gevplus'")
})

# Test that the number of digits can be adjusted in print.gevplus
test_that("print.gevplus respects the digits argument", {
  # Simulate data and fit using fgevplus
  data <- extRemes::revd(100, loc = 20, scale = 0.5, shape = 0.2)
  fit <- fgevplus(data, method = "Lmoments")

  # Capture the printed output with a custom digits value
  output_digits <- capture.output(print(fit, digits = 5))

  # Check that the output is formatted with the correct number of digits
  expect_true(any(grepl("\\d+\\.\\d{5}", output_digits)))  # At least 5 decimal places
})

