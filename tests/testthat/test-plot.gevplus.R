set.seed(123)
gev_data <- extRemes::revd(50, loc = 20, scale = 5, shape = 0.2)
gumbel_data <- rnorm(50, mean = 20, sd = 5)
rp <- c(2, 10, 50)

# Draw to a null device so tests produce no files or windows
with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  force(code)
}

test_that("plot.gevplus draws all panels for each method and branch", {
  fits <- list(
    fgevplus(gev_data, method = "Lmoments"),
    fgevplus(gev_data, method = "MLE"),
    fgevplus(gumbel_data, method = "Lmoments"),
    fgevplus(gumbel_data, method = "MLE")
  )
  for (fit in fits) {
    expect_null(with_null_device(plot(fit, q = rp)))
  }
})

test_that("plot.gevplus draws each single panel type", {
  fit <- fgevplus(gev_data, method = "Lmoments")
  for (type in c("rl", "qq", "pp", "hist")) {
    expect_null(with_null_device(plot(fit, q = rp, type = type)))
  }
  expect_error(with_null_device(plot(fit, type = "foo")), "'arg' should be one of")
})

test_that("plot.gevplus draws confidence intervals", {
  # bootstrap path (Lmoments, GEV and Gumbel branch)
  for (x in list(gev_data, gumbel_data)) {
    fit <- fgevplus(x, method = "Lmoments")
    expect_null(with_null_device(
      plot(fit, q = rp, ci = TRUE, type = "rl", R = 20, alpha = 0.1)
    ))
  }
  # extRemes path (MLE)
  fit <- fgevplus(gev_data, method = "MLE")
  expect_null(with_null_device(
    plot(fit, q = rp, ci = TRUE, type = "rl", alpha = 0.1)
  ))
})

test_that("plot.gevplus restores the panel layout after type = 'all'", {
  fit <- fgevplus(gev_data, method = "Lmoments")
  with_null_device({
    plot(fit, q = rp)
    expect_equal(par("mfrow"), c(1, 1))
  })
})

test_that("plot.gevplus rejects non-gevplus objects", {
  expect_error(plot.gevplus(list()), "x must be object of class 'gevplus'")
})
