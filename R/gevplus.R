# package gevplus

#' @importFrom stats setNames quantile na.omit
#' @importFrom graphics abline hist legend lines par points title
#' @importFrom dplyr tibble
#' @importFrom Lmoments Lmoments
#' @importFrom extRemes fevd strip revd return.level qevd
NULL


#' Generalized extreme value distribution constrained to positive shape parameters (GEV+)
#'
#' Negative GEV shape parameters are considered artifacts from sample
#' uncertainties. The GEV+ provides a fit to data series with
#' shape parameters >= 0. If the L-skewness smaller or equal to
#' ln(9/8)/ln2, a Gumbel distribution is fit, otherwise a GEV.
#' Adapted from 10.1016/j.ejrh.2021.100906
#'
#' @param x A vector of block maxima.
#' @param method A character string defining the method
#' used for parameter estimation. Defaults to \code{Lmoments}.
#'
#' @return A list with components:
#' \item{type}{ The type of the distribution.}
#' \item{method}{ Method used to fit the distribution.}
#' \item{params}{ A named vector with the fitted parameters loc. sacle and shape.}
#' \item{x}{ The data the distribution was fitted to.}
#' \item{fit}{ The fitting object, which is either NULL or of class \code{fevd}.}
#'
#' @export
#'
#' @references 	Moccia, B. et al. (2021) 'Spatial variability of precipitation extremes over Italy using a fine-resolution gridded product', Journal of Hydrology: Regional Studies. Elsevier Ltd, 37 (October). doi: 10.1016/j.ejrh.2021.100906.
#' @examples
#' z <- extRemes::revd(100, loc = 20, scale = 0.5, shape = -0.2)
#' fit <- fgevplus(z)
#' fit$params[3]
#'
#' z <- extRemes::revd(100, loc = 20, scale = 0.5, shape = 0.2)
#' fit <- fgevplus(z)
#' fit$params[3]
#'
fgevplus <- function(x, method = c("Lmoments", "MLE")) {

  method <- match.arg(method, several.ok = FALSE)

  tau3_thresh <- log(9/8)/log(2)

  gamma <- -digamma(1) # Euler-Mascheroni constant
  lmom <-  Lmoments::Lmoments(x, returnobject = TRUE)
  tau3 <- lmom$ratios[3]

  # Gumbel when L-Skewness <= tau3_thresh
  if (tau3 <= tau3_thresh) {
    switch(method,
           "Lmoments" = {
             lambda1 <- lmom$lambdas[1]
             lambda2 <- lmom$lambdas[2]
             scale <- lambda2 / log(2)
             loc <- lambda1 - gamma * scale
             fit <- NULL
           },
           "MLE" = {
             fit <- extRemes::fevd(x, type = "Gumbel", method = "MLE")
             params <- extRemes::strip(fit)
             scale = params[2]
             loc = params[1]
           })
    params <- dplyr::tibble(loc = loc, scale = scale, shape = 0)
  } else {
    fit <- extRemes::fevd(x, type = "GEV", method = method)
    params <- extRemes::strip(fit)
    params <- dplyr::tibble(loc = params[1], scale = params[2], shape = params[3])
  }
  res <- list(type = "GEV+",
       params = setNames(c(params$loc, params$scale, params$shape),
                         c("loc", "scale", "shape")),
       method = method,
       x = x,
       fit = fit)
  class(res) <- "gevplus"
  res
}


#' Return levels and confidence intervals for the GEV+
#'
#' This function computes return levels and confidence intervals
#' for the GEV+ distribution. It handles both cases when the input fit object
#' is of class 'fevd' and when it's not. Confidence intervals can be constructed
#' with a parametric bootstrap approach.
#'
#' @param x A list as returned from function \code{\link{fgevplus}},
#' i.e. an object of class \code{gevplus}.
#' @param return.period A numeric vector of return periods.
#' @param ... Additional parameters alpha and R for the bootstrap procedure.
#' @export
return.level <- function(x, return.period = c(2, 20, 100), ...) {
  UseMethod("return.level", x)
}


#' Return levels and confidence intervals for the GEV+
#'
#' This function computes return levels and confidence intervals
#' for the GEV+ distribution. It handles both cases when the input fit object
#' is of class 'fevd' and when it's not. Confidence intervals can be constructed
#' with a parametric bootstrap approach.
#'
#' @param x A list as returned from function \code{\link{fgevplus}},
#' i.e. an object of class \code{gevplus}.
#' @param return.period A numeric vector of return periods.
#' @param ... Additional parameters alpha and R for the bootstrap procedure.
#' @param do.ci When \code{do.ci=TRUE}, either a parametric bootstrap method is used
#' to calculate confidence intervals of level alpha, or the intrinsic ci functions of
#' package extRemes are used (see also \code{\link[extRemes]{fevd}})
#'
#' @return A tibble containing either only the estimates or additionally the
#' the lower and upper bounds of the confidence interval.
#'
#' @method return.level gevplus
#' @rdname return.level.gevplus
#' @export
#'
#' @examples
#' z <- extRemes::revd(100, loc = 20, scale = 0.5, shape = -0.2)
#' fit <- fgevplus(z)
#' return.level(fit)
#' return.level(fit, do.ci = TRUE)
#' return.level(fit, return.period = c(2, 10, 20, 50, 100), do.ci = TRUE)
#'
#'
return.level.gevplus <- function(x, return.period = c(2, 20, 100), ...,
                                 do.ci = FALSE) {
  #UseMethod("return.level")

  if(!inherits(x, "gevplus"))
    stop("x must be object of class 'gevplus'")

  if (!x$method %in% c("Lmoments", "MLE"))
    stop("invalid method")

  if (any(is.na(return.period)))
    stop("return period must not be NA")
  if (!all(return.period >= 0, na.rm = TRUE))
    stop("return period must be positive")
  if (any(return.period <= 1))
    stop("return period must be > 1")

  opts <- list(...)
  if (is.null(opts$alpha)) {
    alpha <- 0.05
  } else {
    alpha <- opts$alpha
  }

  if (is.null(opts$R)) {
    R <- 502
  } else {
    R <- opts$R
  }

  if (alpha < 0 | alpha > 1)
    stop("alpha must be in [0, 1]")
  if (R <= 0)
    stop("R must be > 0")

  if (inherits(x$fit, "fevd") & x$method != "Lmoments") {
    rls <- extRemes::return.level(x$fit, return.period = return.period)
    res <- as.numeric(rls)

    if (do.ci) {
      rls <- extRemes::return.level(x$fit, return.period = return.period, do.ci = TRUE, alpha = alpha)
      colnames(rls) <- c("ci_l", "rl", "ci_u")
      t <- dplyr::tibble(
        ci_l = rls[, 1],
        rl = rls[, 2],
        ci_u = rls[, 3]
      )
      res <- t
    }

  } else {
    est <- extRemes::qevd(1 - 1 / return.period, x$params[1], x$params[2], x$params[3])
    res <- est

    if (do.ci) {
      # parametric bootstrap
      sam <- list()
      for (i in 1:R){
        nd <- extRemes::revd(n = length(x$x), x$params[1], x$params[2], x$params[3])
        fit <- extRemes::fevd(nd, type = "GEV", method = x$method)
        pars <- extRemes::strip(fit)
        sam[[i]] <- extRemes::qevd(1 - 1 / return.period, pars[1], pars[2], pars[3])
      }
      sam <- do.call(cbind, sam)
      cis <- apply(sam, 1, quantile, probs = c(alpha / 2, 1 - alpha / 2))
      rls <- dplyr::tibble(ci_l = cis[1, ], rl = est, ci_u = cis[2, ])
      res <- rls
    }
  }
  return(res)
}


#' Print method for object of class gevplus
#'
#' Print nicely formatted output of the fit to the GEV+
#'
#' @param x Object of class \code{gevplus}, fitted with GEV+.
#' @param digits Number of digits.
#' @param ... Additional parameters.
#'
#' @return A nicely formatted output of the fitting results.
#'
#' @export
#'
#' @examples
#' z <- extRemes::revd(100, loc = 20, scale = 0.5, shape = -0.2)
#' fit <- fgevplus(z)
#' print(fit)
print.gevplus <- function(x, digits = max(3, getOption("digits") - 3), ...){
  if(!inherits(x, "gevplus"))
    stop("x must be object of class 'gevplus'")

  cat("GEV+ fitting\n\n")
  cat(paste0("Type: ", x$type, "\n"))
  cat(paste0("Estimator: ", x$method, "\n"))
  cat("\nEstimated parameters:\n")
  params <- setNames(x$params, c("location", "scale", "shape"))
  print.default(format(params, digits = digits), print.gap = 2, quote = FALSE)

  invisible(x)
}


#' Plot graphs of GEV+ fit
#'
#' A return level plot, qq-plot, pp-plot and a histogram
#' with the fitted density is produced
#'
#' @param x An object of class\code{gevplus}
#' @param q vector of return periods, \eqn{q > 1}.
#' @param ci if \code{ci=TRUE}, confidence intervals will be computed.
#' @param type if omitted a panel with a return level plot (\code{type='rl'},
#' a density plot (\code{type='hist'}), a qq-plot (\code{type='qq'}) and a
#' probability plot (\code{tpe='pp'}) are shown.
#' @param ... Further parameters for plotting or confidence interval calculation (alpha, R)
#' may also be supplied as arguments.
#' See e.g. \link[base]{plot}.
#'
#' @method plot gevplus
#' @return No return value, only a plot is produced.
#' @export
#'
#' @examples
#' z <- extRemes::revd(100, loc = 20, scale = 10, shape = -0.2)
#' fit <- fgevplus(z)
#' plot(fit, ci = TRUE)
#'
#' # Compare standard GEV with GEV+
#' z <- rnorm(100)
#' rp <- c(2, 10, 20, 30, 50, 75, 100, 150, 200)
#' fit_gev <- extRemes::fevd(z)
#' fit_gevplus <- fgevplus(z)
#' plot(fit_gev, rperiods = rp, type = "rl", ylim = c(0, max(z) * 1.2),
#' main = paste0(fit_gev$type, " / ", fit_gev$method))
#' plot(fit_gevplus, q = rp, ci = TRUE, type = "rl")
#'
plot.gevplus <- function(x, q = c(2, 10, 20, 30, 50, 75, 100, 150, 200),
                         ci = FALSE,
                         type = c("all", "rl", "qq", "pp", "hist"), ...){

  if(!inherits(x, "gevplus"))
    stop("x must be object of class 'gevplus'")

  type <- match.arg(type)

  opts <- list(...)
  if (is.null(opts$alpha)) {
    alpha <- 0.05
  } else {
    alpha <- opts$alpha
  }

  if (is.null(opts$R)) {
    R <- 502
  } else {
    R <- opts$R
  }

  obs.y <- na.omit(x$x)
  n <- length(obs.y)
  Fi <- (1:n) / (n + 1)
  obs.x <- 1 / (1 - Fi)

  if(type == "all"){
    #oldpar <- par(no.readonly = TRUE)
    #on.exit(par(oldpar))
    op <- par()
    par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  }

  if (is.element(type, c("all", "rl"))) {
    if (inherits(x$fit, "fevd") & x$method != "Lmoments") {
      if (ci) {
        rls <- extRemes::return.level(x$fit, return.period = q, do.ci = TRUE, alpha = alpha)
        colnames(rls) <- c("ci_l", "rl", "ci_u")
        t <- dplyr::tibble(
          ci_l = rls[, 1],
          rl = rls[, 2],
          ci_u = rls[, 3]
        )
        res <- t
        plot(q, res$rl, type = "l", xlab = "Return Period (years)", ylab = "Return Level", ylim = c(0, max(res) * 1.2), log = "x")
        points(obs.x, sort(obs.y))
        lines(q, res$ci_l, lty = 2, col = "grey")
        lines(q, res$ci_u, lty = 2, col = "grey")
      } else {
        rls <- extRemes::return.level(x$fit, return.period = q)
        est <- as.numeric(rls)
        plot(q, est, type = "l", xlab = "Return Period (years)", ylab = "Return Level", ylim = c(0, max(est) * 1.2), log = "x")
        points(obs.x, sort(obs.y))
      }

    } else {
      est <- extRemes::qevd(1 - 1 / q, x$params[1], x$params[2], x$params[3])

      if (ci) {
        # parametric bootstrap
        sam <- list()
        for (i in 1:R){
          nd <- extRemes::revd(n = length(x$x), x$params[1], x$params[2], x$params[3])
          fit <- extRemes::fevd(nd, type = "GEV", method = x$method)
          pars <- extRemes::strip(fit)
          sam[[i]] <- extRemes::qevd(1 - 1 / q, pars[1], pars[2], pars[3])
        }
        sam <- do.call(cbind, sam)
        cis <- apply(sam, 1, quantile, probs = c(alpha / 2, 1 - alpha / 2))
        rls <- dplyr::tibble(ci_l = cis[1, ], rl = est, ci_u = cis[2, ])
        plot(q, rls$rl, type = "l", xlab = "Return Period (years)", ylab = "Return Level", ylim = c(0, max(rls) * 1.2), log = "x")
        points(obs.x, sort(obs.y))
        lines(q, rls$ci_l, lty = 2, col = "grey")
        lines(q, rls$ci_u, lty = 2, col = "grey")
      } else {
        plot(q, est, type = "l", xlab = "Return Period (years)", ylab = "Return Level", ylim = c(0, max(est) * 1.2), log = "x")
        points(obs.x, sort(obs.y))
      }
    }

  }

  if (is.element(type, c("all", "hist"))) {
    # density plot
    h <- hist(obs.y, plot = FALSE)
    dens.x <- seq(min(h$breaks), min(max(h$breaks)), length = 100)
    dens.y <- extRemes::devd(dens.x, loc = x$params[1], scale = x$params[2], shape = x$params[3], type = "GEV")
    hist(obs.y, prob = TRUE, main = "Observed yearly maxima", xlab = paste("N =", length(obs.y)))
    lines(dens.x, dens.y, lty = 2, col = "blue", lwd = 1.5)
    legend("topright", legend = c("Empirical", "Modeled"),
           col = c("black", "blue"), lty = c(1, 2), lwd = c(1, 1.5), bty = "n")
  }


  if (is.element(type, c("all", "qq"))) {
    # quantile plot
    q.m <- extRemes::revd(n, loc = x$params[1], scale = x$params[2], shape = x$params[3], type = "GEV")
    plot(q.m, sort(obs.y), main = "Quantile plot", xlab = "Model Quantiles", ylab = "Empirical Quantiles")
    abline(c(0,1))
  }

  if (is.element(type, c("all", "pp"))) {
    # probability plot
    p.m <- extRemes::pevd(sort(obs.y), loc = x$params[1], scale = x$params[2], shape = x$params[3], type = "GEV")
    plot(Fi, p.m, main = "Probability plot", xlab = "Model Probabilities", ylab = "Empirical Probabilities", xlim = c(0, 1), ylim = c(0, 1))
    abline(c(0,1))
  }

  title_strg <- paste(x$type, "/", x$method)
  title(title_strg, outer = TRUE)

  if (type == "all")
    par(mfrow = op$mfrow, oma = op$oma)
  invisible()

}

