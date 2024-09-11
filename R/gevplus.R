# package gevplus

#' @importFrom stats setNames quantile
#' @importFrom dplyr tibble
#' @importFrom Lmoments Lmoments
#' @importFrom extRemes fevd strip revd return.level qevd
NULL


#' GEV+: GEV constrained to positive shape parameters
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
  list(type = "GEV+",
       params = setNames(c(params$loc, params$scale, params$shape),
                         c("loc", "scale", "shape")),
       method = method,
       x = x,
       fit = fit)
}


#' Compute return levels and confidence intervals for
#' Generalized Extreme Value distribution with restriction to positive
#' shape parameters (GEV+) aka the gevplus distribution
#'
#' This function computes return levels and confidence intervals
#' for the GEV+ distribution. It handles both cases when the input fit object
#' is of class 'fevd' and when it's not. Confidence intervals can be constructed
#' with a parametric bootstrap approach.
#'
#' @param x A list as returned from function \code{\link{fgevplus}}.
#' @param return.period A numeric vector of return periods.
#' @param do.ci When \code{do.ci=TRUE}, either a parametric bootstrap method is used
#' to calculate confidence intervals of level alpha, or the intrinsic ci functions of
#' package extRemes are used (see also \code{\link[extRemes]{fevd}})
#' @param alpha A numeric value representing the significance level for
#' the confidence interval. Default is 0.05.
#' @param R An integer representing the number of bootstrap samples.
#' Default is 502.
#'
#' @return A tibble containing either only the estimates or additionally the
#' the lower and upper bounds of the confidence interval.
#'
#' @export
#'
#' @examples
#' z <- extRemes::revd(100, loc = 20, scale = 0.5, shape = -0.2)
#' fit <- fgevplus(z)
#' return.levels.gevplus(fit)
#' return.levels.gevplus(fit, do.ci = TRUE)
#' return.levels.gevplus(fit, return.period = c(2, 10, 20, 50, 100), do.ci = TRUE)
#'
#'
return.levels.gevplus <- function(x, return.period = c(2, 20, 100),
                                 do.ci = FALSE, alpha = 0.05, R = 502) {

  if (!x$method %in% c("Lmoments", "MLE")){
    stop("invalid method")
  }

  if (any(is.na(return.period))) {
    stop("return period must not be NA")
  }
  if (!all(return.period >= 0, na.rm = TRUE)) {
    stop("return period must be positive")
  }
  if (any(return.period <= 1)) {
    stop("return period must be > 1")
  }

  if (alpha < 0 | alpha > 1) {
    stop("alpha must be in [0, 1]")
  }

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
