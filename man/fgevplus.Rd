% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gevplus.R
\name{fgevplus}
\alias{fgevplus}
\title{Generalized extreme value distribution constrained to positive shape parameters (GEV+)}
\usage{
fgevplus(x, method = c("Lmoments", "MLE"))
}
\arguments{
\item{x}{A vector of block maxima.}

\item{method}{A character string defining the method
used for parameter estimation. Defaults to \code{Lmoments}.}
}
\value{
A list with components:
\item{type}{ The type of the distribution.}
\item{method}{ Method used to fit the distribution.}
\item{params}{ A named vector with the fitted parameters loc. sacle and shape.}
\item{x}{ The data the distribution was fitted to.}
\item{fit}{ The fitting object, which is either NULL or of class \code{fevd}.}
}
\description{
Negative GEV shape parameters are considered artifacts from sample
uncertainties. The GEV+ provides a fit to data series with
shape parameters >= 0. If the L-skewness smaller or equal to
ln(9/8)/ln2, a Gumbel distribution is fit, otherwise a GEV.
Adapted from 10.1016/j.ejrh.2021.100906
}
\examples{
z <- extRemes::revd(100, loc = 20, scale = 0.5, shape = -0.2)
fit <- fgevplus(z)
fit$params[3]

z <- extRemes::revd(100, loc = 20, scale = 0.5, shape = 0.2)
fit <- fgevplus(z)
fit$params[3]

}
\references{
\if{html}{\out{<div class="sourceCode">}}\preformatted{Moccia, B. et al. (2021) 'Spatial variability of precipitation extremes over Italy using a fine-resolution gridded product', Journal of Hydrology: Regional Studies. Elsevier Ltd, 37 (October). doi: 10.1016/j.ejrh.2021.100906.
}\if{html}{\out{</div>}}
}
