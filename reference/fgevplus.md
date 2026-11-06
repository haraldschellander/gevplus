# Generalized extreme value distribution constrained to positive shape parameters (GEV+)

Negative GEV shape parameters are considered artifacts from sample
uncertainties. The GEV+ provides a fit to data series with shape
parameters \>= 0. If the L-skewness smaller or equal to ln(9/8)/ln2, a
Gumbel distribution is fit, otherwise a GEV. Adapted from
10.1016/j.ejrh.2021.100906

## Usage

``` r
fgevplus(x, method = c("Lmoments", "MLE"))
```

## Arguments

- x:

  A vector of block maxima.

- method:

  A character string defining the method used for parameter estimation.
  Defaults to `Lmoments`.

## Value

A list with components:

- type:

  The type of the distribution.

- method:

  Method used to fit the distribution.

- params:

  A named vector with the fitted parameters loc. sacle and shape.

- x:

  The data the distribution was fitted to.

- fit:

  The fitting object, which is either NULL or of class `fevd`.

## References

    Moccia, B. et al. (2021) 'Spatial variability of precipitation extremes over Italy using a fine-resolution gridded product', Journal of Hydrology: Regional Studies. Elsevier Ltd, 37 (October). doi: 10.1016/j.ejrh.2021.100906.

## Examples

``` r
z <- extRemes::revd(100, loc = 20, scale = 0.5, shape = -0.2)
fit <- fgevplus(z)
fit$params[3]
#> shape 
#>     0 

z <- extRemes::revd(100, loc = 20, scale = 0.5, shape = 0.2)
fit <- fgevplus(z)
fit$params[3]
#>      shape 
#> 0.03285595 
```
