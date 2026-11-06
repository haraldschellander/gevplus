# Return levels and confidence intervals for the GEV+

This function computes return levels and confidence intervals for the
GEV+ distribution. It handles both cases when the input fit object is of
class 'fevd' and when it's not. Confidence intervals can be constructed
with a parametric bootstrap approach.

## Usage

``` r
# S3 method for class 'gevplus'
return.level(x, return.period = c(2, 20, 100), ..., do.ci = FALSE)
```

## Arguments

- x:

  A list as returned from function
  [`fgevplus`](https://haraldschellander.github.io/gevplus/reference/fgevplus.md),
  i.e. an object of class `gevplus`.

- return.period:

  A numeric vector of return periods.

- ...:

  Additional parameters alpha and R for the bootstrap procedure.

- do.ci:

  When `do.ci=TRUE`, either a parametric bootstrap method is used to
  calculate confidence intervals of level alpha, or the intrinsic ci
  functions of package extRemes are used (see also
  [`fevd`](https://rdrr.io/pkg/extRemes/man/fevd.html))

## Value

A tibble containing either only the estimates or additionally the the
lower and upper bounds of the confidence interval.

## Examples

``` r
z <- extRemes::revd(100, loc = 20, scale = 0.5, shape = -0.2)
fit <- fgevplus(z)
return.level(fit)
#> [1] 20.03384 21.11975 21.79955
return.level(fit, do.ci = TRUE)
#> # A tibble: 3 × 3
#>    ci_l    rl  ci_u
#>   <dbl> <dbl> <dbl>
#> 1  19.9  20.0  20.1
#> 2  20.9  21.1  21.4
#> 3  21.5  21.8  22.5
return.level(fit, return.period = c(2, 10, 20, 50, 100), do.ci = TRUE)
#> # A tibble: 5 × 3
#>    ci_l    rl  ci_u
#>   <dbl> <dbl> <dbl>
#> 1  19.9  20.0  20.1
#> 2  20.6  20.8  21.0
#> 3  20.9  21.1  21.4
#> 4  21.2  21.5  22.0
#> 5  21.5  21.8  22.5

```
