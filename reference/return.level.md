# Return levels and confidence intervals for the GEV+

This function computes return levels and confidence intervals for the
GEV+ distribution. It handles both cases when the input fit object is of
class 'fevd' and when it's not. Confidence intervals can be constructed
with a parametric bootstrap approach.

## Usage

``` r
return.level(x, return.period = c(2, 20, 100), ...)
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
