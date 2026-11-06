# Print method for object of class gevplus

Print nicely formatted output of the fit to the GEV+

## Usage

``` r
# S3 method for class 'gevplus'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- x:

  Object of class `gevplus`, fitted with GEV+.

- digits:

  Number of digits.

- ...:

  Additional parameters.

## Value

A nicely formatted output of the fitting results.

## Examples

``` r
z <- extRemes::revd(100, loc = 20, scale = 0.5, shape = -0.2)
fit <- fgevplus(z)
print(fit)
#> GEV+ fitting
#> 
#> Type: GEV+
#> Estimator: Lmoments
#> 
#> Estimated parameters:
#> location     scale     shape  
#>  19.9263    0.4182    0.0000  
```
