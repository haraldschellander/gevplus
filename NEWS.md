# gevplus 1.0.1

## Bug fixes

* Bootstrap confidence intervals in `return.level()` (with `do.ci = TRUE`) and
  `plot()` (with `ci = TRUE`) now respect the GEV+ constraint. Each bootstrap
  replicate is refitted with `fgevplus()` instead of an unconstrained GEV, so
  replicates can no longer have a negative shape parameter, and the switch to
  a Gumbel fit uses the same L-skewness threshold as the original fit. This
  affects fits with `method = "Lmoments"`. For data where the fitted shape is
  0, lower confidence bounds of long return periods are noticeably higher
  (previously too wide). MLE fits are unchanged.

## Internal

* The bootstrap code shared by `return.level()` and `plot()` now lives in a
  single internal helper.
* Test coverage increased to 100%, including new tests for `plot()`.
* Documentation website built with pkgdown and published at
  <https://haraldschellander.github.io/gevplus/>, with test coverage reported to
  Codecov.

# gevplus 1.0.0

* First release: `fgevplus()` fits the GEV distribution constrained to
  non-negative shape parameters (Moccia et al., 2021) using L-moments or MLE,
  with `return.level()`, `print()` and `plot()` methods.
