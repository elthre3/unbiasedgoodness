# unbiasedgoodness

<!-- badges: start -->
[![R-CMD-check](https://github.com/joftius/unbiasedgoodness/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/joftius/unbiasedgoodness/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Goodness-of-fit tests for linear models chosen by model selection.

Diagnostics computed on the same data used to select a model are biased: the
selected model was chosen *because* it fits those data well. This package
tests a selected model against a data-driven enlargement of it. If the
selected model has `s` variables, the alternative adds variables along the
same selection path until it has about `ceiling(m * s + b)`.

- `gof_groupfs()` selects a model by forward stepwise (optionally with groups
  of variables) and a penalized stopping rule, and computes a **selective
  F (or chi) test** that is valid conditional on the selected and enlarged
  models (Loftus and Taylor, 2015).
- `gof_RPtest()` selects a model with the cross-validated lasso and applies
  the **residual prediction tests** of Shah and Bühlmann (2018) from the
  [RPtests](https://CRAN.R-project.org/package=RPtests) package. These are not
  adjusted for selection; the package articles study their selective error
  rates by simulation.
- `simulate_gof()` estimates selective type 1 error and power of either test
  on simulated high-dimensional regressions.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("joftius/unbiasedgoodness")
```

## Example

``` r
library(unbiasedgoodness)
set.seed(1)
n <- 100
p <- 60
sim_data <- hdi::rXb(n, p, s0 = 3, xtype = "toeplitz", x.par = 1 / 3)
x <- sim_data$x
y <- as.numeric(x %*% sim_data$beta + rnorm(n))

gof_groupfs(x, y, m = 1.2, b = 5)
gof_RPtest(x, y, m = 1.2, b = 5)
```

See the [package website](http://joshualoftus.com/unbiasedgoodness/) for
articles with simulation results.

## References

- Loftus, J. R. and Taylor, J. E. (2015). Selective inference in regression
  models with groups of variables. arXiv:1511.01478.
- Shah, R. D. and Bühlmann, P. (2018). Goodness-of-fit tests for high
  dimensional linear models. *Journal of the Royal Statistical Society:
  Series B*, 80(1), 113–135.
