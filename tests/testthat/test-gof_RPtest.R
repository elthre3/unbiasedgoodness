test_that("gof_RPtest returns p-values for a disjoint larger model", {
  skip_on_cran()
  set.seed(1)
  n <- 80
  p <- 100
  sim <- hdi::rXb(n, p, 3, xtype = "equi.corr", x.par = 1/20, verbose = FALSE)
  y <- as.numeric(sim$x %*% sim$beta + rnorm(n))
  out <- suppressWarnings(gof_RPtest(sim$x, y, m = 1.2, b = 4, B = 9L))
  pvals <- out[startsWith(names(out), "pval")]
  expect_named(pvals, c("pval_OLS", "pval_lasso", "pval_lasso_beta",
                        "pval_full", "pval_full_lasso"))
  expect_true(all(unlist(pvals) >= 0 & unlist(pvals) <= 1))
  expect_length(intersect(out$selected_support, out$alt_support), 0)
  expect_gte(length(out$selected_support), 2)
  expect_gt(length(out$alt_support), 0)
})

test_that("gof_RPtest validates m and b", {
  expect_error(gof_RPtest(matrix(rnorm(100), 10), rnorm(10), b = 0), "b > 0")
})
