test_that("simulate_gof_instance passes m and b to the test", {
  set.seed(1)
  sizes <- sapply(c(2, 6), function(b) {
    out <- simulate_gof_instance(m = 1, b = b, n = 100, p = 30, s0 = 3,
                                 method = "groupfs", k = 2 * log(30))
    c(out$s_null, out$s_alt)
  })
  expect_equal(sizes[2, ], pmax(ceiling(sizes[1, ] + c(2, 6)), sizes[1, ] + 1) - sizes[1, ])
})

test_that("simulate_gof_instance reports the true support correctly", {
  set.seed(2)
  out <- simulate_gof_instance(m = 1, b = 3, n = 100, p = 30, s0 = 2,
                               method = "groupfs", xtype = "equi.corr",
                               btype = "bfix5", k = 2 * log(30))
  expect_equal(out$null_intersection, 2)
  expect_true(out$null)
  expect_true(all(c("pval_selective", "pval_naive") %in% names(out)))
})

test_that("simulate_gof is reproducible with a seed", {
  run <- function() simulate_gof(3, m = 1, b = 2, n = 60, p = 20, s0 = 2,
                                 method = "groupfs", k = 2 * log(20), seed = 5)
  a <- run()
  b <- run()
  expect_equal(nrow(a), 3)
  expect_identical(a, b)
})

test_that("simulate_hdr passes its arguments through", {
  out <- simulate_hdr(2, n = 40, p = 30, s0 = 4, xtype = "equi.corr",
                      btype = "bfix1", seed = 1)
  expect_length(out, 2)
  expect_equal(out[[1]]$true_beta, c(rep(1, 4), rep(0, 26)))
  expect_equal(length(out[[1]]$beta_hat), 31)
})
