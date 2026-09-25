anova_F <- function(x, y, null, full) {
  f0 <- if (length(null)) lm(y ~ x[, null, drop = FALSE]) else lm(y ~ 1)
  a <- anova(f0, lm(y ~ x[, full, drop = FALSE]))
  list(F = a$F[2], df1 = a$Df[2], df2 = a$Res.Df[2], p = a$`Pr(>F)`[2])
}

test_that("untruncated statistic and p-value agree with anova()", {
  set.seed(3)
  n <- 100
  p <- 10
  x <- matrix(rnorm(n * p), n)
  y <- 3 * x[, 1] + rnorm(n)
  fit <- groupfs(x, y, index = 1:p, maxsteps = 4, aicstop = 0)
  for (s in 0:2) {
    null <- fit$action[seq_len(s)]
    alt <- setdiff(fit$action, null)
    out <- gof_groupfs_core(fit, null, alt, aic_steps = 0, condition = FALSE)
    ref <- anova_F(fit$x, fit$y, null, c(null, alt))
    expect_equal(out$stat, ref$F)
    expect_equal(c(out$df1, out$df2), c(ref$df1, ref$df2))
    expect_equal(out$pv, ref$p)
    expect_equal(out$pv_naive, ref$p)
  }
})

test_that("chi statistic is the norm of the added projection", {
  set.seed(4)
  n <- 60
  x <- matrix(rnorm(n * 8), n)
  y <- x[, 2] + rnorm(n)
  fit <- groupfs(x, y, index = 1:8, maxsteps = 3, sigma = 1, aicstop = 0)
  null <- fit$action[1]
  alt <- fit$action[2:3]
  out <- gof_groupfs_core(fit, null, alt, aic_steps = 0, condition = FALSE)
  rss0 <- sum(resid(lm(fit$y ~ fit$x[, null] - 1))^2)
  rss1 <- sum(resid(lm(fit$y ~ fit$x[, c(null, alt)] - 1))^2)
  expect_equal(out$type, "TC")
  expect_equal(out$stat, sqrt(rss0 - rss1))
  expect_true(is.na(out$df2))
})

test_that("gof_groupfs enlarges the selected model by m and b", {
  set.seed(1)
  n <- 100
  p <- 40
  x <- matrix(rnorm(n * p), n)
  y <- 2 * x[, 1] - 2 * x[, 2] + rnorm(n)
  out <- gof_groupfs(x, y, m = 1.5, b = 3, k = 2 * log(p))
  s <- length(out$null)
  expect_s3_class(out, "gof_groupfs")
  expect_length(out$alt, max(ceiling(1.5 * s + 3), s + 1) - s)
  expect_length(intersect(out$null, out$alt), 0)
  expect_true(out$pv >= 0 && out$pv <= 1)
  expect_true(any(out$support[, 1] <= out$stat & out$stat <= out$support[, 2]))
  expect_output(print(out), "Selective goodness-of-fit")
  ref <- anova_F(x, y, out$null, c(out$null, out$alt))
  expect_equal(out$stat, ref$F)
})

test_that("a single selected group is projected out (regression test)", {
  set.seed(2)
  n <- 80
  p <- 10
  x <- matrix(rnorm(n * p), n)
  y <- 4 * x[, 5] + rnorm(n)
  out <- gof_groupfs(x, y, b = 2, k = 3 * log(p))
  expect_equal(out$null, 5)
  ref <- anova_F(x, y, out$null, c(out$null, out$alt))
  expect_equal(out$stat, ref$F)
  expect_equal(c(out$df1, out$df2), c(ref$df1, ref$df2))
})

test_that("selective p-values are uniform under the global null", {
  skip_on_cran()
  set.seed(10)
  n <- 50
  p <- 10
  pv <- replicate(200, {
    x <- matrix(rnorm(n * p), n)
    y <- rnorm(n)
    # NA when the stopping rule never stops, leaving nothing to test against
    tryCatch(gof_groupfs(x, y, m = 1, b = 2, k = 2)$pv, error = function(e) NA)
  })
  expect_lt(mean(is.na(pv)), 0.05)
  expect_gt(ks.test(pv[!is.na(pv)], "punif")$p.value, 0.001)
})

test_that("gof_groupfs validates its arguments", {
  x <- matrix(rnorm(200), 20)
  y <- rnorm(20)
  expect_error(gof_groupfs(x, y, m = 0.5), "m >= 1")
  expect_error(gof_groupfs(x, y, aicstop = 0), "aicstop")
})
