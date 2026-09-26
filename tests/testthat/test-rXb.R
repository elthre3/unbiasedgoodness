test_that("rXb returns a design and a sparse coefficient vector", {
  set.seed(1)
  sim <- rXb(20, 10, 3, xtype = "toeplitz", x.par = 1 / 3, verbose = FALSE)
  expect_equal(dim(sim$x), c(20, 10))
  expect_length(sim$beta, 10)
  expect_true(all(sim$beta[1:3] != 0))
  expect_true(all(sim$beta[4:10] == 0))
  expect_true(all(abs(sim$beta) <= 2))
})

test_that("rXb supports all design and coefficient types", {
  for (xtype in c("toeplitz", "exp.decay", "equi.corr")) {
    sim <- rXb(15, 8, 2, xtype = xtype, btype = "bfix5", permuted = TRUE,
               verbose = FALSE)
    expect_equal(dim(sim$x), c(15, 8))
    expect_equal(sim$beta, c(5, 5, rep(0, 6)))
  }
  expect_equal(dim(rXb(1, 4, 1, verbose = FALSE)$x), c(1, 4))
  expect_error(rXb(10, 5, 1, btype = "N(0,1)", verbose = FALSE), "btype")
  expect_error(rXb(10, 5, 6, verbose = FALSE))
})

test_that("rXb with iteration is reproducible and restores the RNG state", {
  set.seed(42)
  before <- .Random.seed
  a <- rXb(10, 6, 2, iteration = 3, verbose = FALSE)
  expect_identical(.Random.seed, before)
  b <- rXb(10, 6, 2, iteration = 3, verbose = FALSE)
  expect_identical(a, b)
})
