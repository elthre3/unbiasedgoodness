#' Generate a design matrix and sparse coefficient vector
#'
#' @description
#' `rXb()` generates the reference high-dimensional linear regression designs
#' of Dezeure, Bühlmann, Meier and Meinshausen (2015): an \eqn{n \times p}
#' Gaussian design with a structured covariance and a coefficient vector whose
#' first `s0` entries are nonzero.
#'
#' This is a port of `rXb()` from the `hdi` package (version 0.1-10, by Ruben
#' Dezeure and Martin Maechler, GPL), which has been archived on CRAN. Given
#' the same seed it returns identical results to `hdi::rXb()`.
#'
#' @details
#' The rows of `x` are drawn from \eqn{N_p(0, \Sigma)} with [MASS::mvrnorm()],
#' where \eqn{\Sigma} depends on `xtype`:
#' * `"toeplitz"`: \eqn{\Sigma_{jk} = \rho^{|j-k|}} with \eqn{\rho} = `x.par`.
#' * `"equi.corr"`: \eqn{\Sigma_{jj} = 1} and \eqn{\Sigma_{jk} = \rho} for
#'   \eqn{j \ne k}, with \eqn{\rho} = `x.par`.
#' * `"exp.decay"`: \eqn{\Sigma^{-1}_{jk} = a^{|j-k|/b}} with
#'   `x.par = c(a, b)`.
#'
#' The first `s0` coefficients are drawn from a uniform distribution when
#' `btype` has the form `"U[lower,upper]"`, or all equal to the number `k` when
#' `btype` is `"bfixk"`; the remaining `p - s0` are zero. With
#' `permuted = TRUE` the columns of `x` are permuted, not `beta`, so the true
#' support is always `1:s0`.
#'
#' @param n Sample size, number of rows in the design matrix.
#' @param p Number of predictor variables, columns in the design matrix.
#' @param s0 Number of nonzero coefficients.
#' @param xtype Covariance structure of the design, see Details.
#' @param btype Distribution of the nonzero coefficients, see Details.
#' @param permuted Whether to randomly permute the columns of `x`.
#' @param iteration If a number, the seed is set to `iteration + 2` for the
#'   data generation and the previous random number generator state is
#'   restored afterwards. Iterations `1:50` give the data sets of Dezeure et
#'   al. (2015).
#' @param do2S Whether to compute the Toeplitz and equicorrelation covariance
#'   as `solve(solve(Sigma))`. This changes the result only through rounding
#'   and exists to reproduce the reference data sets exactly.
#' @param x.par Parameter of the design covariance, see Details.
#' @param verbose Whether to print a message when `iteration` is used.
#' @return A list with the design matrix `x` and coefficient vector `beta`.
#' @references Dezeure, R., Bühlmann, P., Meier, L. and Meinshausen, N.
#'   (2015). High-dimensional inference: confidence intervals, p-values and
#'   R-software hdi. *Statistical Science*, 30(4), 533--558.
#' @export
#' @examples
#' set.seed(1)
#' sim_data <- rXb(n = 50, p = 100, s0 = 5, xtype = "toeplitz", x.par = 1 / 3)
#' dim(sim_data$x)
#' which(sim_data$beta != 0)
rXb <- function(n, p, s0,
                xtype = c("toeplitz", "exp.decay", "equi.corr"),
                btype = "U[-2,2]", permuted = FALSE,
                iteration = NA, do2S = TRUE,
                x.par = switch(xtype,
                               "toeplitz"  = 0.9,
                               "equi.corr" = 0.8,
                               "exp.decay" = c(0.4, 5)),
                verbose = TRUE) {
  xtype <- match.arg(xtype)
  stopifnot(is.character(btype), length(btype) == 1,
            n == as.integer(n), length(n) == 1, n >= 1,
            p == as.integer(p), length(p) == 1, p >= 1,
            s0 == as.integer(s0), length(s0) == 1, 0 <= s0, s0 <= p)
  do_seed <- is.numeric(iteration) && !is.na(iteration)

  if (do_seed) {
    if (verbose) {
      message("Setting the seed to iteration + 2; the previous random ",
              "number generator state is restored afterwards.")
    }
    # As in stats:::simulate.lm; .Random.seed may not exist yet.
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      stats::runif(1)
    }
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
    seed <- iteration + 2
    set.seed(seed)
  }

  x <- rX(n = n, p = p, xtype = xtype, permuted = permuted, do2S = do2S,
          par = x.par)
  if (do_seed && s0 > 0 && !grepl("^bfix", btype)) set.seed(seed)
  beta <- rb(p = p, s0 = s0, btype = btype)
  list(x = x, beta = beta)
}

rX <- function(n, p, xtype, permuted, do2S = TRUE, par) {
  xtype <- tolower(xtype)
  stopifnot(is.numeric(par), is.finite(par))
  Sigma <- switch(xtype,
    "toeplitz" = {
      stopifnot(length(par) == 1)
      cov <- par^abs(stats::toeplitz(0:(p - 1)))
      # The double inversion only changes rounding; kept so that results
      # match the reference data sets of hdi exactly.
      if (do2S) solve(solve(cov)) else cov
    },
    "equi.corr" = {
      stopifnot(length(par) == 1)
      cov <- matrix(par, p, p)
      diag(cov) <- 1
      if (do2S) solve(solve(cov)) else cov
    },
    "exp.decay" = {
      stopifnot(length(par) == 2)
      solve(par[1]^(abs(stats::toeplitz(0:(p - 1))) / par[2]))
    },
    stop("Invalid 'xtype': must be one of 'toeplitz', 'equi.corr' or ",
         "'exp.decay'")
  )

  x <- MASS::mvrnorm(n, rep(0, p), Sigma)
  # mvrnorm() drops to a vector when n == 1
  if (!is.matrix(x)) x <- matrix(x, nrow = n, ncol = p)
  if (permuted) x[, sample.int(p), drop = FALSE] else x
}

rb <- function(p, s0, btype) {
  stopifnot(s0 <= p, p >= 0, length(s0) == 1, length(p) == 1,
            is.character(btype), length(btype) == 1)
  invalid <- paste("Invalid 'btype': use 'bfix*' for a fixed value or",
                   "'U[*,*]' for a uniform distribution with the given",
                   "lower and upper bounds.")

  if (grepl("^U\\[.*\\]$", btype) &&
        lengths(regmatches(btype, gregexpr(",", btype))) == 1) {
    bounds <- strsplit(sub("^U\\[", "", sub("\\]$", "", btype)), ",")[[1]]
    lower <- suppressWarnings(as.numeric(bounds[1]))
    upper <- suppressWarnings(as.numeric(bounds[2]))
    if (is.na(lower) || is.na(upper)) stop(invalid)
    b <- stats::runif(s0, lower, upper)
  } else if (grepl("^bfix", btype)) {
    b <- suppressWarnings(as.numeric(sub("^bfix", "", btype)))
    if (is.na(b)) stop(invalid)
    b <- rep(b, s0)
  } else {
    stop(invalid)
  }

  c(b, rep(0, p - s0))
}
