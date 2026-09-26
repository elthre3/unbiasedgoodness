#' Simulate goodness-of-fit tests after model selection
#'
#' @description
#' `simulate_gof_instance()` generates one high-dimensional linear regression
#' data set with [rXb()], selects and enlarges a model, and computes
#' goodness-of-fit p-values with [gof_RPtest()] or [gof_groupfs()].
#' `simulate_gof()` repeats this `niters` times, optionally in parallel, and
#' returns one row per iteration.
#'
#' @details
#' The true coefficient vector returned by [rXb()] is nonzero exactly in
#' its first `s0` entries, so the true support is `1:s0` (with
#' `permuted = TRUE` the columns of the design are permuted, not `beta`).
#' The column `null` is `TRUE` when the selected model contains the whole true
#' support, i.e. when the goodness-of-fit null hypothesis holds for the
#' selected model. Rejection rates computed within `null == TRUE` therefore
#' estimate \eqn{P(\mathrm{reject} \mid \mathrm{selected\ model\ is\ correct})},
#' a selective type I error averaged over the selected models.
#'
#' Iterations in which the test function throws an error (for example because
#' the model could not be enlarged) are dropped; their number is stored in
#' the attribute `"n_failed"` of the result of `simulate_gof()`.
#'
#' @param niters Number of simulation iterations.
#' @param m Multiplier of the selected model size.
#' @param b Additive increase in the size of the alternative model.
#' @param n Sample size, number of rows in the design matrix.
#' @param p Number of predictor variables, columns in the design matrix.
#' @param s0 Number of nonzero coefficients in the true linear model.
#' @param method Which goodness-of-fit test to use.
#' @param xtype,btype,permuted Data-generating process, passed to
#'   [rXb()].
#' @param x.par Parameter of the design correlation, passed to [rXb()].
#'   The default is `1/3` for `"toeplitz"`, `1/20` for `"equi.corr"` and
#'   `c(0.4, 5)` for `"exp.decay"`, weaker correlation than the defaults of
#'   [rXb()].
#' @param sigma Error standard deviation of the data-generating process.
#' @param cores Number of cores. Uses forked processes via
#'   [parallel::mclapply()] where available, otherwise a socket cluster.
#' @param seed Optional seed. When given, the `"L'Ecuyer-CMRG"` generator is
#'   used so results are reproducible for a fixed value of `cores`.
#' @param ... Further arguments passed to [gof_RPtest()] or [gof_groupfs()].
#' @return A data frame with one row per successful iteration, containing
#'   design diagnostics (`max_corr`, `min_beta`), selected and added model
#'   sizes (`s_null`, `s_alt`), the number of true variables in the selected
#'   and enlarged models (`null_intersection`, `alt_intersection`), whether
#'   each contains the true support (`null`, `alt`), and one column per
#'   p-value (names starting with `pval`).
#' @export
#' @examples
#' \donttest{
#' simulate_gof_instance(m = 1.2, b = 5, n = 100, p = 200, s0 = 5,
#'                       method = "groupfs")
#' }
simulate_gof_instance <- function(m, b, n, p, s0,
                                  method = c("RPtest", "groupfs"),
                                  xtype = c("toeplitz", "exp.decay", "equi.corr"),
                                  btype = "U[-2,2]", permuted = TRUE,
                                  x.par = NULL, sigma = 1, ...) {
  method <- match.arg(method)
  xtype <- match.arg(xtype)
  if (is.null(x.par)) x.par <- default_x_par(xtype)

  sim_data <- rXb(n = n, p = p, s0 = s0, xtype = xtype, btype = btype,
                       permuted = permuted, x.par = x.par, verbose = FALSE)
  x <- sim_data$x
  beta <- sim_data$beta
  y <- as.numeric(x %*% beta + sigma * stats::rnorm(n))

  true_support <- seq_len(s0)
  has_both <- s0 > 0 && s0 < p
  diagnostics <- data.frame(
    max_corr = if (has_both) max(abs(stats::cor(x[, true_support], x[, -true_support]))) else NA,
    min_beta = if (s0 > 0) min(abs(beta[true_support])) else NA
  )

  if (method == "RPtest") {
    fit <- gof_RPtest(x, y, m = m, b = b, ...)
    selected <- fit$selected_support
    added <- fit$alt_support
    pvals <- fit[startsWith(names(fit), "pval")]
  } else {
    fit <- gof_groupfs(x, y, m = m, b = b, ...)
    selected <- fit$null
    added <- fit$alt
    pvals <- list(pval_selective = fit$pv, pval_naive = fit$pv_naive)
  }

  null_intersection <- length(intersect(selected, true_support))
  alt_intersection <- length(intersect(union(selected, added), true_support))
  cbind(diagnostics,
        data.frame(s_null = length(selected),
                   s_alt = length(added),
                   null_intersection = null_intersection,
                   alt_intersection = alt_intersection,
                   null = null_intersection == s0,
                   alt = alt_intersection == s0),
        as.data.frame(pvals))
}

#' @rdname simulate_gof_instance
#' @export
simulate_gof <- function(niters, m, b, n, p, s0, ..., cores = 1L, seed = NULL) {
  one <- function(i) {
    tryCatch(simulate_gof_instance(m = m, b = b, n = n, p = p, s0 = s0, ...),
             error = function(e) NULL)
  }
  results <- run_iterations(niters, one, cores = cores, seed = seed)
  failed <- vapply(results, is.null, logical(1))
  if (any(failed)) {
    warning(sum(failed), " of ", niters, " iterations failed and were dropped")
  }
  out <- do.call(rbind, results[!failed])
  attr(out, "n_failed") <- sum(failed)
  out
}

#' Generic simulation harness for high-dimensional linear regression
#'
#' @description
#' `instance_hdr()` generates data with [rXb()], computes a response
#' with `y_fun`, fits a model with `fit_fun`, and post-processes the fit with
#' `post_fun`. `simulate_hdr()` repeats this `niters` times and returns a
#' list of the results.
#'
#' @inheritParams simulate_gof_instance
#' @param y_fun Function `(x, beta, ...)` returning the response.
#' @param y_args List of further arguments to `y_fun`.
#' @param fit_fun Function `(x, y, beta, ...)` returning a fitted object.
#' @param fit_args List of further arguments to `fit_fun`.
#' @param post_fun Function `(fit_obj, x, y, beta, ...)` returning the
#'   quantities of interest.
#' @param post_args List of further arguments to `post_fun`.
#' @return For `instance_hdr()`, the value of `post_fun`; for
#'   `simulate_hdr()`, a list of such values.
#' @export
#' @examples
#' out <- simulate_hdr(2, n = 50, p = 100, s0 = 3, seed = 1)
#' out[[1]]$lambda
simulate_hdr <- function(niters, n, p, s0,
                         xtype = c("toeplitz", "exp.decay", "equi.corr"),
                         btype = "U[-2,2]", permuted = TRUE, x.par = NULL,
                         y_fun = y_standard_linear, y_args = list(),
                         fit_fun = fit_glmnet_cv, fit_args = list(),
                         post_fun = post_glmnet_cv, post_args = list(),
                         cores = 1L, seed = NULL) {
  xtype <- match.arg(xtype)
  one <- function(i) {
    instance_hdr(n, p, s0, xtype = xtype, btype = btype, permuted = permuted,
                 x.par = x.par, y_fun = y_fun, y_args = y_args,
                 fit_fun = fit_fun, fit_args = fit_args,
                 post_fun = post_fun, post_args = post_args)
  }
  run_iterations(niters, one, cores = cores, seed = seed)
}

#' @rdname simulate_hdr
#' @export
instance_hdr <- function(n, p, s0,
                         xtype = c("toeplitz", "exp.decay", "equi.corr"),
                         btype = "U[-2,2]", permuted = TRUE, x.par = NULL,
                         y_fun = y_standard_linear, y_args = list(),
                         fit_fun = fit_glmnet_cv, fit_args = list(),
                         post_fun = post_glmnet_cv, post_args = list()) {
  xtype <- match.arg(xtype)
  if (is.null(x.par)) x.par <- default_x_par(xtype)
  sim_data <- rXb(n = n, p = p, s0 = s0, xtype = xtype, btype = btype,
                       permuted = permuted, x.par = x.par, verbose = FALSE)
  x <- sim_data$x
  beta <- sim_data$beta
  y <- do.call(y_fun, c(list(x, beta), y_args))
  fit_obj <- do.call(fit_fun, c(list(x, y, beta), fit_args))
  do.call(post_fun, c(list(fit_obj, x, y, beta), post_args))
}

#' @describeIn simulate_hdr Linear model response with Gaussian errors of
#'   standard deviation `sigma`.
#' @param x,beta,y,fit_obj Design, coefficients, response and fitted object.
#' @export
y_standard_linear <- function(x, beta, sigma = 1) {
  as.numeric(x %*% beta + sigma * stats::rnorm(nrow(x)))
}

#' @describeIn simulate_hdr Cross-validated lasso; `...` is passed to
#'   [glmnet::cv.glmnet()].
#' @export
fit_glmnet_cv <- function(x, y, beta, ...) {
  glmnet::cv.glmnet(x, y, ...)
}

#' @describeIn simulate_hdr True and estimated coefficients at `lambda`
#'   (`"lambda.1se"` or `"lambda.min"`).
#' @param lambda Which `cv.glmnet` choice of `lambda` to use.
#' @export
post_glmnet_cv <- function(fit_obj, x, y, beta, lambda = "lambda.1se") {
  list(true_beta = beta,
       lambda = fit_obj[[lambda]],
       beta_hat = stats::coef(fit_obj, s = lambda))
}

default_x_par <- function(xtype) {
  switch(xtype,
         "toeplitz" = 1 / 3,
         "equi.corr" = 1 / 20,
         "exp.decay" = c(0.4, 5))
}

# Apply fun to seq_len(niters), in parallel if cores > 1, with reproducible
# parallel random number streams when seed is given.
run_iterations <- function(niters, fun, cores = 1L, seed = NULL) {
  if (!is.null(seed)) {
    old_kind <- RNGkind()
    on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]), add = TRUE)
    set.seed(seed, kind = "L'Ecuyer-CMRG")
  }
  if (cores <= 1) return(lapply(seq_len(niters), fun))
  if (.Platform$OS.type != "windows") {
    return(parallel::mclapply(seq_len(niters), fun, mc.cores = cores,
                              mc.set.seed = TRUE))
  }
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterSetRNGStream(cl, if (is.null(seed)) NULL else seed)
  parallel::clusterEvalQ(cl, library(unbiasedgoodness))
  parallel::parLapply(cl, seq_len(niters), fun)
}
