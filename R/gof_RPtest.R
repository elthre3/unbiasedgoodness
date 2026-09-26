#' Goodness-of-fit tests after lasso selection using residual prediction
#'
#' @description
#' `gof_RPtest()` selects a model with the cross-validated lasso, enlarges it
#' by moving further along the lasso path, and uses [RPtests::RPtest()] to
#' test whether the variables in the larger model can predict the residuals
#' of the selected model.
#'
#' @details
#' The selected (null) model is the lasso support at `cv_fit[[lambda]]`, of
#' size `s`. The alternative adds the variables that enter the lasso path by
#' the first `lambda` at which the support has at least `ceiling(m * s + b)`
#' variables, or at the end of the path if that size is never reached.
#'
#' `RPtest` with `resid_type = "Lasso"` fits a lasso to the null variables, so
#' the null model needs at least two variables. If fewer are selected, the
#' null model is taken from the first `lambda` on the path with at least two.
#'
#' The p-values are not adjusted for selection. Whether they control error
#' rates conditional on the selected model is the question studied in the
#' package articles.
#'
#' @param x A matrix of predictors (n by p).
#' @param y A numeric vector of outcomes.
#' @param m Multiplier of the selected model size.
#' @param b Additive increase in the size of the alternative model.
#' @param lambda Which `cv.glmnet` choice of `lambda` defines the selected
#'   model.
#' @param B Number of bootstrap samples used by each `RPtest` call.
#' @param ... Further arguments passed to [RPtests::RPtest()].
#' @return A list with components
#'   \describe{
#'     \item{selected_support}{Columns of `x` in the selected model.}
#'     \item{alt_support}{Columns added in the alternative model.}
#'     \item{pval_OLS, pval_lasso}{Tests of the added columns using OLS or
#'       square-root lasso residuals of the selected model.}
#'     \item{pval_lasso_beta}{As `pval_lasso`, but the residual bootstrap is
#'       generated from the lasso coefficients of the selected model instead
#'       of refitting `cv.glmnet` inside `RPtest`.}
#'     \item{pval_full, pval_full_lasso}{Tests using all unselected columns
#'       as the alternative, with OLS or square-root lasso residuals.}
#'   }
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 100
#' p <- 200
#' s0 <- 5
#' sim_data <- rXb(n, p, s0, xtype = "equi.corr", x.par = 1 / 20)
#' x <- sim_data$x
#' y <- as.numeric(x %*% sim_data$beta + rnorm(n))
#' gof_RPtest(x, y, B = 19L)
#' }
gof_RPtest <- function(x, y, m = 1.2, b = 5,
                       lambda = c("lambda.1se", "lambda.min"), B = 49L, ...) {
  lambda <- match.arg(lambda)
  x <- as.matrix(x)
  y <- as.numeric(y)
  if (m < 1 || b <= 0) stop("Need m >= 1 and b > 0 so the alternative is larger")

  cv_fit <- glmnet::cv.glmnet(x, y)
  path <- cv_fit$glmnet.fit

  lambda_null <- cv_fit[[lambda]]
  beta_null <- stats::coef(path, s = lambda_null)
  if (sum(beta_null[-1] != 0) < 2) {
    warning("Fewer than 2 variables selected by cv.glmnet; using the first ",
            "lambda with at least 2 for compatibility with RPtest")
    lambda_null <- path$lambda[which(path$df >= 2)[1]]
    beta_null <- stats::coef(path, s = lambda_null)
  }
  support_null <- which(beta_null[-1] != 0)
  s <- length(support_null)

  target <- ceiling(m * s + b)
  idx <- which(path$df >= target)
  if (length(idx) == 0) {
    idx <- length(path$lambda)
    warning("Lasso path never reaches ", target, " variables; using the end ",
            "of the path (", path$df[idx], " variables)")
  }
  lambda_alt <- path$lambda[idx[1]]
  beta_alt <- stats::coef(path, s = lambda_alt)
  support_new <- setdiff(which(beta_alt[-1] != 0), support_null)
  if (length(support_new) == 0) {
    stop("Failed to include additional variables in the larger model")
  }

  x_null <- x[, support_null, drop = FALSE]
  rp <- function(resid_type, x_alt, ...) {
    quiet_sqrt_lasso(RPtests::RPtest(x_null, y, resid_type = resid_type,
                                     test = "group", x_alt = x_alt, B = B, ...))
  }
  x_new <- x[, support_new, drop = FALSE]
  x_rest <- x[, -support_null, drop = FALSE]

  list(selected_support = support_null,
       alt_support = support_new,
       lambda_null = lambda_null,
       lambda_alt = lambda_alt,
       pval_OLS = rp("OLS", x_new, ...),
       pval_lasso = rp("Lasso", x_new, ...),
       pval_lasso_beta = rp("Lasso", x_new,
                            beta_est = as.numeric(beta_null[-1])[support_null],
                            ...),
       pval_full = rp("OLS", x_rest, ...),
       pval_full_lasso = rp("Lasso", x_rest, ...))
}

# RPtest warns "Smallest lambda chosen by sqrt_lasso" whenever the null model
# is low-dimensional, which is always the case here.
quiet_sqrt_lasso <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("sqrt_lasso", conditionMessage(w))) invokeRestart("muffleWarning")
  })
}
