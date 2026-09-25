#' Selective goodness-of-fit test after forward stepwise selection
#'
#' @description
#' `gof_groupfs()` selects a linear model by forward stepwise with an AIC-type
#' stopping rule, enlarges it by continuing forward stepwise, and tests whether
#' the selected model fits as well as the enlarged one. The test statistic is
#' the usual \eqn{F} statistic (or \eqn{\chi} statistic if `sigma` is known),
#' truncated to the selection event so that the p-value is valid conditional
#' on the models that were selected.
#'
#' @details
#' The procedure is:
#'
#' 1. Run forward stepwise (with groups given by `index`) until the penalized
#'    criterion \eqn{n \log(\mathrm{RSS}/n) + k \cdot \mathrm{df}} (or
#'    \eqn{\mathrm{RSS}/\sigma^2 + k \cdot \mathrm{df}} when `sigma` is given)
#'    increases `aicstop` times in a row. The first `s` groups form the
#'    selected (null) model.
#' 2. Continue forward stepwise up to `M = max(ceiling(m * s + b), s + aicstop)`
#'    groups, capped so that the enlarged model is not saturated. Groups
#'    `s + 1, ..., M` form the alternative.
#' 3. Test \eqn{H_0: \mu \in \mathrm{span}(1, X_{\mathrm{null}})} against
#'    \eqn{\mu \in \mathrm{span}(1, X_{\mathrm{null}}, X_{\mathrm{alt}})}.
#'    The statistic is truncated to the event that forward stepwise took the
#'    observed path through step `M` and that the successive penalized
#'    criterion values had the observed signs through step `s + aicstop`.
#'    Conditioning on this (finer) event is sufficient for selective validity.
#'
#' Testing against *all* unselected variables is not possible when
#' \eqn{p > n}, because the full model is saturated; this is why the
#' alternative is a data-driven enlargement of the selected model.
#'
#' @param x Matrix of predictors (n by p).
#' @param y Numeric response vector of length n.
#' @param index Group membership of each column of `x`. The default puts each
#'   column in its own group.
#' @param m Multiplier of the selected model size used to choose the size of
#'   the alternative model.
#' @param b Additive increase in size of the alternative model.
#' @param k Penalty per degree of freedom in the stopping criterion. The
#'   default `log(ncol(x))` is between AIC (`k = 2`) and RIC (`k = 2 * log(p)`).
#' @param sigma Known error standard deviation. If `NULL` (the default) a
#'   truncated \eqn{F} test is used, otherwise a truncated \eqn{\chi} test.
#' @param aicstop Number of consecutive increases of the criterion before
#'   forward stepwise stops.
#' @param maxsteps Maximum number of forward stepwise steps when selecting the
#'   null model. Defaults to `min(G, n - intercept - 2)` for `G` groups, the
#'   largest number of single-column steps that leaves at least one residual
#'   degree of freedom. With `sigma = NULL` and \eqn{p > n} the criterion can
#'   keep decreasing until then, in which case there is no model left to test
#'   against and an error is returned.
#' @param intercept Include an (unpenalized) intercept?
#' @param center,normalize Center and scale the columns of `x`, by group?
#' @return An object of class `"gof_groupfs"`, a list with components
#'   \describe{
#'     \item{pv}{Selective p-value.}
#'     \item{pv_naive}{Unadjusted p-value from the same statistic, which
#'       ignores selection.}
#'     \item{stat, type}{Observed statistic and its type, `"TF"` or `"TC"`.}
#'     \item{df1, df2}{Degrees of freedom (`df2` is `NA` for `"TC"`).}
#'     \item{sigma}{Error standard deviation used, if any.}
#'     \item{null, alt}{Group labels in the selected model and the groups
#'       added in the alternative.}
#'     \item{support}{Truncation region of the statistic, a union of
#'       intervals.}
#'   }
#' @export
#' @examples
#' set.seed(1)
#' n <- 100
#' p <- 60
#' s0 <- 3
#' sim_data <- hdi::rXb(n, p, s0, xtype = "toeplitz", x.par = 1 / 3)
#' x <- sim_data$x
#' y <- as.numeric(x %*% sim_data$beta + rnorm(n))
#' gof_groupfs(x, y, m = 1.2, b = 5)
gof_groupfs <- function(x, y, index = seq_len(ncol(x)), m = 1.2, b = 5,
                        k = log(ncol(x)), sigma = NULL, aicstop = 1,
                        maxsteps, intercept = TRUE, center = TRUE,
                        normalize = TRUE) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- nrow(x)
  G <- length(unique(index))
  if (m < 1 || b <= 0) stop("Need m >= 1 and b > 0 so the alternative is larger")
  if (aicstop < 1) stop("aicstop must be at least 1")
  if (missing(maxsteps)) maxsteps <- min(G, n - intercept - 2)

  fit_sel <- groupfs(x, y, index = index, maxsteps = maxsteps, sigma = sigma,
                     k = k, intercept = intercept, center = center,
                     normalize = normalize, aicstop = aicstop)
  if (!attr(fit_sel, "stopped")) {
    stop("The stopping rule did not stop within maxsteps = ", maxsteps,
         " steps, so there is no larger model to test against.")
  }
  s <- length(fit_sel$action) - aicstop
  null_groups <- fit_sel$action[seq_len(s)]

  # Size of the alternative, capped to keep the full model unsaturated
  M_target <- max(ceiling(m * s + b), s + aicstop)
  M_cap <- min(G, n - intercept - 2)
  M <- max(min(M_target, M_cap), s + aicstop)
  if (M < M_target) {
    warning("Alternative model size reduced from ", M_target, " to ", M,
            " groups to avoid a saturated model.")
  }

  fit_alt <- groupfs(x, y, index = index, maxsteps = M, sigma = sigma, k = k,
                     intercept = intercept, center = center,
                     normalize = normalize, aicstop = 0)
  path_len <- length(fit_sel$action)
  if (length(fit_alt$action) < path_len) {
    stop("Could not continue forward stepwise past the selected model ",
         "without saturating it.")
  }
  if (!identical(fit_alt$action[seq_len(path_len)], fit_sel$action)) {
    stop("Forward stepwise paths disagree; this should not happen.")
  }
  M <- length(fit_alt$action)
  if (M <= s) stop("No room to enlarge the selected model.")
  alt_groups <- fit_alt$action[(s + 1):M]

  out <- gof_groupfs_core(fit_alt, null_groups, alt_groups,
                          aic_steps = path_len)
  out$call <- match.call()
  out
}

# Truncated F / chi test of null_groups against null_groups + alt_groups,
# conditional on the forward stepwise path stored in obj (all of its steps)
# and the signs of successive criterion differences for steps 1..aic_steps.
# condition = FALSE skips the truncation (used to check the statistic against
# anova()).
gof_groupfs_core <- function(obj, null_groups, alt_groups, aic_steps,
                             condition = TRUE) {

  n <- nrow(obj$x)
  x <- obj$x
  y <- obj$y
  index <- obj$index
  sigma <- obj$sigma

  if (length(intersect(null_groups, alt_groups)) > 0) {
    stop("null_groups and alt_groups must be disjoint")
  }
  x_null <- x[, index %in% null_groups, drop = FALSE]
  x_full <- x[, index %in% c(null_groups, alt_groups), drop = FALSE]
  U_null <- if (ncol(x_null) > 0) svdu_thresh(x_null) else matrix(0, n, 0)
  U_full <- svdu_thresh(x_full)
  df_null <- ncol(U_null)
  df_full <- ncol(U_full)
  df1 <- df_full - df_null
  if (df1 < 1) stop("Alternative does not increase the column space")

  # Projection of y onto the null model (zero if the null model is empty;
  # y is already centered when obj$intercept is TRUE)
  Z <- U_null %*% crossprod(U_null, y)

  if (!is.null(sigma)) {
    type <- "TC"
    df2 <- NA_real_
    # Component of y in span(full) orthogonal to span(null)
    U_alt <- svdu_thresh(x_full - U_null %*% crossprod(U_null, x_full))
    R <- crossprod(U_alt, y)
    stat <- sqrt(sum(R^2))
    eta <- U_alt %*% R / stat
    Z <- y - eta * stat
    ydecomp <- list(Z = Z, eta = eta)
    pv_naive <- pchisq((stat / sigma)^2, df1, lower.tail = FALSE)
  } else {
    type <- "TF"
    df2 <- n - obj$intercept - df_full
    if (df2 < 1) stop("Full model is saturated (df2 < 1)")
    C <- df1 / df2
    R1 <- y - Z
    R2 <- y - U_full %*% crossprod(U_full, y)
    R1sq <- sum(R1^2)
    R2sq <- sum(R2^2)
    R <- sqrt(R1sq)
    delta <- R1 - R2
    Vdelta <- delta / sqrt(sum(delta^2))
    V2 <- R2 / sqrt(R2sq)
    stat <- (R1sq - R2sq) / (C * R2sq)
    ydecomp <- list(R = R, Z = Z, Vd = Vdelta, V2 = V2, C = C)
    pv_naive <- pf(stat, df1, df2, lower.tail = FALSE)
  }

  if (condition) {
    intervallist <- c(truncationRegion(obj, ydecomp, type),
                      aic_intervals(obj, ydecomp, type, aic_steps))
    region <- do.call(interval_union, intervallist)
    region <- interval_union(region, Intervals(c(-Inf, 0)))
    E <- interval_complement(region, check_valid = FALSE)
  } else {
    E <- Intervals(c(0, Inf))
  }
  if (length(E) == 0) stop("Empty truncation region")
  if (!any(E[, 1] <= stat & stat <= E[, 2])) {
    warning("Observed statistic lies outside the computed truncation region; ",
            "the p-value may be numerically unreliable.")
  }

  pv <- if (type == "TC") TC_surv(stat, sigma, df1, E) else TF_surv(stat, df1, df2, E)
  if (is.nan(pv)) {
    pv <- 0
    warning("P-value of the form 0/0 converted to 0. This typically occurs ",
            "for numerical reasons when the signal-to-noise ratio is large.")
  }

  out <- list(pv = pv, pv_naive = pv_naive, stat = stat, type = type,
              df1 = df1, df2 = df2, sigma = sigma,
              null = null_groups, alt = alt_groups, support = E)
  class(out) <- "gof_groupfs"
  out
}

# Intervals of the statistic excluded by the signs of successive differences
# of the stopping criterion along the path, for steps 1..aic_steps.
aic_intervals <- function(obj, ydecomp, type, aic_steps) {
  n <- nrow(obj$x)
  k <- attr(obj, "k")
  Z <- ydecomp$Z
  if (type == "TC") {
    eta <- ydecomp$eta
    aic_begin <- sum(obj$y^2) / obj$sigma^2 - n + k * obj$intercept
    pen0 <- k * obj$intercept
  } else {
    Vdelta <- ydecomp$Vd
    V2 <- ydecomp$V2
    R <- ydecomp$R
    C <- ydecomp$C
    aic_begin <- n * (log(2 * pi) + log(mean(obj$y^2))) + k * (1 + n + obj$intercept)
    pen0 <- exp(k * (1 + obj$intercept) / n)
  }
  AICs <- c(aic_begin, obj$AIC)

  ulist <- c(list(matrix(0, n, 1)), obj$maxprojs)
  penlist <- c(pen0, obj$maxpens)
  zlist <- vector("list", aic_steps + 1)
  zlist[[1]] <- zlist[[2]] <- Z
  if (type == "TC") {
    etalist <- vector("list", aic_steps + 1)
    etalist[[1]] <- etalist[[2]] <- eta
  } else {
    vdlist <- v2list <- vector("list", aic_steps + 1)
    vdlist[[1]] <- vdlist[[2]] <- Vdelta
    v2list[[1]] <- v2list[[2]] <- V2
  }
  if (aic_steps > 1) {
    for (step in 1:(aic_steps - 1)) {
      cproj <- obj$cumprojs[[step]]
      zlist[[step + 2]] <- cproj %*% Z
      if (type == "TC") {
        etalist[[step + 2]] <- cproj %*% eta
      } else {
        vdlist[[step + 2]] <- cproj %*% Vdelta
        v2list[[step + 2]] <- cproj %*% V2
      }
    }
  }

  lapply(seq_len(aic_steps), function(step) {
    # Compare the criterion at step with step - 1; the roots functions
    # assume g indexes the smaller value
    peng <- penlist[[step + 1]]
    Ug <- ulist[[step + 1]]
    Uh <- ulist[[step]]
    Zg <- zlist[[step + 1]]
    Zh <- zlist[[step]]
    if (type == "TC") {
      coeffs <- quadratic_coefficients(obj$sigma, Ug, Uh, peng, 0,
                                       etalist[[step + 1]], etalist[[step]],
                                       Zg, Zh)
      if (AICs[step] < AICs[step + 1]) coeffs <- lapply(coeffs, function(cf) -cf)
      quadratic_roots(coeffs$A, coeffs$B, coeffs$C, tol = 1e-15)
    } else {
      coeffs <- TF_coefficients(R, Ug, Uh, peng, 1, Zg, Zh,
                                vdlist[[step + 1]], vdlist[[step]],
                                v2list[[step + 1]], v2list[[step]])
      if (AICs[step] < AICs[step + 1]) coeffs <- lapply(coeffs, function(cf) -cf)
      TF_roots(R, C, coeffs)
    }
  })
}

#' @export
print.gof_groupfs <- function(x, ...) {
  cat("Selective goodness-of-fit test after forward stepwise\n\n")
  cat("Selected groups (", length(x$null), "): ",
      paste(x$null, collapse = ", "), "\n", sep = "")
  cat("Added in alternative (", length(x$alt), "): ",
      paste(x$alt, collapse = ", "), "\n\n", sep = "")
  tab <- data.frame(stat = x$stat, df1 = x$df1, df2 = x$df2,
                    trunc_min = min(x$support), trunc_max = max(x$support),
                    p_selective = x$pv, p_naive = x$pv_naive)
  names(tab)[1] <- x$type
  print(tab, row.names = FALSE, digits = 4)
  invisible(x)
}
