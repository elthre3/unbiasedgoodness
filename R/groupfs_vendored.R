# Forward stepwise selection with groups and the truncation-region
# machinery used for selective F / chi tests.
#
# Adapted from the selectiveInference package, version 1.2.5
# (R/funs.groupfs.R, R/funs.quadratic.R, R/funs.common.R), by
# Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Joshua Loftus,
# Stephen Reid and Jelena Markovic. The groupfs code was written by
# Joshua Loftus. Only the functions needed by gof_groupfs() are kept.
#
# Changes from upstream:
# - scaleGroups() centers each column rather than subtracting a single
#   group-wide mean, so that centered columns are orthogonal to the
#   intercept for groups with more than one column.
# - TF_roots() loops over seq_along(negative) instead of 1:length(negative),
#   and returns (-Inf, 0] when no region is excluded; upstream produced
#   spurious NA intervals in that case.
# - add1.groupfs() renamed add1_groupfs() so it is not taken for an S3
#   method of stats::add1().
# - Roxygen blocks removed; all functions here are internal.

groupfs <- function(x, y, index, maxsteps, sigma = NULL, k = 2, intercept = TRUE, center = TRUE, normalize = TRUE, aicstop = 0, verbose = FALSE) {

  if (missing(index)) stop("Missing argument: index.")
  p <- ncol(x)
  n <- nrow(x)

  # Group labels
  labels <- unique(index)
  G <- length(labels)
  inactive <- labels
  active <- c()

  if (missing(maxsteps) || maxsteps >= min(n, G)) maxsteps <- min(n-1, G)
  checkargs.xy(x=x, y=y)
  checkargs.groupfs(x, index, maxsteps)
  if (maxsteps > G) stop("maxsteps is larger than number of groups")
  gsizes <- sort(rle(sort(index))$lengths, decreasing = TRUE)
  if (sum(gsizes[1:maxsteps]) >= nrow(x)) {
      maxsteps <- max(which(cumsum(gsizes) < nrow(x)))
      warning(paste("If the largest groups are included the model will be saturated/overdetermined. To prevent this maxsteps has been changed to", maxsteps))
  }

  # Initialize copies of data for loop
  by <- mean(y)
  y.update <- y
  if (intercept) y.update <- y - by
  y.last <- y.update

  # Center and scale design matrix
  xscaled <- scaleGroups(x, index, center, normalize)
  xm <- xscaled$xm
  xs <- xscaled$xs
  x.update <- xscaled$x

  x.begin <- x.update
  y.begin <- y.update
  stopped <- FALSE
  # Store all projections computed along the path
  terms = projections = maxprojs = aicpens = maxpens = cumprojs = vector("list", maxsteps)

  # Store other information from each step
  path.info <- data.frame(imax=integer(maxsteps), df=integer(maxsteps), AIC=numeric(maxsteps), RSS=numeric(maxsteps), RSSdrop=numeric(maxsteps), chisq=numeric(maxsteps))

  modelrank <- as.numeric(intercept)
  if (is.null(sigma)) {
      modelrank <- modelrank + 1
      aic.begin <- aic.last <- n*(log(2*pi) + log(mean(y.update^2))) + k * (n + modelrank)
  } else {
      aic.begin <- aic.last <- sum(y.update^2)/sigma^2 - n + k * modelrank
  }
  if (verbose) print(paste0("Start:  AIC=", round(aic.begin, 3)), quote = FALSE)

  # Begin main loop
  for (step in 1:maxsteps) {

    added <- add1_groupfs(x.update, y.update, index, labels, inactive, k, sigma)

    # Group to be added
    imax <- added$imax
    inactive <- setdiff(inactive, imax)
    active <- union(active, imax)
    inactive.inds <- which(!index %in% active)

    # Rank of group
    modelrank <- modelrank + added$df

    # Stop without adding if model has become saturated
    if (modelrank >= n) {
        stop("Saturated model. Abandon ship!")
    }

    # Regress added group out of y and inactive x
    P.imax <- added$maxproj %*% t(added$maxproj)
    P.imax <- diag(rep(1, n)) - P.imax
    y.update <- P.imax %*% y.update
    x.update[, inactive.inds] <- P.imax %*% x.update[, inactive.inds]

    # Compute AIC
    if (is.null(sigma)) {
        added$AIC <- n * log(added$maxterm/n) - k * added$df + n*log(2*pi) + k * (n + modelrank)
    } else {
        added$AIC <- sum(y.update^2)/sigma^2 - n + k * modelrank
    }

    projections[[step]] <- added$projections
    maxprojs[[step]] <- added$maxproj
    aicpens[[step]] <- added$aicpens
    maxpens[[step]] <- added$maxpen
    if (step == 1) cumprojs[[step]] <- P.imax
    if (step > 1) cumprojs[[step]] <- P.imax %*% cumprojs[[step-1]]
    terms[[step]] <- added$terms

    # Compute RSS for unadjusted chisq p-values
    added$RSS <- sum(y.update^2)
    scale.chisq <- 1

    added$RSSdrop <- sum((y.last - y.update)^2)
    added$chisq <- pchisq(added$RSSdrop/scale.chisq, lower.tail=FALSE, df = added$df)
    y.last <- y.update

    # Projections are stored separately
    step.info <- data.frame(added[-c(3:(length(added)-4))])
    path.info[step, ] <- step.info

    if (verbose) print(round(step.info, 3))

    if (aicstop > 0 && step < maxsteps && step >= aicstop && aic.last < added$AIC) {
        if (all(diff(c(aic.begin, path.info$AIC)[(step+1-aicstop):(step+1)]) > 0)) {

            if (is.null(sigma)) {
                added$AIC <- n * log(added$maxterm/n) - k * added$df + n + n*log(2*pi) + k * modelrank
            } else {
                added$AIC <- sum(y.update^2)/sigma^2 - n + k * modelrank
            }

            path.info <- path.info[1:step, ]
            projections[(step+1):maxsteps] <- NULL
            maxprojs[(step+1):maxsteps] <- NULL
            aicpens[(step+1):maxsteps] <- NULL
            maxpens[(step+1):maxsteps] <- NULL
            cumprojs[(step+1):maxsteps] <- NULL
            terms[(step+1):maxsteps] <- NULL
            maxsteps <- step
            stopped <- TRUE
            break
        }
    }
    aic.last <- added$AIC
  }

  # Is there a better way of doing this?
  # Use some projections already computed?
  beta <- coef(lm(y.begin ~ x.begin[,index %in% path.info$imax]-1))
  names(beta) <- index[index %in% path.info$imax]

  # Create output object
  value <- list(action = path.info$imax, L = path.info$L, AIC = path.info$AIC, projections = projections, maxprojs = maxprojs, aicpens = aicpens, maxpens = maxpens, cumprojs = cumprojs, log = path.info, index = index, y = y.begin, x = x.begin, coefficients = beta, bx = xm, by = by, sx = xs, sigma = sigma, intercept = intercept, call = match.call(), terms = terms)

  class(value) <- "groupfs"
  attr(value, "center") <- center
  attr(value, "normalize") <- normalize
  attr(value, "labels") <- labels
  attr(value, "maxsteps") <- maxsteps
  attr(value, "sigma") <- sigma
  attr(value, "k") <- k
  attr(value, "aicstop") <- aicstop
  attr(value, "stopped") <- stopped
  if (is.null(attr(x, "varnames"))) {
    attr(value, "varnames") <- colnames(x)
  } else {
    attr(value, "varnames") <- attr(x, "varnames")
  }
  return(value)
}

add1_groupfs <- function(xr, yr, index, labels, inactive, k, sigma = NULL) {

  # Use characters to avoid issues where
  # list() populates NULL lists in the positions
  # of the active variables
  ### Question for later: does this slow down lapply?
  keys = as.character(inactive)
  n <- nrow(xr)

  # Compute sums of squares to determine which group is added
  # penalized by rank of group if k > 0
  projections = aicpens = terms = vector("list", length(keys))
  names(projections) = names(terms) = names(aicpens) = keys
  for (key in keys) {
      inds <- which(index == key)
      xi <- xr[,inds]
      ui <- svdu_thresh(xi)
      dfi <- ncol(ui)
      projections[[key]] <- ui
      uy <- t(ui) %*% yr
      if (is.null(sigma)) {
          aicpens[[key]] <- exp(k*dfi/n)
          terms[[key]] <- (sum(yr^2) - sum(uy^2)) * aicpens[[key]]
      } else {
          aicpens[[key]] <- sigma^2 * k * dfi
          terms[[key]] <- (sum(yr^2) - sum(uy^2)) + aicpens[[key]]
      }
  }

  # Maximizer = group to be added
  terms.optind <- which.min(terms)
  imax <- inactive[terms.optind]
  optkey <- which(keys == imax)
  maxproj <- projections[[optkey]]
  maxpen <- aicpens[[optkey]]
  maxterm <- terms[[optkey]]
  projections[[optkey]] <- NULL
  aicpens[[optkey]] <- NULL

  return(list(imax=imax, df = ncol(maxproj), projections = projections, maxproj = maxproj, aicpens = aicpens, maxpen = maxpen, maxterm = maxterm, terms = terms))
}


TC_surv <- function(TC, sigma, df, E) {
    if (length(E) == 0) {
        stop("Empty TC support")
    }

    # Sum truncated cdf over each part of E
    denom <- do.call(sum, lapply(1:nrow(E), function(v) {
      tchi_interval(E[v,1], E[v,2], sigma, df)
    }))

    # Sum truncated cdf from observed value to max of
    # truncation region
    numer <- do.call(sum, lapply(1:nrow(E), function(v) {
      lower <- E[v,1]
      upper <- E[v,2]
      if (upper > TC) {
        # Observed value is left of this interval's right endpoint
        if (lower < TC) {
          # Observed value is in this interval
          return(tchi_interval(TC, upper, sigma, df))
        } else {
          # Observed value is not in this interval
          return(tchi_interval(lower, upper, sigma, df))
        }
      } else {
        # Observed value is right of this entire interval
        return(0)
      }
    }))

    # Survival function
    value <- numer/denom
    # Force p-value to lie in the [0,1] interval
    # in case of numerical issues
    value <- max(0, min(1, value))
    value
}

tchi_interval <- function(lower, upper, sigma, df) {
  a <- (lower/sigma)^2
  b <- (upper/sigma)^2
  if (b == Inf) {
      integral <- pchisq(a, df, lower.tail = FALSE)
  } else {
      integral <- pchisq(b, df) - pchisq(a, df)
  }
  if ((integral < .Machine$double.eps) && (b < Inf)) {
      integral <- num_int_chi(a, b, df)
  }
  return(integral)
}

num_int_chi <- function(a, b, df, nsamp = 10000) {
  grid <- seq(from=a, to=b, length.out=nsamp)
  integrand <- dchisq(grid, df)
  return((b-a)*mean(integrand))
}

TF_surv <- function(TF, df1, df2, E) {
    if (length(E) == 0) {
        stop("Empty TF support")
    }

    # Sum truncated cdf over each part of E
    denom <- do.call(sum, lapply(1:nrow(E), function(v) {
      TF_interval(E[v,1], E[v,2], df1, df2)
    }))

    # Sum truncated cdf from observed value to max of
    # truncation region
    numer <- do.call(sum, lapply(1:nrow(E), function(v) {
      lower <- E[v,1]
      upper <- E[v,2]
      if (upper > TF) {
        # Observed value is left of this interval's right endpoint
        if (lower < TF) {
          # Observed value is in this interval
          return(TF_interval(TF, upper, df1, df2))
        } else {
          # Observed value is not in this interval
          return(TF_interval(lower, upper, df1, df2))
        }
      } else {
        # Observed value is right of this entire interval
        return(0)
      }
    }))

    # Survival function
    value <- numer/denom
    # Force p-value to lie in the [0,1] interval
    # in case of numerical issues
    #value <- max(0, min(1, value))
    value
}

TF_interval <- function(lower, upper, df1, df2) {
  a <- lower
  b <- upper
  if (b == Inf) {
      integral <- pf(a, df1, df2, lower.tail = FALSE)
  } else {
      integral <- pf(b, df1, df2) - pf(a, df1, df2)
  }
  if ((integral < .Machine$double.eps) && (b < Inf)) {
      integral <- num_int_F(a, b, df1, df2)
  }
  return(integral)
}

num_int_F <- function(a, b, df1, df2, nsamp = 10000) {
  grid <- seq(from=a, to=b, length.out=nsamp)
  integrand <- df(grid, df1, df2)
  return((b-a)*mean(integrand))
}

scaleGroups <- function(x, index, center = TRUE, normalize = TRUE) {
  keys <- unique(index)
  xm <- rep(0, ncol(x))
  xs <- rep(1, ncol(x))

  for (j in keys) {
    inds <- which(index == j)
    if (center) {
        xmj <- colMeans(x[, inds, drop = FALSE])
        xm[inds] <- xmj
        x[, inds] <- sweep(x[, inds, drop = FALSE], 2, xmj)
    }
    normsq <- sum(x[, inds]^2)
    xsj <- sqrt(normsq)
    xs[inds] <- xsj
    if (xsj > 0) {
        if (normalize) x[, inds] <- x[, inds] / xsj
    } else {
        stop(paste("Design matrix contains identically zero group of variables:", j))
    }
  }
  return(list(x=x, xm=xm, xs=xs))
}

svdu_thresh <- function(x) {
    svdx <- svd(x)
    inds <- svdx$d > svdx$d[1] * sqrt(.Machine$double.eps)
    return(svdx$u[, inds, drop = FALSE])
}

checkargs.groupfs <- function(x, index, maxsteps) {
    if (length(index) != ncol(x)) stop("Length of index does not match number of columns of x")
    if ((round(maxsteps) != maxsteps) || (maxsteps <= 0)) stop("maxsteps must be an integer > 0")
}

checkargs.xy <- function(x, y) {
  if (missing(x)) stop("x is missing")
  if (is.null(x) || !is.matrix(x)) stop("x must be a matrix")
  if (missing(y)) stop("y is missing")
  if (is.null(y) || !is.numeric(y)) stop("y must be numeric")
  if (ncol(x) == 0) stop("There must be at least one predictor [must have ncol(x) > 0]")
  if (checkcols(x)) stop("x cannot have duplicate columns")
  if (length(y) == 0) stop("There must be at least one data point [must have length(y) > 0]")
  if (length(y)!=nrow(x)) stop("Dimensions don't match [length(y) != nrow(x)]")
}

checkcols <- function(A) {
  b = rnorm(nrow(A))
  a = sort(t(A)%*%b)
  if (any(diff(a)==0)) return(TRUE)
  return(FALSE)
}

truncationRegion <- function(obj, ydecomp, type, tol = 1e-15) {

  n <- nrow(obj$x)
  Z <- ydecomp$Z
  if (type == "TC") {
      eta <- ydecomp$eta
  } else {
      Vd <- ydecomp$Vd
      V2 <- ydecomp$V2
      C <- ydecomp$C
      R <- ydecomp$R
  }
  L <- lapply(1:length(obj$action), function(s) {

    Ug <- obj$maxprojs[[s]]
    peng <- obj$maxpens[[s]]
    if (s > 1) {
        Zs <- obj$cumprojs[[s-1]] %*% Z
        if (type == "TC") {
            etas <- obj$cumprojs[[s-1]] %*% eta
        } else {
            Vds <- obj$cumprojs[[s-1]] %*% Vd
            V2s <- obj$cumprojs[[s-1]] %*% V2
        }
    } else {
        Zs <- Z
        if (type == "TC") {
            etas <- eta
        } else {
            Vds <- Vd
            V2s <- V2
        }
    }

    num.projs <- length(obj$projections[[s]])
    if (num.projs == 0) {
        return(list(Intervals(c(-Inf,0))))
    } else {
      lapply(1:num.projs, function(l) {

          Uh <- obj$projections[[s]][[l]]
          penh <- obj$aicpens[[s]][[l]]
          # The quadratic form corresponding to
          # (t*U + Z)^T %*% Q %*% (t*U + Z) \geq 0
          # we find the roots in t, if there are any
          # and return the interval of potential t
          if (type == "TC") {
              coeffs <- quadratic_coefficients(obj$sigma, Ug, Uh, peng, penh, etas, etas, Zs, Zs)
              quadratic_roots(coeffs$A, coeffs$B, coeffs$C, tol)
          } else {
              coeffs <- TF_coefficients(R, Ug, Uh, peng, penh, Zs, Zs, Vds, Vds, V2s, V2s)
              roots <- TF_roots(R, C, coeffs)
              return(roots)
          }
      })
    }
    # LL is a list of intervals
  })
  # L is now a list of lists of intervals
  return(unlist(L, recursive = FALSE, use.names = FALSE))
}

quadratic_coefficients <- function(sigma, Ug, Uh, peng, penh, etag, etah, Zg, Zh) {
    # g indexes minimizer, h the comparison
    Uheta <- t(Uh) %*% etah
    Ugeta <- t(Ug) %*% etag
    UhZ <- t(Uh) %*% Zh
    UgZ <- t(Ug) %*% Zg
    etaZh <- t(etah) %*% Zh
    etaZg <- t(etag) %*% Zg
    if (is.null(sigma)) {
        A <- penh * (sum(etah^2) - sum(Uheta^2)) - peng * (sum(etag^2) - sum(Ugeta^2))
        B <- 2 * penh * (etaZh - t(Uheta) %*% UhZ) - 2 * peng * (etaZg - t(Ugeta) %*% UgZ)
        C <- penh * (sum(Zh^2) - sum(UhZ^2)) - peng * (sum(Zg^2) - sum(UgZ^2))
    } else {
        A <- (sum(etah^2) - sum(Uheta^2)) - (sum(etag^2) - sum(Ugeta^2))
        B <- 2 * (etaZh - t(Uheta) %*% UhZ) - 2 * (etaZg - t(Ugeta) %*% UgZ)
        C <- (sum(Zh^2) - sum(UhZ^2) + penh) - (sum(Zg^2) - sum(UgZ^2) + peng)
    }
    return(list(A = A, B = B, C = C))
}

quadratic_roots <- function(A, B, C, tol) {
    disc <- B^2 - 4*A*C
    b2a <- -B/(2*A)

    if (disc > tol) {
        # Real roots
        pm <- sqrt(disc)/(2*A)
        endpoints <- sort(c(b2a - pm, b2a + pm))

    } else {
        # No real roots
        if (A > -tol) {
          # Quadratic form always positive
            return(Intervals(c(-Inf,0)))
        } else {
          # Quadratic form always negative
            stop("Empty TC support is infeasible")
        }
    }

    if (A > tol) {
        # Parabola opens upward
        if (min(endpoints) > 0) {
          # Both roots positive, union of intervals
            return(Intervals(rbind(c(-Inf,0), endpoints)))
        } else {
          # At least one negative root
            return(Intervals(c(-Inf, max(0, endpoints[2]))))
        }
    } else {
        if (A < -tol) {
          # Parabola opens downward
            if (endpoints[2] < 0) {
            # Positive quadratic form only when t negative
                stop("Negative TC support is infeasible")
            } else {
            # Part which is positive
                if (endpoints[1] > 0) {
                    return(Intervals(rbind(c(-Inf, endpoints[1]), c(endpoints[2], Inf))))
                } else {
                    return(Intervals(c(endpoints[2], Inf)))
                }
            }
        } else {
          # a is too close to 0, quadratic is actually linear
            if (abs(B) > tol) {
                if (B > 0) {
                    return(Intervals(c(-Inf, max(0, -C/B))))
                } else {
                    if (-C/B < 0) stop("Infeasible linear equation")
                    return(Intervals(rbind(c(-Inf, 0), c(-C/B, Inf))))
                }
            } else {
                warning("Ill-conditioned quadratic")
                return(Intervals(c(-Inf,0)))
            }
        }
    }
}

# Helper functions for TF roots
roots_to_checkpoints <- function(roots) {
    checkpoints <- unique(sort(c(0, roots)))
    return(c(0, (checkpoints + c(checkpoints[-1], 200 + checkpoints[length(checkpoints)]))/2))
}
roots_to_partition <- function(roots) {
    checkpoints <- unique(sort(c(0, roots)))
    return(list(endpoints = c(checkpoints, Inf), midpoints = (checkpoints + c(checkpoints[-1], 200 + checkpoints[length(checkpoints)]))/2))
}

# Efficiently compute coefficients of one-dimensional TF slice function
TF_coefficients <- function(R, Ug, Uh, peng, penh, Zg, Zh, Vdg, Vdh, V2g, V2h) {

    UhZ <- t(Uh) %*% Zh
    UgZ <- t(Ug) %*% Zg
    UhVd <- t(Uh) %*% Vdh
    UgVd <- t(Ug) %*% Vdg
    UhV2 <- t(Uh) %*% V2h
    UgV2 <- t(Ug) %*% V2g
    VdZh <- sum(Vdh*Zh)
    VdZg <- sum(Vdg*Zg)
    V2Zh <- sum(V2h*Zh)
    V2Zg <- sum(V2g*Zg)

    x0 <- penh * (sum(Zh^2) - sum(UhZ^2)) - peng * (sum(Zg^2) - sum(UgZ^2))
    x1 <- 2*R*(penh * (VdZh - sum(UhZ*UhVd)) - peng * (VdZg - sum(UgZ*UgVd)))
    x2 <- 2*R*(penh * (V2Zh - sum(UhZ*UhV2)) - peng * (V2Zg - sum(UgZ*UgV2)))
    x12 <- 2*R^2*(penh * (sum(Vdh*V2h) - sum(UhVd*UhV2)) - peng * (sum(Vdg*V2g) - sum(UgVd*UgV2)))
    x11 <- R^2*(penh * (sum(Vdh^2) - sum(UhVd^2)) - peng * (sum(Vdg^2) - sum(UgVd^2)))
    x22 <- R^2*(penh * (sum(V2h^2) - sum(UhV2^2)) - peng * (sum(V2g^2) - sum(UgV2^2)))

    return(list(x11=x11, x22=x22, x12=x12, x1=x1, x2=x2, x0=x0))
}

# Numerically solve for roots of TF slice using
# hybrid polyroot/uniroot approach

TF_roots <- function(R, C, coeffs, tol = 1e-8, tol2 = 1e-6) {

    x11 <- coeffs$x11
    x22 <- coeffs$x22
    x12 <- coeffs$x12
    x1 <- coeffs$x1
    x2 <- coeffs$x2
    x0 <- coeffs$x0

    g1 <- function(t) sqrt(C*t/(1+C*t))
    g2 <- function(t) 1/sqrt(1+C*t)
    I <- function(t) x11*g1(t)^2 + x12*g1(t)*g2(t) + x22*g2(t)^2 + x1*g1(t) + x2*g2(t) + x0

    z4 <- complex(real = -x11 + x22, imaginary = -x12)/4
    z3 <- complex(real = x2, imaginary = -x1)/2
    z2 <- complex(real = x11/2+x22/2+x0)
    z1 <- Conj(z3)
    z0 <- Conj(z4)

    zcoefs <- c(z0, z1, z2, z3, z4)
    croots <- polyroot(zcoefs)
    thetas <- Arg(croots)
    # Can't specify polyroot precision :(
    modinds <- Mod(croots) <= 1 + tol2 & Mod(croots) >= 1 - tol2
    angleinds <- thetas >=0 & thetas <= pi/2
    roots <- unique(thetas[which(modinds & angleinds)])
    troots <- tan(roots)^2/C

    checkpoints <- c()
    if (length(troots) > 0) checkpoints <- roots_to_checkpoints(troots)
    checkpoints <- sort(
        c(checkpoints, 0, tol, tol2,
                        seq(from = sqrt(tol2), to = 1, length.out = 50),
                        seq(from = 1.2, to=50, length.out = 20),
                        100, 1000, 10000))
    ## if (length(troots) == 0) {
    ##     # Polyroot didn't catch any roots
    ##     # ad-hoc check:
    ##     checkpoints <- c(0, tol, tol2,
    ##                      seq(from = sqrt(tol2), to = 1, length.out = 50),
    ##                      seq(from = 1.2, to=50, length.out = 20),
    ##                      100, 1000, 10000)
    ## } else {
    ##     checkpoints <- roots_to_checkpoints(troots)
    ## }

    signs <- sign(I(checkpoints))
    diffs <- c(0, diff(signs))
    changeinds <- which(diffs != 0)

    if (length(changeinds) > 0) {

        roots <- unlist(lapply(changeinds, function(ind) {
            uniroot(I, lower = checkpoints[ind-1], upper = checkpoints[ind], tol = tol)$root
        }))

        partition <- roots_to_partition(roots)
        negative <- which(I(partition$midpoints) < 0)
        if (length(negative) == 0) return(Intervals(c(-Inf,0)))

        intervals <- matrix(NA, ncol=2)
        for (i in seq_along(negative)) {
            ind <- negative[i]
            if ((i > 1) && (ind == negative[i-1] + 1)) {
                # There was not a sign change at end of previous interval
                intervals[nrow(intervals), 2] <- partition$endpoints[ind+1]
            } else {
                intervals <- rbind(intervals, c(partition$endpoints[ind], partition$endpoints[ind+1]))
            }
        }

        return(Intervals(intervals[-1,]))
    }

    # Apparently no roots, always positive
    if (I(0) < 0) stop("Infeasible constraint!")
    return(Intervals(c(-Inf,0)))
}

