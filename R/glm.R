#' Survey-weighted generalised linear models
#'
#' Same estimator and linearised (sandwich) variance as `survey::svyglm()`.
#' The IRLS cross-products are computed in parallel C++; families other than
#' gaussian, binomial, quasibinomial, poisson, quasipoisson and Gamma fall
#' back to `stats::glm.fit()`. The model is fitted in the current domain
#' (`filter()`); groups are ignored.
#'
#' @param .data An `fsvy` design.
#' @param formula Model formula.
#' @param family A family object, as in `glm()`.
#' @param control `glm.control()` settings.
#' @return An object of class `fs_glm` with `coef()`, `vcov()`, `summary()`,
#'   `confint()` and `tidy()` methods.
#' @examples
#' set.seed(1)
#' d <- data.frame(psu = rep(1:40, each = 5), str = rep(1:8, each = 25),
#'                 w = runif(200, 1, 3), x = rnorm(200))
#' d$y <- rbinom(200, 1, plogis(-0.5 + d$x))
#' des <- as_fsvy(d, ids = psu, strata = str, weights = w)
#' fit <- fs_glm(des, y ~ x, family = quasibinomial())
#' tidy(fit, conf.int = TRUE, exponentiate = TRUE)
#' @export
fs_glm <- function(.data, formula, family = stats::gaussian(), control = stats::glm.control()) {
  design <- .data
  if (is.function(family)) family <- family()
  if (is.character(family)) family <- get(family, mode = "function")()
  dom <- domain_mask(design)
  w_all <- design$w
  # survey rescales the weights to mean 1 within the (sub)design
  pw <- w_all / mean(w_all[dom])

  mf <- stats::model.frame(formula, design$data[dom, , drop = FALSE], na.action = stats::na.omit)
  used <- which(dom)
  nas <- attr(mf, "na.action")
  if (length(nas)) used <- used[-nas]
  mt <- attr(mf, "terms")
  X <- stats::model.matrix(mt, mf)
  y <- stats::model.response(mf, "any")
  if (is.factor(y)) y <- as.double(y != levels(y)[1L])
  offset <- stats::model.offset(mf)
  if (is.null(offset)) offset <- rep(0, nrow(X))
  wts <- pw[used]

  fit <- irls_fit(X, as.double(y), wts, offset, family, control)

  # estimating functions, linearised over the design (svy.varcoef)
  p <- fit$rank
  keep <- fit$keep
  Ainv <- fit$Ainv
  e <- numeric(length(w_all))
  e[used] <- fit$resid * fit$wt
  gid <- rep(-1L, length(w_all))
  gid[used] <- 0L
  zrow <- rep(0L, length(w_all))
  zrow[used] <- seq_along(used) - 1L
  kr <- kernel_num(design, X[, keep, drop = FALSE], gid, 1L, full = TRUE, mult = e, zrow = zrow)
  meat <- matrix(kr$V[1, ], p, p)
  V <- Ainv %*% meat %*% Ainv

  cf <- rep(NA_real_, ncol(X))
  names(cf) <- colnames(X)
  cf[keep] <- fit$coef
  dimnames(V) <- list(colnames(X)[keep], colnames(X)[keep])

  sub <- design
  sub$domain <- seq_along(w_all) %in% used
  df_res <- fs_degf(sub) + 1 - p

  structure(list(coefficients = cf, vcov = V, df.residual = df_res, family = family,
                 formula = formula, iter = fit$iter, converged = fit$converged,
                 n = length(used), deviance = fit$dev, fitted = fit$mu, rank = p),
            class = "fs_glm")
}

fast_families <- c("gaussian", "binomial", "quasibinomial", "poisson", "quasipoisson", "Gamma")

# glm.fit's IRLS, with X'WX done in C++ threads. Returns the pieces svyglm
# uses: final working residuals, the working weights of the last iteration and
# (X'WX)^-1 from that iteration.
irls_fit <- function(X, y, w, offset, family, control) {
  n <- nrow(X)
  if (!family$family %in% fast_families) return(glmfit_fallback(X, y, w, offset, family, control))
  linkinv <- family$linkinv
  variance <- family$variance
  mu.eta <- family$mu.eta
  dev.resids <- family$dev.resids
  nobs <- n
  weights <- w
  etastart <- NULL; mustart <- NULL
  eval(family$initialize)
  eta <- family$linkfun(mustart)
  mu <- linkinv(eta)
  devold <- sum(dev.resids(y, mu, weights))
  keep <- seq_len(ncol(X))
  coef <- NULL
  conv <- FALSE
  nth <- fs_threads()

  for (iter in seq_len(control$maxit)) {
    mev <- mu.eta(eta)
    good <- weights > 0 & mev != 0
    z <- (eta - offset) + (y - mu) / mev
    W2 <- weights * mev^2 / variance(mu)
    W2[!good] <- 0
    z[!good] <- 0
    cp <- fs_xtwx(X, W2, z, nth)
    A <- cp$XtWX
    b <- cp$XtWz
    if (iter == 1L) {
      keep <- alias_check(X, W2, A)
      if (length(keep) < ncol(X)) {
        X <- X[, keep, drop = FALSE]
        A <- A[keep, keep, drop = FALSE]
        b <- b[keep]
      }
    }
    R <- chol(A)
    start <- backsolve(R, forwardsolve(t(R), b))
    eta_new <- drop(X %*% start) + offset
    mu_new <- linkinv(eta_new)
    dev <- sum(dev.resids(y, mu_new, weights))
    # step halving on invalid fits, as glm.fit
    if (!is.finite(dev) || !(family$valideta(eta_new) && family$validmu(mu_new))) {
      if (is.null(coef)) stop("No valid set of coefficients: please supply starting values", call. = FALSE)
      ii <- 1
      while (!is.finite(dev) || !(family$valideta(eta_new) && family$validmu(mu_new))) {
        if (ii > control$maxit) stop("inner loop; cannot correct step size", call. = FALSE)
        ii <- ii + 1
        start <- (start + coef) / 2
        eta_new <- drop(X %*% start) + offset
        mu_new <- linkinv(eta_new)
        dev <- sum(dev.resids(y, mu_new, weights))
      }
    }
    eta <- eta_new
    mu <- mu_new
    Ainv <- chol2inv(R)
    wt <- W2
    if (abs(dev - devold) / (abs(dev) + 0.1) < control$epsilon) {
      conv <- TRUE
      coef <- start
      break
    }
    devold <- dev
    coef <- start
  }
  if (!conv) warning("fs_glm: algorithm did not converge", call. = FALSE)
  resid <- (y - mu) / mu.eta(eta)
  list(coef = drop(coef), Ainv = Ainv, resid = resid, wt = wt, keep = keep,
       rank = length(keep), iter = iter, converged = conv, dev = dev, mu = mu)
}

# Aliased columns are dropped the way glm's pivoting QR does (later columns
# go first); the QR is only run when the Cholesky pivots suggest a problem.
alias_check <- function(X, W2, A) {
  # all-zero columns (e.g. unused factor levels) are always aliased
  nz <- which(diag(A) > 0)
  B <- A[nz, nz, drop = FALSE]
  d <- diag(B)
  ch <- suppressWarnings(chol(B, pivot = TRUE, tol = -1))
  rk <- attr(ch, "rank")
  piv_min <- min(abs(diag(ch))[seq_len(rk)]^2 / d[attr(ch, "pivot")][seq_len(rk)])
  if (rk == ncol(B) && isTRUE(piv_min > 1e-8)) return(nz)
  q <- qr(X[, nz, drop = FALSE] * sqrt(W2), tol = 1e-11, LAPACK = FALSE)
  nz[sort(q$pivot[seq_len(q$rank)])]
}

glmfit_fallback <- function(X, y, w, offset, family, control) {
  f <- stats::glm.fit(X, y, weights = w, offset = offset, family = family, control = control)
  keep <- f$qr$pivot[seq_len(f$rank)]
  Ainv <- chol2inv(f$qr$qr[seq_len(f$rank), seq_len(f$rank), drop = FALSE])
  ord <- order(keep)
  keep <- keep[ord]
  Ainv <- Ainv[ord, ord, drop = FALSE]
  eta <- f$linear.predictors
  list(coef = f$coefficients[keep], Ainv = Ainv,
       resid = (y - f$fitted.values) / family$mu.eta(eta), wt = f$weights, keep = keep,
       rank = f$rank, iter = f$iter, converged = f$converged, dev = f$deviance,
       mu = f$fitted.values)
}

#' @export
coef.fs_glm <- function(object, na.rm = TRUE, ...) {
  cf <- object$coefficients
  if (na.rm) cf[!is.na(cf)] else cf
}

#' @export
vcov.fs_glm <- function(object, ...) object$vcov

coef_table <- function(object) {
  cf <- object$coefficients[!is.na(object$coefficients)]
  se <- sqrt(diag(object$vcov))
  tv <- cf / se
  df <- object$df.residual
  pv <- if (df > 0) 2 * stats::pt(-abs(tv), df) else NaN
  cbind(Estimate = cf, `Std. Error` = se, `t value` = tv, `Pr(>|t|)` = pv)
}

#' @export
confint.fs_glm <- function(object, parm, level = 0.95, ...) {
  cf <- object$coefficients[!is.na(object$coefficients)]
  se <- sqrt(diag(object$vcov))
  ci <- ci_t(cf, se, level, object$df.residual)
  a <- (1 - level) / 2
  dimnames(ci) <- list(names(cf), paste(format(100 * c(a, 1 - a), trim = TRUE), "%"))
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' @export
summary.fs_glm <- function(object, ...) {
  structure(list(call = object$formula, family = object$family, coefficients = coef_table(object),
                 df.residual = object$df.residual, n = object$n), class = "summary.fs_glm")
}

#' @export
print.summary.fs_glm <- function(x, ...) {
  cat("Survey-weighted GLM (rapidsvy), family:", x$family$family, "link:", x$family$link, "\n")
  cat("Formula:", deparse(x$call), "\n")
  cat("n =", x$n, " design df =", x$df.residual, "\n\n")
  stats::printCoefmat(x$coefficients, ...)
  invisible(x)
}

#' @export
print.fs_glm <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' Tidy an fs_glm fit
#' @param x An `fs_glm` object.
#' @param conf.int Add confidence intervals.
#' @param conf.level Confidence level.
#' @param exponentiate Exponentiate estimates and intervals (odds/rate ratios).
#' @param ... Unused.
#' @export
tidy.fs_glm <- function(x, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE, ...) {
  ct <- coef_table(x)
  out <- tibble::tibble(term = rownames(ct), estimate = ct[, 1], std.error = ct[, 2],
                        statistic = ct[, 3], p.value = ct[, 4])
  if (conf.int) {
    ci <- ci_t(ct[, 1], ct[, 2], conf.level, x$df.residual)
    out$conf.low <- ci[, 1]
    out$conf.high <- ci[, 2]
  }
  if (exponentiate) {
    out$estimate <- exp(out$estimate)
    if (conf.int) {
      out$conf.low <- exp(out$conf.low)
      out$conf.high <- exp(out$conf.high)
    }
  }
  out
}
