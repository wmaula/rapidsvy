# Core estimators shared by the summarise() statistics, fs_tab(), fs_chisq()
# and fs_glm(). `gid` is a 0-based domain index per row (-1 = not in any
# domain), `G` the number of domains.

.fs_ctx <- new.env(parent = emptyenv())

cur_ctx <- function() {
  if (is.null(.fs_ctx$design))
    stop("rapidsvy statistics (fs_mean(), fs_prop(), ...) can only be used inside ",
         "summarise() on an fsvy design.", call. = FALSE)
  .fs_ctx
}

check_lonely <- function(design, res) {
  h <- res$fail_stratum
  if (h > 0) {
    lab <- design$str_label(design$str_first_row[h])
    stop(sprintf(paste0("Stratum (%s) has only one PSU.\n",
                        "Set lonely_psu = \"adjust\", \"average\", \"remove\" or \"certainty\" ",
                        "in as_fsvy() (same meaning as options(survey.lonely.psu))."), lab),
         call. = FALSE)
  }
  res
}

kernel_num <- function(design, Z, gid, G, full = FALSE, mult = NULL, zrow = NULL) {
  Z <- as.matrix(Z)
  if (storage.mode(Z) != "double") storage.mode(Z) <- "double"
  if (anyNA(Z)) Z[is.na(Z)] <- 0
  if (!is.null(mult) && anyNA(mult)) mult[is.na(mult)] <- 0
  res <- fs_var_num(Z, mult, zrow, as.integer(gid), as.integer(G), design$strat_ptr,
                    design$psu_ptr, design$rows, design$fpc_f,
                    lonely_code(design$lonely_psu), full, fs_threads())
  check_lonely(design, res)
}

kernel_cat <- function(design, ycat, gid, G, L, P, Wg, full = FALSE) {
  P[is.na(P)] <- 0
  res <- fs_var_cat(as.integer(ycat), design$w, as.integer(gid), as.integer(G), as.integer(L),
                    P, as.double(Wg), design$strat_ptr, design$psu_ptr, design$rows,
                    design$fpc_f, lonely_code(design$lonely_psu), full, fs_threads())
  check_lonely(design, res)
}

grp_sum <- function(X, gid, G) {
  X <- as.matrix(X)
  if (storage.mode(X) != "double") storage.mode(X) <- "double"
  fs_grp_sum(X, as.integer(gid), as.integer(G))
}

as_num <- function(x, what = "x") {
  if (is.factor(x) || is.character(x))
    stop(sprintf("`%s` is categorical; use fs_prop() or fs_tab() for proportions.", what),
         call. = FALSE)
  if (is.logical(x)) return(as.double(x))
  as.double(vctrs::vec_data(x))
}

# Drop rows with missing values from the domain (na.rm = TRUE) or flag the
# domains that contain them (na.rm = FALSE gives NA, as survey does)
handle_na <- function(gid, G, na_rows, na.rm) {
  bad <- logical(G)
  if (any(na_rows)) {
    hit <- gid[na_rows]
    if (!na.rm) bad[unique(hit[hit >= 0]) + 1L] <- TRUE
    gid[na_rows] <- -1L
  }
  list(gid = gid, bad = bad)
}

# Mean, total and ratio of totals for every domain.
ratio_engine <- function(design, gid, G, y, x = NULL, type = c("mean", "total", "ratio"),
                         na.rm = FALSE, deff = FALSE) {
  type <- match.arg(type)
  w <- design$w
  na_rows <- is.na(y)
  if (!is.null(x)) na_rows <- na_rows | is.na(x)
  h <- handle_na(gid, G, na_rows, na.rm)
  gid <- h$gid
  ok <- gid >= 0L
  gi <- gid + 1L
  y[!ok] <- 0
  den <- if (type == "ratio") { x[!ok] <- 0; x } else rep(1, length(y))

  S <- grp_sum(cbind(w * y, w * den, w, as.double(w != 0)), gid, G)
  Ty <- S[, 1]; Tx <- S[, 2]; W <- S[, 3]; nobs <- S[, 4]
  est <- switch(type, mean = Ty / W, total = Ty, ratio = Ty / Tx)

  z <- numeric(length(y))
  if (type == "total") {
    z[ok] <- w[ok] * y[ok]
  } else {
    z[ok] <- w[ok] * (y[ok] - est[gi[ok]] * den[ok]) / Tx[gi[ok]]
  }
  v <- kernel_num(design, z, gid, G)$V[, 1]

  d <- NULL
  if (!isFALSE(deff)) {
    # survey: svyvar() times the SRS factor; ratios go through svytotal(r)
    u <- switch(type, mean = y, total = y, ratio = (y - est[pmax(gi, 1L)] * den) / Tx[pmax(gi, 1L)])
    u[!ok] <- 0
    ubar <- grp_sum(w * u, gid, G)[, 1] / W
    ss <- numeric(length(u))
    ss[ok] <- w[ok] * (u[ok] - ubar[gi[ok]])^2
    s2 <- grp_sum(ss, gid, G)[, 1] / W * nobs / (nobs - 1)
    fac <- if (identical(deff, "replace")) 1 / nobs else (W - nobs) / (W * nobs)
    vsrs <- if (type == "mean") s2 * fac else s2 * W^2 * fac
    if (!identical(deff, "replace") && any(W < nobs & nobs > 0))
      warning("Sample size greater than population size: are weights correctly scaled?", call. = FALSE)
    vsrs[W < nobs] <- NA
    d <- v / vsrs
    d[h$bad] <- NA
  }
  est[h$bad] <- NA
  v[h$bad] <- NA
  list(est = est, se = sqrt(v), deff = d, W = W, n = nobs)
}

# Proportions of categorical y (0-based codes, -1 = missing) within domains.
prop_engine <- function(design, gid, G, ycat, L, full = FALSE) {
  w <- design$w
  gid[ycat < 0L] <- -1L
  ok <- gid >= 0L
  cell <- ifelse(ok, gid * L + ycat, -1L)
  num <- matrix(grp_sum(w, cell, G * L)[, 1], G, L, byrow = TRUE)
  S <- grp_sum(cbind(w, as.double(w != 0)), gid, G)
  W <- S[, 1]
  P <- num / W
  ycat[!ok] <- 0L
  res <- kernel_cat(design, ycat, gid, G, L, P, W, full = full)
  list(P = P, V = res$V, W = W, n = S[, 2], res = res)
}

ci_t <- function(est, se, level, df) {
  a <- (1 - level) / 2
  fac <- ifelse(is.finite(df), stats::qt(1 - a, df), stats::qnorm(1 - a))
  cbind(est - fac * se, est + fac * se)
}

check_vartype <- function(vartype) {
  if (is.null(vartype)) return(character())
  match.arg(vartype, c("se", "ci", "var", "cv"), several.ok = TRUE)
}

# Assemble srvyr-style result columns; "coef" is renamed by summarise()
make_out <- function(est, se, vartype, level = 0.95, df = Inf, deff = NULL, ci = NULL, suffix = "") {
  out <- list()
  out[[if (nzchar(suffix)) suffix else "coef"]] <- unname(est)
  for (v in vartype) {
    if (v == "se") out[[paste0(suffix, "_se")]] <- unname(se)
    if (v == "ci") {
      if (is.null(ci)) ci <- ci_t(est, se, level, df)
      out[[paste0(suffix, "_low")]] <- unname(ci[, 1])
      out[[paste0(suffix, "_upp")]] <- unname(ci[, 2])
    }
    if (v == "var") out[[paste0(suffix, "_var")]] <- unname(se^2)
    if (v == "cv") out[[paste0(suffix, "_cv")]] <- unname(se / est)
  }
  if (!is.null(deff)) out[[paste0(suffix, "_deff")]] <- unname(deff)
  tibble::new_tibble(out, nrow = length(est))
}

# Confidence intervals for proportions, following survey::svyciprop()
prop_ci <- function(p, se, level, df, method, n_rows = NULL) {
  alpha <- 1 - level
  tq <- ifelse(is.finite(df), stats::qt(1 - alpha / 2, df), stats::qnorm(1 - alpha / 2))
  ci <- switch(method,
    mean = cbind(p - tq * se, p + tq * se),
    logit = , xlogit = {
      eta <- stats::qlogis(p)
      s <- se / (p * (1 - p))
      cbind(stats::plogis(eta - tq * s), stats::plogis(eta + tq * s))
    },
    asin = {
      eta <- asin(sqrt(p))
      s <- se / (2 * sqrt(p * (1 - p)))
      cbind(sin(eta - tq * s)^2, sin(eta + tq * s)^2)
    },
    beta = {
      neff <- p * (1 - p) / se^2
      neff <- neff * (stats::qt(alpha / 2, n_rows - 1) / stats::qt(alpha / 2, df))^2
      cbind(stats::qbeta(alpha / 2, neff * p, neff * (1 - p) + 1),
            stats::qbeta(1 - alpha / 2, neff * p + 1, neff * (1 - p)))
    },
    wilson = {
      neff <- p * (1 - p) / se^2
      den <- 1 + tq^2 / neff
      rt <- sqrt(4 * neff * p * (1 - p) + tq^2)
      cbind((p + tq^2 / (2 * neff) - tq / (2 * neff) * rt) / den,
            (p + tq^2 / (2 * neff) + tq / (2 * neff) * rt) / den)
    },
    stop("Unknown ci_method: ", method, call. = FALSE))
  # degenerate cells (p = 0 or 1, or zero variance): interval collapses to p
  deg <- !is.na(p) & (p <= 0 | p >= 1 | se == 0)
  ci[deg, ] <- p[deg]
  ci
}

quantile_suffix <- function(q) {
  paste0("_q", gsub(".", "", formatC(q * 100, width = 2, flag = "0"), fixed = TRUE))
}

# Weighted quantiles with Woodruff confidence intervals (survey::svyquantile
# with qrule = "math", interval.type = "mean" or "xlogit")
quantile_engine <- function(design, gid, G, x, probs, na.rm = FALSE, level = 0.95,
                            interval_type = "mean", df = NULL) {
  w <- design$w
  h <- handle_na(gid, G, is.na(x), na.rm)
  gid <- h$gid
  ok <- gid >= 0L
  Q <- length(probs)

  sel <- which(ok & w != 0)
  o <- order(gid[sel], x[sel], method = "radix")
  idx <- sel[o]
  ptr <- c(0L, cumsum(tabulate(gid[idx] + 1L, G)))
  xs <- x[idx]
  ws <- w[idx]
  qhat <- fs_wquantile(xs, ws, as.integer(ptr), matrix(probs, G, Q, byrow = TRUE))

  # estimated CDF at each quantile, and its linearised variance
  gi <- gid + 1L
  Wg <- grp_sum(w, gid, G)[, 1]
  Z <- matrix(0, length(x), Q)
  phat <- matrix(NA_real_, G, Q)
  for (j in seq_len(Q)) {
    ind <- numeric(length(x))
    ind[ok] <- as.double(x[ok] <= qhat[gi[ok], j])
    phat[, j] <- grp_sum(w * ind, gid, G)[, 1] / Wg
    Z[ok, j] <- w[ok] * (ind[ok] - phat[gi[ok], j]) / Wg[gi[ok]]
  }
  kr <- kernel_num(design, Z, gid, G)
  sep <- sqrt(kr$V)
  if (is.null(df)) df <- kr$npsu - kr$nstrata
  df <- rep_len(df, G)
  alpha <- round(1 - level, 7)
  tq <- ifelse(is.finite(df), stats::qt(1 - alpha / 2, df), stats::qnorm(1 - alpha / 2))

  plo <- phat - tq * sep
  pup <- phat + tq * sep
  if (interval_type == "xlogit") {
    eta <- stats::qlogis(phat)
    s <- sep / (phat * (1 - phat))
    plo <- stats::plogis(eta - tq * s)
    pup <- stats::plogis(eta + tq * s)
  }
  plo_q <- ifelse(is.nan(plo) | plo < 0, NA, plo)
  pup_q <- ifelse(is.nan(pup) | pup > 1, NA, pup)
  lo <- fs_wquantile(xs, ws, as.integer(ptr), plo_q)
  up <- fs_wquantile(xs, ws, as.integer(ptr), pup_q)
  lo[is.na(plo_q)] <- NaN
  up[is.na(pup_q)] <- NaN
  se <- (up - lo) / (2 * tq)

  bad <- h$bad
  qhat[bad, ] <- NA; lo[bad, ] <- NA; up[bad, ] <- NA; se[bad, ] <- NA
  list(q = qhat, se = se, low = lo, upp = up)
}
