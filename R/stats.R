#' Design-based statistics inside summarise()
#'
#' Use these inside `summarise()` on an [as_fsvy()] design, optionally after
#' `group_by()` and `filter()`. Each grouping combination is a domain;
#' variances use the full design (PSUs outside the domain count as zeros), as
#' in `survey::svyby()` and `srvyr`.
#'
#' * `fs_mean()`, `fs_total()`, `fs_ratio()`: `svymean()`, `svytotal()`,
#'   `svyratio()`.
#' * `fs_prop()`: proportion of each level of the **last** grouping variable
#'   within the other grouping variables (srvyr's `survey_prop()`).
#' * `fs_quantile()`, `fs_median()`: `svyquantile()` with `qrule = "math"`
#'   and Woodruff intervals.
#' * `fs_n()`: unweighted number of rows.
#'
#' @param x,numerator,denominator Numeric (or logical) variables.
#' @param na.rm Drop missing values from the domain. With `FALSE`, domains
#'   containing missing values return `NA`.
#' @param vartype Any of "se", "ci", "var", "cv"; `NULL` for estimates only.
#' @param level Confidence level.
#' @param deff `TRUE` for design effects (`"replace"` for the with-replacement
#'   SRS reference), as in survey.
#' @param df Degrees of freedom for intervals. Default: `fs_degf()` of the
#'   (filtered) design, as srvyr. Use `Inf` for normal intervals.
#' @param ci_method Interval for proportions, as `svyciprop()`: "logit"
#'   (default, as srvyr's `survey_prop()`), "mean", "xlogit", "asin", "beta",
#'   "wilson".
#' @param quantiles Probabilities.
#' @param interval_type "mean" (default) or "xlogit" Woodruff interval.
#' @return A tibble with one row per group; `summarise()` names the columns
#'   `name`, `name_se`, `name_low`, `name_upp`, `name_var`, `name_cv`,
#'   `name_deff`.
#' @name fs_stats
#' @examples
#' set.seed(1)
#' d <- data.frame(psu = rep(1:40, each = 5), str = rep(1:8, each = 25),
#'                 w = runif(200, 1, 3), y = rnorm(200),
#'                 sex = sample(c("F", "M"), 200, TRUE),
#'                 reg = sample(c("A", "B"), 200, TRUE))
#' des <- as_fsvy(d, ids = psu, strata = str, weights = w)
#' des |>
#'   dplyr::group_by(reg) |>
#'   dplyr::summarise(m = fs_mean(y, vartype = c("se", "ci")),
#'                    med = fs_median(y), n = fs_n())
#' des |> dplyr::group_by(reg, sex) |> dplyr::summarise(p = fs_prop(vartype = "ci"))
NULL

#' @rdname fs_stats
#' @export
fs_mean <- function(x, na.rm = FALSE, vartype = "se", level = 0.95, deff = FALSE, df = NULL) {
  ctx <- cur_ctx()
  r <- ratio_engine(ctx$design, ctx$gid, ctx$G, as_num(x), type = "mean", na.rm = na.rm, deff = deff)
  make_out(r$est, r$se, check_vartype(vartype), level, df %||% ctx$df, r$deff)
}

#' @rdname fs_stats
#' @export
fs_total <- function(x, na.rm = FALSE, vartype = "se", level = 0.95, deff = FALSE, df = NULL) {
  ctx <- cur_ctx()
  y <- if (missing(x)) rep(1, nrow(ctx$design$data)) else as_num(x)
  r <- ratio_engine(ctx$design, ctx$gid, ctx$G, y, type = "total", na.rm = na.rm, deff = deff)
  make_out(r$est, r$se, check_vartype(vartype), level, df %||% ctx$df, r$deff)
}

#' @rdname fs_stats
#' @export
fs_ratio <- function(numerator, denominator, na.rm = FALSE, vartype = "se", level = 0.95,
                     deff = FALSE, df = NULL) {
  ctx <- cur_ctx()
  r <- ratio_engine(ctx$design, ctx$gid, ctx$G, as_num(numerator, "numerator"),
                    as_num(denominator, "denominator"), type = "ratio", na.rm = na.rm, deff = deff)
  make_out(r$est, r$se, check_vartype(vartype), level, df %||% ctx$df, r$deff)
}

#' @rdname fs_stats
#' @export
fs_prop <- function(vartype = "se", level = 0.95,
                    ci_method = c("logit", "mean", "xlogit", "asin", "beta", "wilson"),
                    deff = FALSE, df = NULL) {
  ctx <- cur_ctx()
  ci_method <- match.arg(ci_method)
  gv <- ctx$design$groups
  if (!length(gv))
    stop("fs_prop() gives proportions of the last grouping variable; call group_by() first.",
         call. = FALSE)
  keys <- ctx$keys
  last <- gv[length(gv)]
  pid <- if (length(gv) > 1) vctrs::vec_group_id(keys[gv[-length(gv)]]) else rep(1L, ctx$G)
  cid <- vctrs::vec_group_id(keys[last])
  Gp <- max(pid)
  L <- max(cid)
  gid <- ctx$gid
  ok <- gid >= 0L
  row_p <- rep(-1L, length(gid))
  row_c <- rep(0L, length(gid))
  row_p[ok] <- pid[gid[ok] + 1L] - 1L
  row_c[ok] <- cid[gid[ok] + 1L] - 1L

  pe <- prop_engine(ctx$design, row_p, Gp, row_c, L)
  at <- cbind(pid, cid)
  p <- pe$P[at]
  se <- sqrt(pe$V[at])
  df <- df %||% ctx$df
  vartype <- check_vartype(vartype)
  ci <- if ("ci" %in% vartype) prop_ci(p, se, level, df, ci_method, n_rows = ctx$n_domain) else NULL
  d <- NULL
  if (!isFALSE(deff)) {
    n <- pe$n[pid]
    W <- pe$W[pid]
    s2 <- p * (1 - p) * n / (n - 1)
    vsrs <- if (identical(deff, "replace")) s2 / n else s2 * (W - n) / (W * n)
    d <- se^2 / vsrs
  }
  make_out(p, se, vartype, level, df, d, ci)
}

#' @rdname fs_stats
#' @export
fs_quantile <- function(x, quantiles, na.rm = FALSE, vartype = "se", level = 0.95,
                        interval_type = c("mean", "xlogit"), df = NULL) {
  ctx <- cur_ctx()
  interval_type <- match.arg(interval_type)
  q <- quantile_engine(ctx$design, ctx$gid, ctx$G, as_num(x), quantiles, na.rm = na.rm,
                       level = level, interval_type = interval_type, df = df)
  vartype <- check_vartype(vartype)
  outs <- lapply(seq_along(quantiles), function(j) {
    make_out(q$q[, j], q$se[, j], vartype, level, ci = cbind(q$low[, j], q$upp[, j]),
             suffix = quantile_suffix(quantiles[j]))
  })
  tibble::as_tibble(do.call(c, lapply(outs, as.list)))
}

#' @rdname fs_stats
#' @export
fs_median <- function(x, na.rm = FALSE, vartype = "se", level = 0.95,
                      interval_type = c("mean", "xlogit"), df = NULL) {
  out <- fs_quantile(x, 0.5, na.rm = na.rm, vartype = vartype, level = level,
                     interval_type = match.arg(interval_type), df = df)
  names(out) <- sub("^_q50", "", names(out))
  names(out)[names(out) == ""] <- "coef"
  out
}

#' @rdname fs_stats
#' @export
fs_n <- function() {
  ctx <- cur_ctx()
  g <- ctx$gid
  tabulate(g[g >= 0L] + 1L, ctx$G)
}
