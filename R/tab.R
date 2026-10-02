#' Weighted frequency table with proportions
#'
#' One call for the usual descriptive table: unweighted n, estimated
#' population count, and proportion of each level of `var` within the
#' groups in `by` (or the current `group_by()` groups).
#'
#' @param .data An `fsvy` design.
#' @param var Categorical variable to tabulate.
#' @param by Grouping columns (tidyselect). Default: current groups.
#' @param vartype,level,ci_method,deff,df See [fs_stats].
#' @param percent Multiply proportions (and their intervals) by 100.
#' @param na.rm Drop rows with missing `var` before tabulating.
#' @return A tibble.
#' @examples
#' set.seed(1)
#' d <- data.frame(psu = rep(1:40, each = 5), str = rep(1:8, each = 25),
#'                 w = runif(200, 1, 3), sex = sample(c("F", "M"), 200, TRUE),
#'                 reg = sample(c("A", "B"), 200, TRUE))
#' des <- as_fsvy(d, ids = psu, strata = str, weights = w)
#' fs_tab(des, sex, by = reg)
#' @export
fs_tab <- function(.data, var, by = NULL, vartype = c("se", "ci"), level = 0.95,
                   ci_method = c("logit", "mean", "xlogit", "asin", "beta", "wilson"),
                   deff = FALSE, df = NULL, percent = FALSE, na.rm = FALSE) {
  ci_method <- match.arg(ci_method)
  var <- names(tidyselect::eval_select(enquo(var), .data$data))
  byq <- enquo(by)
  by <- if (quo_is_null(byq)) .data$groups else names(tidyselect::eval_select(byq, .data$data))
  if (na.rm) .data <- filter.fsvy(.data, !is.na(.data[[var]]))
  .data <- group_by.fsvy(ungroup.fsvy(.data), dplyr::across(dplyr::all_of(c(by, var))))
  out <- summarise.fsvy(.data,
                        n = fs_n(),
                        N = fs_total(vartype = NULL),
                        prop = fs_prop(vartype = vartype, level = level, ci_method = ci_method,
                                       deff = deff, df = df))
  if (percent) {
    pc <- intersect(c("prop", "prop_se", "prop_low", "prop_upp"), names(out))
    out[pc] <- lapply(out[pc], function(v) 100 * v)
    out$prop_var <- if ("prop_var" %in% names(out)) out$prop_var * 1e4 else NULL
  }
  out
}

#' Rao-Scott test of association for a two-way table
#'
#' Same statistic as `survey::svychisq()` with `statistic = "F"` (default)
#' or `"Chisq"`. Uses the current domain (`filter()`); groups are ignored.
#'
#' @param .data An `fsvy` design.
#' @param row,col Categorical variables.
#' @param statistic "F" (second-order Rao-Scott) or "Chisq" (first-order).
#' @return A one-row tibble with the statistic, degrees of freedom and p-value.
#' @export
fs_chisq <- function(.data, row, col, statistic = c("F", "Chisq")) {
  statistic <- match.arg(statistic)
  design <- .data
  rv <- names(tidyselect::eval_select(enquo(row), design$data))
  cv <- names(tidyselect::eval_select(enquo(col), design$data))
  dom <- domain_mask(design)
  r <- design$data[[rv]]
  cc <- design$data[[cv]]
  N <- sum(dom)
  nu <- fs_degf(design)

  ok <- dom & !is.na(r) & !is.na(cc)
  rl <- sort(unique(r[ok]))
  cl <- sort(unique(cc[ok]))
  nr <- length(rl)
  nc <- length(cl)
  if (nr < 2 || nc < 2) stop("Both variables need at least two observed levels.", call. = FALSE)
  # cells ordered with rows varying fastest, as interaction(rows, cols)
  ycat <- rep(-1L, length(r))
  ycat[ok] <- (match(r[ok], rl) - 1L) + nr * (match(cc[ok], cl) - 1L)
  gid <- ifelse(dom, 0L, -1L)
  pe <- prop_engine(design, gid, 1L, ycat, nr * nc, full = TRUE)
  mean2 <- pe$P[1, ]
  V <- matrix(pe$V[1, ], nr * nc, nr * nc)

  mf1 <- expand.grid(rows = 1:nr, cols = 1:nc)
  X1 <- stats::model.matrix(~ factor(rows) + factor(cols), mf1)
  X12 <- stats::model.matrix(~ factor(rows) * factor(cols), mf1)
  Cmat <- qr.resid(qr(X1), X12[, -(1:(nr + nc - 1)), drop = FALSE])
  iD <- diag(ifelse(mean2 == 0, 0, 1 / mean2))
  denom <- t(Cmat) %*% (iD / N) %*% Cmat
  numr <- t(Cmat) %*% iD %*% V %*% iD %*% Cmat
  Delta <- solve(denom, numr)
  d0 <- sum(diag(Delta))^2 / sum(diag(Delta %*% Delta))

  tab <- matrix(mean2 * N, nr, nc)
  E <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  X2 <- sum((tab - E)^2 / E)

  if (statistic == "F") {
    stat <- X2 / sum(diag(Delta))
    tibble::tibble(statistic = stat, ndf = d0, ddf = d0 * nu,
                   p.value = stats::pf(stat, d0, d0 * nu, lower.tail = FALSE),
                   method = "Pearson's X^2: Rao & Scott adjustment (F)")
  } else {
    # survey reports the unadjusted X^2; the p-value uses X^2 / mean(eigenvalues)
    stat <- X2 / mean(diag(Delta))
    tibble::tibble(statistic = X2, adj_statistic = stat, df = NCOL(Delta),
                   p.value = stats::pchisq(stat, NCOL(Delta), lower.tail = FALSE),
                   method = "Pearson's X^2: Rao & Scott adjustment (Chisq)")
  }
}
