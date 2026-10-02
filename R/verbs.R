#' dplyr verbs for fsvy designs
#'
#' `filter()` defines a domain: rows outside it are kept in the design so
#' variances stay correct (as `subset()` in survey). `group_by()` sets the
#' domains estimated by `summarise()`. `mutate()`, `select()` and `rename()`
#' work on the data and never change the design.
#'
#' @param .data,x An `fsvy` object.
#' @param ... Passed to the dplyr verb.
#' @param .add,.groups,.preserve,.by As in dplyr.
#' @name fsvy-verbs
NULL

#' @rdname fsvy-verbs
#' @importFrom dplyr summarise
#' @export
summarise.fsvy <- function(.data, ..., .by = NULL, .groups = NULL) {
  design <- .data
  by <- enquo(.by)
  if (!quo_is_null(by)) design <- group_by.fsvy(design, dplyr::across(!!by))
  dots <- enquos(..., .named = TRUE)
  n <- nrow(design$data)
  dom <- domain_mask(design)
  gv <- design$groups

  if (length(gv)) {
    gdf <- design$data[gv]
    keys <- vctrs::vec_unique(gdf[dom, , drop = FALSE])
    keys <- keys[do.call(order, c(unname(as.list(keys)), list(method = "radix"))), , drop = FALSE]
    gid <- vctrs::vec_match(gdf, keys) - 1L
    gid[is.na(gid) | !dom] <- -1L
  } else {
    keys <- tibble::new_tibble(list(), nrow = 1L)
    gid <- ifelse(dom, 0L, -1L)
  }
  G <- nrow(keys)
  if (G == 0) stop("No rows in the current domain.", call. = FALSE)

  old <- as.list(.fs_ctx)
  on.exit({
    rm(list = ls(.fs_ctx), envir = .fs_ctx)
    if (length(old)) list2env(old, .fs_ctx)
  }, add = TRUE)
  .fs_ctx$design <- design
  .fs_ctx$gid <- as.integer(gid)
  .fs_ctx$G <- G
  .fs_ctx$keys <- keys
  .fs_ctx$df <- fs_degf(design)
  .fs_ctx$n_domain <- sum(dom)

  mask <- rlang::as_data_mask(design$data)
  res <- list()
  for (nm in names(dots)) {
    val <- eval_tidy(dots[[nm]], data = mask)
    if (is.data.frame(val)) {
      if (nrow(val) != G) stop(sprintf("`%s` returned %d rows, expected %d.", nm, nrow(val), G), call. = FALSE)
      cn <- names(val)
      cn <- ifelse(cn == "coef", nm, paste0(nm, cn))
      res[cn] <- as.list(val)
    } else if (length(val) %in% c(1L, G)) {
      res[[nm]] <- rep_len(val, G)
    } else {
      stop(sprintf("`%s` must be an fs_*() statistic or have length 1 or %d.", nm, G), call. = FALSE)
    }
  }
  out <- tibble::as_tibble(c(as.list(keys), res))
  if (!is.null(.groups) && length(gv)) {
    keep <- switch(.groups, drop = character(), drop_last = gv[-length(gv)], keep = gv,
                   stop("Unsupported .groups", call. = FALSE))
    if (length(keep)) out <- dplyr::group_by(out, dplyr::across(dplyr::all_of(keep)))
  }
  out
}

#' @rdname fsvy-verbs
#' @importFrom dplyr group_by
#' @export
group_by.fsvy <- function(.data, ..., .add = FALSE) {
  d <- .data$data
  if (.add && length(.data$groups)) d <- dplyr::group_by(d, dplyr::across(dplyr::all_of(.data$groups)))
  d <- dplyr::group_by(d, ..., .add = .add)
  .data$groups <- dplyr::group_vars(d)
  .data$data <- dplyr::ungroup(d)
  .data
}

#' @rdname fsvy-verbs
#' @importFrom dplyr ungroup
#' @export
ungroup.fsvy <- function(x, ...) {
  x$groups <- character()
  x
}

#' @rdname fsvy-verbs
#' @importFrom dplyr group_vars
#' @export
group_vars.fsvy <- function(x) x$groups

#' @rdname fsvy-verbs
#' @importFrom dplyr filter
#' @export
filter.fsvy <- function(.data, ..., .preserve = FALSE) {
  dots <- enquos(...)
  mask <- rlang::as_data_mask(.data$data)
  n <- nrow(.data$data)
  keep <- domain_mask(.data)
  for (q in dots) {
    v <- eval_tidy(q, data = mask)
    if (!is.logical(v)) stop("filter() conditions must be logical.", call. = FALSE)
    keep <- keep & !is.na(v) & rep_len(v, n)
  }
  .data$domain <- keep
  .data
}

#' @rdname fsvy-verbs
#' @importFrom dplyr mutate
#' @export
mutate.fsvy <- function(.data, ...) {
  d <- .data$data
  if (length(.data$groups)) d <- dplyr::group_by(d, dplyr::across(dplyr::all_of(.data$groups)))
  d <- dplyr::ungroup(dplyr::mutate(d, ...))
  if (nrow(d) != nrow(.data$data)) stop("mutate() must not change the number of rows.", call. = FALSE)
  .data$data <- d
  .data
}

#' @rdname fsvy-verbs
#' @importFrom dplyr select
#' @export
select.fsvy <- function(.data, ...) {
  pos <- tidyselect::eval_select(rlang::expr(c(...)), .data$data)
  d <- .data$data[pos]
  names(d) <- names(pos)
  miss <- setdiff(.data$groups, names(d))
  if (length(miss)) d <- dplyr::bind_cols(.data$data[miss], d)
  .data$data <- d
  .data
}

#' @rdname fsvy-verbs
#' @importFrom dplyr rename
#' @export
rename.fsvy <- function(.data, ...) {
  pos <- tidyselect::eval_rename(rlang::expr(c(...)), .data$data)
  old <- names(.data$data)[pos]
  names(.data$data)[pos] <- names(pos)
  hit <- .data$groups %in% old
  .data$groups[hit] <- names(pos)[match(.data$groups[hit], old)]
  .data
}

#' @importFrom dplyr pull
#' @export
pull.fsvy <- function(.data, var = -1, name = NULL, ...) {
  dplyr::pull(.data$data, {{ var }}, {{ name }}, ...)
}
