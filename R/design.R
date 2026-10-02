#' Create a fast survey design object
#'
#' Declares a stratified one-stage cluster design (PSUs with replacement),
#' the same design as `survey::svydesign(ids = ~psu, strata = ~strata,
#' weights = ~w, nest = ...)` or `srvyr::as_survey_design()`.
#'
#' @param .data A data frame, a `survey.design2` object, or a `tbl_svy`.
#' @param ids PSU (cluster) column(s), tidyselect or one-sided formula. Use
#'   `NULL` or `1` when each row is its own PSU.
#' @param strata Stratum column(s); several columns are crossed.
#' @param weights Sampling weight column.
#' @param fpc Optional finite population correction: population number of
#'   PSUs per stratum, or the sampling fraction (values <= 1).
#' @param nest If `TRUE`, PSU ids are only unique within strata.
#' @param lonely_psu Handling of strata with a single PSU, as
#'   `options(survey.lonely.psu)`: "fail" (survey default), "adjust",
#'   "average", "remove" or "certainty".
#' @param ... Unused.
#' @return An object of class `fsvy`.
#' @examples
#' d <- data.frame(psu = rep(1:20, each = 5), str = rep(1:4, each = 25),
#'                 w = runif(100, 1, 3), y = rnorm(100))
#' des <- as_fsvy(d, ids = psu, strata = str, weights = w)
#' des
#' @export
as_fsvy <- function(.data, ...) UseMethod("as_fsvy")

#' @rdname as_fsvy
#' @export
as_fsvy.data.frame <- function(.data, ids = NULL, strata = NULL, weights = NULL,
                               fpc = NULL, nest = FALSE,
                               lonely_psu = getOption("rapidsvy.lonely_psu", "fail"),
                               ...) {
  lonely_psu <- match.arg(lonely_psu, c("fail", "adjust", "average", "remove", "certainty"))
  data <- tibble::as_tibble(.data)
  n <- nrow(data)
  if (n == 0) stop("data has no rows", call. = FALSE)

  ids_cols <- select_cols(enquo(ids), data)
  str_cols <- select_cols(enquo(strata), data)
  w_cols <- select_cols(enquo(weights), data)
  fpc_cols <- select_cols(enquo(fpc), data)

  w <- if (length(w_cols)) as.double(vctrs::vec_data(data[[w_cols[1]]])) else rep(1, n)
  if (anyNA(w)) stop("weights contain missing values", call. = FALSE)
  if (any(w < 0)) stop("weights must be non-negative", call. = FALSE)

  strat_raw <- if (length(str_cols)) vctrs::vec_group_id(data[str_cols]) else rep(1L, n)
  if (length(str_cols) && anyNA(data[str_cols])) stop("strata contain missing values", call. = FALSE)
  if (length(ids_cols)) {
    if (anyNA(data[ids_cols])) stop("PSU ids contain missing values", call. = FALSE)
    psu_raw <- vctrs::vec_group_id(data[ids_cols])
  } else {
    psu_raw <- seq_len(n)
  }
  fpc_raw <- if (length(fpc_cols)) as.double(vctrs::vec_data(data[[fpc_cols[1]]])) else NULL

  str_label <- if (length(str_cols)) {
    lab <- data[str_cols]
    function(h) paste(vapply(lab[h, , drop = FALSE], function(v) as.character(v[1]), ""), collapse = "/")
  } else {
    function(h) "1"
  }

  structure(
    c(build_design(psu_raw, strat_raw, w, fpc_raw, nest),
      list(data = data, domain = NULL, groups = character(),
           lonely_psu = lonely_psu, str_label = str_label,
           vars = list(ids = ids_cols, strata = str_cols, weights = w_cols, fpc = fpc_cols))),
    class = "fsvy"
  )
}

#' @rdname as_fsvy
#' @export
as_fsvy.survey.design2 <- function(.data, lonely_psu = getOption("survey.lonely.psu", "fail"), ...) {
  if (NCOL(.data$cluster) > 1 && !isTRUE(getOption("survey.ultimate.cluster")))
    message("Only the first stage is used (ultimate cluster variance), as in rapidsvy designs.")
  if (!is.null(.data$postStrata)) stop("Post-stratified or calibrated designs are not supported yet", call. = FALSE)
  data <- tibble::as_tibble(.data$variables)
  w <- 1 / as.vector(.data$prob)
  dom <- is.finite(w) & w > 0
  w[!is.finite(w)] <- 0
  psu_raw <- vctrs::vec_group_id(.data$cluster[[1]])
  strat_raw <- vctrs::vec_group_id(.data$strata[[1]])
  popsize <- .data$fpc$popsize
  fpc_raw <- if (is.null(popsize)) NULL else as.vector(popsize[, 1])
  lonely_psu <- match.arg(lonely_psu, c("fail", "adjust", "average", "remove", "certainty"))
  # survey already made PSU ids unique within strata when nest = TRUE
  pairs <- vctrs::vec_unique_count(data.frame(strat_raw, psu_raw))
  nest <- pairs != vctrs::vec_unique_count(psu_raw)
  strata_vec <- .data$strata[[1]]
  des <- build_design(psu_raw, strat_raw, w, fpc_raw, nest)
  sampsize <- .data$fpc$sampsize
  if (!is.null(sampsize) && any(as.vector(sampsize[, 1]) != diff(des$strat_ptr)[des$strat_id]))
    stop("This design was subset with survey's `[` or subset(), which drops PSUs.\n",
         "Convert the full design and use filter() on the fsvy object instead.", call. = FALSE)
  structure(
    c(des,
      list(data = data, domain = if (all(dom)) NULL else dom, groups = character(),
           lonely_psu = lonely_psu, str_label = function(h) as.character(strata_vec[h]),
           vars = list())),
    class = "fsvy"
  )
}

#' @rdname as_fsvy
#' @export
as_fsvy.tbl_svy <- function(.data, ...) {
  out <- as_fsvy.survey.design2(.data, ...)
  out$data <- tibble::as_tibble(.data$variables)
  g <- dplyr::group_vars(.data)
  if (length(g)) out$groups <- g
  out
}

select_cols <- function(q, data) {
  if (quo_is_null(q)) return(character())
  ex <- rlang::quo_get_expr(q)
  if (is.numeric(ex) && length(ex) == 1 && ex == 1) return(character())
  if (rlang::is_formula(ex)) {
    v <- all.vars(ex)
    if (!length(v)) return(character())
    miss <- setdiff(v, names(data))
    if (length(miss)) stop("Column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
    return(v)
  }
  names(tidyselect::eval_select(q, data))
}

# Sort rows by stratum then PSU and record CSR style offsets for the kernel
build_design <- function(psu_raw, strat_raw, w, fpc_raw, nest) {
  n <- length(psu_raw)
  if (nest) {
    psu_raw <- vctrs::vec_group_id(data.frame(s = strat_raw, p = psu_raw))
  } else {
    n_pairs <- vctrs::vec_unique_count(data.frame(s = strat_raw, p = psu_raw))
    if (n_pairs != vctrs::vec_unique_count(psu_raw))
      stop("Clusters not nested in strata at top level; you may want nest = TRUE.", call. = FALSE)
  }
  ord <- order(strat_raw, psu_raw, method = "radix")
  hs <- strat_raw[ord]
  ps <- psu_raw[ord]
  new_psu <- c(TRUE, hs[-1L] != hs[-n] | ps[-1L] != ps[-n])
  psu_start <- which(new_psu)
  C <- length(psu_start)
  psu_ptr <- c(psu_start - 1L, n)
  psu_strat <- hs[psu_start]
  new_str <- c(TRUE, psu_strat[-1L] != psu_strat[-C])
  str_start <- which(new_str)
  H <- length(str_start)
  strat_ptr <- c(str_start - 1L, C)

  # per-row PSU and stratum index (1-based, in sorted order)
  psu_id <- integer(n)
  psu_id[ord] <- cumsum(new_psu)
  str_of_psu <- cumsum(new_str)
  strat_id <- str_of_psu[psu_id]
  n_h <- diff(strat_ptr)
  # first row of each stratum, used for fpc and error messages
  str_first_row <- ord[psu_ptr[str_start] + 1L]

  fpc_f <- rep(1, H)
  if (!is.null(fpc_raw)) {
    if (anyNA(fpc_raw)) stop("fpc contains missing values", call. = FALSE)
    if (all(fpc_raw <= 1)) {
      popsize <- n_h[strat_id] / fpc_raw
    } else {
      popsize <- fpc_raw
    }
    pop_h <- popsize[str_first_row]
    if (vctrs::vec_unique_count(data.frame(strat_id, popsize)) != H)
      stop("fpc must be constant within strata", call. = FALSE)
    if (any(pop_h < n_h)) stop("fpc smaller than the number of sampled PSUs", call. = FALSE)
    fpc_f <- ifelse(is.infinite(pop_h), 1, (pop_h - n_h) / pop_h)
  }
  list(w = w, rows = ord - 1L, psu_ptr = as.integer(psu_ptr), strat_ptr = as.integer(strat_ptr),
       fpc_f = fpc_f, psu_id = psu_id, strat_id = strat_id, n_psu = C, n_strata = H,
       str_first_row = str_first_row, has_fpc = !is.null(fpc_raw))
}

#' @export
print.fsvy <- function(x, ...) {
  dom <- x$domain
  cat("<fsvy> Stratified 1-stage cluster design (with replacement)\n")
  cat(sprintf("  rows: %s | PSUs: %s | strata: %s | lonely PSU: %s%s\n",
              format(nrow(x$data), big.mark = ","), format(x$n_psu, big.mark = ","),
              format(x$n_strata, big.mark = ","), x$lonely_psu,
              if (x$has_fpc) " | fpc" else ""))
  lonely <- sum(diff(x$strat_ptr) == 1L)
  if (lonely) cat(sprintf("  strata with a single PSU: %s\n", format(lonely, big.mark = ",")))
  if (!is.null(dom)) cat(sprintf("  domain (filter): %s rows\n", format(sum(dom), big.mark = ",")))
  if (length(x$groups)) cat("  groups:", paste(x$groups, collapse = ", "), "\n")
  v <- x$vars
  if (length(v)) {
    cat("  ids:", if (length(v$ids)) paste(v$ids, collapse = " + ") else "1",
        "| strata:", if (length(v$strata)) paste(v$strata, collapse = " x ") else "none",
        "| weights:", if (length(v$weights)) v$weights else "none", "\n")
  }
  cat("Data:\n")
  print(x$data, n = 5)
  invisible(x)
}

#' Degrees of freedom of a design
#'
#' Number of PSUs minus number of strata among rows in the current domain,
#' as `survey::degf()`.
#' @param design An `fsvy` object.
#' @export
fs_degf <- function(design) {
  dom <- design$domain
  if (is.null(dom)) dom <- design$w != 0 else dom <- dom & design$w != 0
  sum(tabulate(design$psu_id[dom], design$n_psu) > 0) -
    sum(tabulate(design$strat_id[dom], design$n_strata) > 0)
}

#' Sampling weights of a design
#' @param design An `fsvy` object.
#' @export
fs_weights <- function(design) design$w

domain_mask <- function(design) {
  dom <- design$domain
  if (is.null(dom)) rep(TRUE, nrow(design$data)) else dom
}

#' @export
dim.fsvy <- function(x) dim(x$data)
