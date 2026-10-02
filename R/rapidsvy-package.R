#' rapidsvy: rapid design-based estimation for complex surveys
#'
#' Taylor linearisation for stratified one-stage cluster designs (PSUs
#' sampled with replacement, optional finite population correction). All
#' domains of a `group_by()` are estimated in one multithreaded pass, so the
#' cost does not grow with the number of groups the way `survey::svyby()` and
#' `srvyr` do.
#'
#' Number of threads: `options(rapidsvy.threads = n)`; default is all cores
#' minus one.
#'
#' @keywords internal
#' @useDynLib rapidsvy, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang .data := enquo enquos quo_is_null eval_tidy as_label
#' @importFrom generics tidy
"_PACKAGE"

fs_threads <- function() {
  n <- getOption("rapidsvy.threads")
  if (is.null(n)) {
    n <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) 1L)
    if (is.na(n)) n <- 1L
    n <- max(1L, n - 1L)
  }
  as.integer(max(1L, n))
}

lonely_code <- function(x) {
  match(x, c("fail", "remove", "certainty", "adjust", "average")) - 1L
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' @export
generics::tidy
