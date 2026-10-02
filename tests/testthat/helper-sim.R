sim_data <- function(n_str = 30, seed = 1, lonely = 3) {
  set.seed(seed)
  # PSUs per stratum vary; a few strata have a single PSU
  npsu <- sample(2:6, n_str, replace = TRUE)
  npsu[seq_len(lonely)] <- 1L
  psu <- unlist(lapply(seq_len(n_str), function(h) paste0(h, "_", seq_len(npsu[h]))))
  str <- as.integer(sub("_.*", "", psu))
  size <- sample(1:5, length(psu), replace = TRUE)
  d <- data.frame(str = rep(str, size), psu = rep(psu, size))
  n <- nrow(d)
  d$w <- round(runif(n, 5, 200), 2)
  d$reg <- sample(c("A", "B", "C", "D"), n, TRUE)
  d$sex <- factor(sample(c("F", "M"), n, TRUE))
  d$cat <- sample(c("x", "y", "z"), n, TRUE, prob = c(.6, .3, .1))
  d$y <- rgamma(n, 2, 0.01)
  d$x <- d$y * runif(n, 0.2, 1) + rnorm(n, 50, 10)
  d$b <- rbinom(n, 1, plogis(-1 + 0.004 * d$y))
  d$y_na <- ifelse(runif(n) < 0.05, NA, d$y)
  d$npop <- npsu[d$str] * 3
  d
}
