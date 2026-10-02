skip_if_not_installed("survey")
library(survey)
library(dplyr)

test_that("fs_tab equals fs_n, fs_total and fs_prop, and matches svytable", {
  d <- sim_data(lonely = 0)
  sd <- svydesign(ids = ~psu, strata = ~str, weights = ~w, data = d)
  fd <- as_fsvy(d, ids = psu, strata = str, weights = w)
  tb <- fs_tab(fd, cat, by = reg)
  expect_named(tb, c("reg", "cat", "n", "N", "prop", "prop_se", "prop_low", "prop_upp"))
  ref <- as.data.frame(svytable(~reg + cat, sd))
  ref <- ref[order(ref$reg, ref$cat), ]
  expect_equal(tb$N, ref$Freq, tolerance = 1e-10)
  expect_equal(tb$n, as.vector(t(table(d$reg, d$cat))))
  expect_equal(sum(tb$prop[tb$reg == "A"]), 1)
  pc <- fs_tab(fd, cat, by = reg, percent = TRUE)
  expect_equal(pc$prop, 100 * tb$prop)
  # groups of the design are used when `by` is missing
  expect_equal(fs_tab(group_by(fd, reg), cat), tb)
})

test_that("fpc given as sampling fraction equals population sizes", {
  d <- sim_data(lonely = 0)
  npsu <- tapply(d$psu, d$str, function(v) length(unique(v)))
  d$frac <- npsu[as.character(d$str)] / d$npop
  a <- as_fsvy(d, ids = psu, strata = str, weights = w, fpc = npop) |> summarise(m = fs_mean(y))
  b <- as_fsvy(d, ids = psu, strata = str, weights = w, fpc = frac) |> summarise(m = fs_mean(y))
  expect_equal(a, b)
  ref <- svymean(~y, svydesign(ids = ~psu, strata = ~str, weights = ~w, fpc = ~frac, data = d))
  expect_equal(b$m_se, as.vector(SE(ref)), tolerance = 1e-8)
})

test_that("formula and tidyselect specifications are equivalent", {
  d <- sim_data(lonely = 0)
  a <- as_fsvy(d, ids = ~psu, strata = ~str, weights = ~w) |> summarise(m = fs_mean(y))
  b <- as_fsvy(d, ids = psu, strata = str, weights = w) |> summarise(m = fs_mean(y))
  expect_equal(a, b)
  # no clusters, no strata
  ref <- svymean(~y, svydesign(ids = ~1, weights = ~w, data = d))
  got <- as_fsvy(d, weights = w) |> summarise(m = fs_mean(y))
  expect_equal(got$m_se, as.vector(SE(ref)), tolerance = 1e-8)
})

test_that("mutate, select, rename and group_by(.add) keep the design", {
  d <- sim_data(lonely = 0)
  fd <- as_fsvy(d, ids = psu, strata = str, weights = w)
  base <- fd |> group_by(reg) |> summarise(m = fs_mean(y))
  got <- fd |>
    mutate(y2 = y / 2) |>
    select(reg, y2) |>
    rename(region = reg) |>
    group_by(region) |>
    summarise(m = fs_mean(y2))
  expect_equal(got$m_se, base$m_se / 2)
  g2 <- fd |> group_by(reg) |> group_by(sex, .add = TRUE)
  expect_equal(group_vars(g2), c("reg", "sex"))
  # .by works like group_by
  expect_equal(summarise(fd, m = fs_mean(y), .by = reg), base)
})

test_that("statistics outside summarise give a clear error", {
  expect_error(fs_mean(1:3), "inside")
})

test_that("factor input to fs_mean is rejected", {
  d <- sim_data(lonely = 0)
  fd <- as_fsvy(d, ids = psu, strata = str, weights = w)
  expect_error(summarise(fd, m = fs_mean(sex)), "fs_prop")
})
