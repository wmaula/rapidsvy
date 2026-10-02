skip_if_not_installed("survey")
library(survey)
library(dplyr)

tol <- 1e-8

for (lp in c("adjust", "average", "remove", "certainty")) {
  test_that(paste("mean, total, ratio, deff by domain match svyby, lonely =", lp), {
    d <- sim_data()
    old <- options(survey.lonely.psu = lp)
    on.exit(options(old))
    sd <- svydesign(ids = ~psu, strata = ~str, weights = ~w, data = d)
    fd <- as_fsvy(d, ids = psu, strata = str, weights = w, lonely_psu = lp)

    ref <- svyby(~y, ~reg, sd, svymean, deff = TRUE)
    got <- fd |> group_by(reg) |> summarise(m = fs_mean(y, deff = TRUE))
    expect_equal(got$m, unname(coef(ref)), tolerance = tol)
    expect_equal(got$m_se, unname(SE(ref)), tolerance = tol)
    expect_equal(got$m_deff, unname(deff(ref)), tolerance = 1e-6)

    ref <- svyby(~y, ~reg + sex, sd, svytotal, deff = TRUE)
    got <- fd |> group_by(reg, sex) |> summarise(t = fs_total(y, deff = TRUE)) |> arrange(sex, reg)
    expect_equal(got$t, unname(coef(ref)), tolerance = tol)
    expect_equal(got$t_se, unname(SE(ref)), tolerance = tol)
    expect_equal(got$t_deff, unname(deff(ref)), tolerance = 1e-6)

    ref <- svyby(~y, ~reg, denominator = ~x, sd, svyratio)
    got <- fd |> group_by(reg) |> summarise(r = fs_ratio(y, x))
    expect_equal(got$r, unname(coef(ref)), tolerance = tol)
    expect_equal(got$r_se, unname(SE(ref)), tolerance = tol)

    # whole sample
    ref <- svymean(~y, sd)
    got <- fd |> summarise(m = fs_mean(y))
    expect_equal(got$m_se, as.vector(SE(ref)), tolerance = tol)
  })
}

test_that("lonely_psu = 'fail' errors like survey", {
  d <- sim_data()
  fd <- as_fsvy(d, ids = psu, strata = str, weights = w)
  expect_error(summarise(fd, m = fs_mean(y)), "only one PSU")
})

test_that("proportions match svyby(svymean) on factors and svyciprop intervals", {
  d <- sim_data()
  old <- options(survey.lonely.psu = "adjust")
  on.exit(options(old))
  sd <- svydesign(ids = ~psu, strata = ~str, weights = ~w, data = d)
  fd <- as_fsvy(d, ids = psu, strata = str, weights = w, lonely_psu = "adjust")

  ref <- svyby(~cat, ~reg, sd, svymean, deff = TRUE)
  got <- fd |> group_by(reg, cat) |> summarise(p = fs_prop(deff = TRUE)) |> arrange(cat, reg)
  expect_equal(got$p, unname(as.vector(coef(ref))), tolerance = tol)
  expect_equal(got$p_se, unname(unlist(SE(ref))), tolerance = tol)
  expect_equal(got$p_deff, unname(unlist(deff(ref))), tolerance = 1e-6)

  df <- degf(sd)
  for (m in c("logit", "mean", "asin", "beta", "xlogit")) {
    sub <- subset(sd, reg == "B")
    ci <- attr(svyciprop(~I(cat == "y"), sub, method = m, df = df), "ci")
    got <- fd |> group_by(reg, cat) |>
      summarise(p = fs_prop(vartype = "ci", ci_method = m)) |>
      filter(reg == "B", cat == "y")
    if (m == "beta") {
      # survey's beta uses nrow() of the design passed in
      ci <- attr(svyciprop(~I(cat == "y"), sd, method = m, df = df), "ci")
      got <- fd |> mutate(yy = cat == "y") |> group_by(yy) |>
        summarise(p = fs_prop(vartype = "ci", ci_method = m)) |> filter(yy)
    }
    expect_equal(c(got$p_low, got$p_upp), unname(ci), tolerance = 1e-6, label = m)
  }
})

test_that("domains via filter() match subset()", {
  d <- sim_data()
  old <- options(survey.lonely.psu = "adjust")
  on.exit(options(old))
  sd <- svydesign(ids = ~psu, strata = ~str, weights = ~w, data = d)
  fd <- as_fsvy(d, ids = psu, strata = str, weights = w, lonely_psu = "adjust")
  ref <- svymean(~y, subset(sd, sex == "F" & y > 100))
  got <- fd |> filter(sex == "F", y > 100) |> summarise(m = fs_mean(y, vartype = c("se", "ci")))
  expect_equal(got$m_se, as.vector(SE(ref)), tolerance = tol)
  expect_equal(fs_degf(filter(fd, sex == "F", y > 100)), degf(subset(sd, sex == "F" & y > 100)))
  ci <- confint(ref, df = degf(subset(sd, sex == "F" & y > 100)))
  expect_equal(c(got$m_low, got$m_upp), as.vector(ci), tolerance = tol)
})

test_that("missing values: na.rm = TRUE matches survey, FALSE gives NA", {
  d <- sim_data()
  old <- options(survey.lonely.psu = "adjust")
  on.exit(options(old))
  sd <- svydesign(ids = ~psu, strata = ~str, weights = ~w, data = d)
  fd <- as_fsvy(d, ids = psu, strata = str, weights = w, lonely_psu = "adjust")
  ref <- svyby(~y_na, ~reg, sd, svymean, na.rm = TRUE)
  got <- fd |> group_by(reg) |> summarise(m = fs_mean(y_na, na.rm = TRUE))
  expect_equal(got$m_se, unname(SE(ref)), tolerance = tol)
  got2 <- fd |> group_by(reg) |> summarise(m = fs_mean(y_na))
  expect_true(all(is.na(got2$m)))
})

test_that("fpc and nest = TRUE match survey", {
  d <- sim_data(lonely = 0)
  d$psu_num <- as.integer(sub(".*_", "", d$psu))   # ids repeat across strata
  sd <- svydesign(ids = ~psu_num, strata = ~str, weights = ~w, fpc = ~npop, nest = TRUE, data = d)
  fd <- as_fsvy(d, ids = psu_num, strata = str, weights = w, fpc = npop, nest = TRUE)
  ref <- svyby(~y, ~reg, sd, svymean)
  got <- fd |> group_by(reg) |> summarise(m = fs_mean(y))
  expect_equal(got$m_se, unname(SE(ref)), tolerance = tol)
  expect_error(as_fsvy(d, ids = psu_num, strata = str, weights = w), "nest = TRUE")
})

test_that("quantiles match svyquantile (qrule math, Woodruff)", {
  d <- sim_data()
  old <- options(survey.lonely.psu = "adjust")
  on.exit(options(old))
  sd <- svydesign(ids = ~psu, strata = ~str, weights = ~w, data = d)
  fd <- as_fsvy(d, ids = psu, strata = str, weights = w, lonely_psu = "adjust")
  got <- fd |> group_by(reg) |>
    summarise(q = fs_quantile(y, c(0.25, 0.5, 0.9), vartype = c("se", "ci")))
  for (i in seq_len(nrow(got))) {
    ref <- svyquantile(~y, subset(sd, reg == got$reg[i]), c(0.25, 0.5, 0.9), ci = TRUE)$y
    expect_equal(c(got$q_q25[i], got$q_q50[i], got$q_q90[i]), unname(ref[, 1]), tolerance = tol)
    expect_equal(c(got$q_q25_low[i], got$q_q50_low[i], got$q_q90_low[i]), unname(ref[, 2]), tolerance = tol)
    expect_equal(c(got$q_q25_upp[i], got$q_q50_upp[i], got$q_q90_upp[i]), unname(ref[, 3]), tolerance = tol)
    expect_equal(c(got$q_q25_se[i], got$q_q50_se[i], got$q_q90_se[i]), unname(ref[, 4]), tolerance = tol)
  }
  med <- fd |> summarise(m = fs_median(y))
  ref <- svyquantile(~y, sd, 0.5)$y
  expect_equal(c(med$m, med$m_se), unname(ref[c(1, 4)]), tolerance = tol)
})

test_that("Rao-Scott test matches svychisq", {
  d <- sim_data()
  old <- options(survey.lonely.psu = "adjust")
  on.exit(options(old))
  sd <- svydesign(ids = ~psu, strata = ~str, weights = ~w, data = d)
  fd <- as_fsvy(d, ids = psu, strata = str, weights = w, lonely_psu = "adjust")
  ref <- svychisq(~reg + cat, sd)
  got <- fs_chisq(fd, reg, cat)
  expect_equal(got$statistic, unname(ref$statistic), tolerance = 1e-6)
  expect_equal(c(got$ndf, got$ddf), unname(ref$parameter), tolerance = 1e-6)
  expect_equal(got$p.value, unname(ref$p.value), tolerance = 1e-6)
  ref <- svychisq(~reg + cat, sd, statistic = "Chisq")
  got <- fs_chisq(fd, reg, cat, statistic = "Chisq")
  expect_equal(got$statistic, unname(ref$statistic), tolerance = 1e-6)
  expect_equal(got$p.value, unname(ref$p.value), tolerance = 1e-6)
})

test_that("fs_glm matches svyglm", {
  d <- sim_data()
  old <- options(survey.lonely.psu = "adjust")
  on.exit(options(old))
  sd <- svydesign(ids = ~psu, strata = ~str, weights = ~w, data = d)
  fd <- as_fsvy(d, ids = psu, strata = str, weights = w, lonely_psu = "adjust")

  ref <- svyglm(b ~ x + sex + reg, sd, family = quasibinomial())
  got <- fs_glm(fd, b ~ x + sex + reg, family = quasibinomial())
  expect_equal(coef(got), coef(ref), tolerance = 1e-6)
  expect_equal(unname(vcov(got)), unname(vcov(ref)), tolerance = 1e-6)
  expect_equal(got$df.residual, ref$df.residual)
  expect_equal(unname(tidy(got)$p.value), unname(summary(ref)$coefficients[, 4]), tolerance = 1e-6)

  ref <- svyglm(y ~ x + reg, sd)
  got <- fs_glm(fd, y ~ x + reg)
  expect_equal(unname(vcov(got)), unname(vcov(ref)), tolerance = 1e-6)

  ref <- svyglm(round(y) ~ x + sex, subset(sd, reg != "A"), family = quasipoisson())
  got <- fs_glm(filter(fd, reg != "A"), round(y) ~ x + sex, family = quasipoisson())
  expect_equal(coef(got), coef(ref), tolerance = 1e-6)
  expect_equal(unname(vcov(got)), unname(vcov(ref)), tolerance = 1e-6)

  # aliased column
  d2 <- d
  d2$x2 <- 2 * d2$x
  sd2 <- svydesign(ids = ~psu, strata = ~str, weights = ~w, data = d2)
  ref <- svyglm(y ~ x + x2 + sex, sd2)
  got <- fs_glm(as_fsvy(d2, ids = psu, strata = str, weights = w, lonely_psu = "adjust"), y ~ x + x2 + sex)
  expect_equal(coef(got), coef(ref), tolerance = 1e-6)
})

test_that("conversion from survey and srvyr objects", {
  skip_if_not_installed("srvyr")
  d <- sim_data()
  old <- options(survey.lonely.psu = "adjust")
  on.exit(options(old))
  sd <- svydesign(ids = ~psu, strata = ~str, weights = ~w, data = d)
  fd <- as_fsvy(sd)
  ref <- svyby(~y, ~reg, sd, svymean)
  got <- fd |> group_by(reg) |> summarise(m = fs_mean(y))
  expect_equal(got$m_se, unname(SE(ref)), tolerance = tol)
  ss <- srvyr::as_survey(sd)
  got2 <- as_fsvy(ss) |> group_by(reg) |> summarise(m = fs_mean(y))
  expect_equal(got2$m_se, got$m_se)
})

test_that("results do not depend on the number of threads", {
  d <- sim_data(n_str = 200, lonely = 0)
  fd <- as_fsvy(d, ids = psu, strata = str, weights = w)
  old <- options(rapidsvy.threads = 1)
  a <- fd |> group_by(reg, cat) |> summarise(p = fs_prop(), m = fs_mean(y))
  options(rapidsvy.threads = 7)
  b <- fd |> group_by(reg, cat) |> summarise(p = fs_prop(), m = fs_mean(y))
  options(old)
  expect_equal(a, b, tolerance = 1e-12)
})
