# Benchmark rapidsvy against survey and srvyr on the BPJS 2020 sample
# (01_kepesertaan). Usage: Rscript bench_bpjs.R <path to .rds or .dta> <out.rds>
suppressMessages({
  library(dplyr); library(survey); library(srvyr); library(rapidsvy); library(haven)
})
args <- commandArgs(TRUE)
path <- args[1]
out_file <- if (length(args) > 1) args[2] else "bench_results.rds"
d <- if (grepl("\\.rds$", path)) readRDS(path) else read_dta(path)
d <- d |>
  mutate(prov = as_factor(PSTV09), seg = as_factor(PSTV08), sex = as_factor(PSTV05),
         kab = as.integer(PSTV10), kelas3 = as.integer(PSTV07 == 3), pbpu = as.integer(PSTV08 == 4),
         age = as.numeric(as.Date("2020-12-31") - PSTV03) / 365.25,
         strat = interaction(PSTV14, seg, drop = TRUE)) |>
  select(PSTV02, PSTV14, PSTV15, prov, seg, sex, kab, kelas3, pbpu, age, strat) |>
  zap_labels()

options(survey.lonely.psu = "adjust")
tm <- function(expr) {
  t <- system.time(val <- force(expr))[["elapsed"]]
  list(time = t, value = val)
}
res <- list()

# survey reads `strata = ~a + b` as stage 1 and stage 2 strata, so cross them first
t_sd <- tm(svydesign(ids = ~PSTV02, strata = ~strat, weights = ~PSTV15, nest = TRUE, data = d))
sd <- t_sd$value
ss <- as_survey(sd)
t_fd <- tm(as_fsvy(d, ids = PSTV02, strata = c(PSTV14, seg), weights = PSTV15, nest = TRUE,
                   lonely_psu = "adjust"))
fd <- t_fd$value
res$design <- c(survey = t_sd$time, rapidsvy = t_fd$time)
cat("design:", res$design, "\n")

# 1. proportion of membership segment within province (34 x 5)
a <- tm(svyby(~seg, ~prov, sd, svymean))
b <- tm(ss |> group_by(prov, seg) |> summarise(p = survey_mean()))
f <- tm(fd |> group_by(prov, seg) |> summarise(p = fs_prop()))
res$prop_seg_prov <- c(survey = a$time, srvyr = b$time, rapidsvy = f$time,
                       max_abs_diff_se = max(abs(sort(unlist(SE(a$value))) - sort(f$value$p_se))))
cat("prop:", res$prop_seg_prov, "\n")

# 2. mean age by district (515 domains)
a <- tm(svyby(~age, ~kab, sd, svymean))
f <- tm(fd |> group_by(kab) |> summarise(m = fs_mean(age, vartype = c("se", "ci"))))
res$mean_age_kab <- c(survey = a$time, rapidsvy = f$time,
                      max_abs_diff_se = max(abs(unname(SE(a$value)) - f$value$m_se)))
cat("mean kab:", res$mean_age_kab, "\n")

# 3. median age by province
a <- tm(svyby(~age, ~prov, sd, svyquantile, quantiles = 0.5, ci = TRUE))
f <- tm(fd |> group_by(prov) |> summarise(m = fs_median(age)))
res$median_age_prov <- c(survey = a$time, rapidsvy = f$time,
                         max_abs_diff_se = max(abs(unname(SE(a$value)) - f$value$m_se)))
cat("median:", res$median_age_prov, "\n")

# 4. Rao-Scott chi-square: segment x sex
a <- tm(svychisq(~seg + sex, sd))
f <- tm(fs_chisq(fd, seg, sex))
res$chisq <- c(survey = a$time, rapidsvy = f$time,
               abs_diff_F = abs(unname(a$value$statistic) - f$value$statistic))
cat("chisq:", res$chisq, "\n")

# 5. logistic regression: PBPU membership
a <- tm(svyglm(pbpu ~ age + sex + prov, sd, family = quasibinomial()))
f <- tm(fs_glm(fd, pbpu ~ age + sex + prov, family = quasibinomial()))
res$glm <- c(survey = a$time, rapidsvy = f$time,
             max_rel_diff_se = max(abs(sqrt(diag(vcov(a$value))) / sqrt(diag(vcov(f$value))) - 1)))
cat("glm:", res$glm, "\n")

res$threads <- getOption("rapidsvy.threads", parallel::detectCores(logical = FALSE) - 1)
res$n <- nrow(d)
saveRDS(res, out_file)
print(res)
