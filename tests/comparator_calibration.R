# tests/comparator_calibration.R -------------------------------------------------
# Traditional calibration of the comparator models (author, 2026-11-01): every
# model that states a probability, out of fold in the MIMIC-IV training stays, on
# the MIMIC-IV test set and at eICU, with all parameters frozen at MIMIC-IV.
#
# The score is NOT in this table by design (author, 2026-11-01): its stated risk
# is the mortality of its reference stratum, whose transport is reported in the
# strata tables. Nothing here reads or computes a calibration of the summed
# weights of evidence.
#
# Models: the gated and ungated joint GAMs, the EBM with and without
# interactions, XGBoost on the constructed variables and on the raw window
# summaries, and the recalibrated APACHE II acute physiology score, APACHE II
# total and SOFA (on the stays where the severity scores can be computed).
#
# Non-fitting consumer. Per-stay predictions come from:
#   - the EBM and joint GAM caches (out/cache/ebm, out/cache/jointgam), read as
#     tests/ebm_strata_bins.R reads them;
#   - out of fold, the `_targets` store (oof_xgb_feat, oof_xgb_raw,
#     severity_oof);
#   - test set and eICU, apply_bundle() and apply_severity() on the frozen bundle
#     named in config/external.yml (application, not fitting; hard rule 9).
# Before any number is written, each source must reproduce its run of record:
# the AUROC of the EBM and joint GAM scoring runs, and the calibration slope of
# the boosters and severity scores in internal_20260909T112643,
# test_20260912T215546 and external_20260909T112824 (tolerance 1e-3). The
# test-set rows describe the frozen models post hoc; they change no design.
#
# Measures, per model and setting: observed and mean predicted mortality, their
# ratio, calibration-in-the-large (the intercept of a logistic regression with
# the predicted log-odds as offset; 0 is calibrated), the calibration slope (1
# is calibrated), the integrated calibration index (Austin and Steyerberg 2019:
# the mean absolute difference between the predicted probability and a smoothed
# observed probability) and the smoothed calibration curve on a fixed grid.
# Intervals: profile likelihood, stays as independent (one stay per patient).
# Aggregates only (hard rule 1).
#
#   Rscript tests/comparator_calibration.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(arrow); library(yaml); library(mgcv)
  library(xgboost); library(qs2) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

TOL  <- 1e-3
EPS  <- 1e-6
GRID <- stats::qlogis(seq(0.005, 0.80, length.out = 80))   # curve grid, predicted log-odds
# Input runs (options added 2026-11-13; each defaults to the run of record, a fresh chain passes its own):
#   --internal --test --external, and the comparator scoring runs --jgamgated --jgamgatedtest --jgamungated
#   --jgamungatedtest --ebmadd --ebmaddtest --ebm --ebmtest.
args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
RUNS <- list(internal = .opt("--internal", "out/runs/internal_20260909T112643"), test = .opt("--test", "out/runs/test_20260912T215546"),
             external = .opt("--external", "out/runs/external_20260909T112824"))
CACHED <- list(
  jgam_gated   = list(fit = "out/cache/jointgam/gated",   score = .opt("--jgamgated",   "out/runs/jgamgated_20261030T151327"),   test = .opt("--jgamgatedtest",   "out/runs/jgamgatedtest_20261030T151338")),
  jgam_ungated = list(fit = "out/cache/jointgam/ungated", score = .opt("--jgamungated", "out/runs/jgamungated_20261030T150052"), test = .opt("--jgamungatedtest", "out/runs/jgamungatedtest_20261030T150101")),
  ebm_additive = list(fit = "out/cache/ebm/fit_additive", score = .opt("--ebmadd",      "out/runs/ebmadd_20261028T234749"),      test = .opt("--ebmaddtest",      "out/runs/ebmaddtest_20261028T234803")),
  ebm          = list(fit = "out/cache/ebm/fit",          score = .opt("--ebm",         "out/runs/ebm_20261028T224038"),         test = .opt("--ebmtest",         "out/runs/ebmtest_20261028T224126")))
SETTINGS <- c("mimic_oof", "mimic_test", "eicu")

check <- function(what, got, want) {
  if (!is.finite(got) || !is.finite(want) || abs(got - want) > TOL)
    stop(sprintf("%s: reproduced %.5f, run of record %.5f", what, got, want), call. = FALSE)
}
clamp <- function(p) pmin(pmax(p, EPS), 1 - EPS)

calib <- function(p, y, model, setting) {
  ok <- is.finite(p) & !is.na(y)
  p <- clamp(p[ok]); y <- as.integer(y[ok]); lp <- stats::qlogis(p)
  f_citl  <- stats::glm(y ~ 1, offset = lp, family = stats::binomial(), na.action = na.fail)
  f_slope <- stats::glm(y ~ lp, family = stats::binomial(), na.action = na.fail)
  ci_c <- suppressMessages(stats::confint(f_citl))
  ci_s <- suppressMessages(stats::confint(f_slope))[2, ]
  sm <- mgcv::gam(y ~ s(lp, k = 10, bs = "ts"), family = stats::binomial(), method = "REML", na.action = na.fail)
  ici <- mean(abs(as.numeric(stats::predict(sm, type = "response")) - p))
  inr <- GRID >= stats::quantile(lp, 0.01) & GRID <= stats::quantile(lp, 0.99)
  pr  <- stats::predict(sm, newdata = data.frame(lp = GRID[inr]), se.fit = TRUE)
  list(table = data.frame(model = model, setting = setting, n = length(y), deaths = sum(y),
                          observed = mean(y), predicted = mean(p), oe_ratio = mean(y) / mean(p),
                          citl = unname(stats::coef(f_citl)[1]), citl_lo = ci_c[1], citl_hi = ci_c[2],
                          slope = unname(stats::coef(f_slope)[2]), slope_lo = ci_s[1], slope_hi = ci_s[2],
                          ici = ici, ci_method = "glm_profile_stay_iid", stringsAsFactors = FALSE),
       curve = data.frame(model = model, setting = setting, predicted = stats::plogis(GRID[inr]),
                          observed = stats::plogis(pr$fit), lo = stats::plogis(pr$fit - 1.96 * pr$se.fit),
                          hi = stats::plogis(pr$fit + 1.96 * pr$se.fit), stringsAsFactors = FALSE))
}
slope_of <- function(p, y) unname(stats::coef(stats::glm(as.integer(y) ~ stats::qlogis(clamp(p)),
                                                          family = stats::binomial(), na.action = na.fail))[2])
ref_slope <- function(run, name) read.csv(file.path(run, "tables", paste0("llr_calibration_", name, ".csv")))$slope

# --- frozen bundle, both sites ------------------------------------------------
rc        <- yaml::read_yaml("config/internal.yml")
cfg_local <- load_config(rc$config)
xcfg      <- yaml::read_yaml("config/external.yml")
bpath     <- cfg_req(xcfg, "bundle")
if (!identical(bpath, rc$test_look$bundle)) stop("test-look and external bundles differ", call. = FALSE)
bundle <- load_bundle(bpath, cfg = cfg_local, strict = TRUE)
cfg_m  <- bundle_cfg(bundle, cfg_local$paths$mimiciv)
tabs_m <- load_tables(cfg_m$paths, cfg_m, site = "mimic", verbose = FALSE)
folds  <- assign_folds(tabs_m$cohort, cfg_m)
test_ids <- folds$stay_id[folds$split == "test"]
y_test <- as.integer(tabs_m$cohort$mortality[match(as.character(test_ids), as.character(tabs_m$cohort$stay_id))])
cfg_e  <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
tabs_e <- load_tables(cfg_e$paths, cfg_e, site = "eicu", verbose = FALSE)
eicu_ids <- tabs_e$cohort$stay_id
y_eicu <- as.integer(tabs_e$cohort$mortality)
if (anyNA(y_test) || anyNA(y_eicu)) stop("an outcome is missing", call. = FALSE)

preds <- list()   # model -> setting -> list(p, y)
put <- function(model, setting, p, y) {
  if (length(p) != length(y)) stop(model, " ", setting, ": prediction and outcome lengths differ", call. = FALSE)
  preds[[model]][[setting]] <<- list(p = as.numeric(p), y = as.integer(y))
}

# --- EBM and joint GAMs: cached predictions --------------------------------------
rows  <- arrow::read_parquet("out/cache/ebm/rows.parquet")
t_ids <- readRDS("out/cache/ebm/test_ids.rds")
e_ids <- readRDS("out/cache/ebm/eicu_ids.rds")
if (!identical(as.character(t_ids), as.character(test_ids))) stop("cached test rows are not in test-set order", call. = FALSE)
if (!identical(as.character(e_ids), as.character(eicu_ids))) stop("cached eICU rows are not in cohort order", call. = FALSE)
for (m in names(CACHED)) {
  a  <- CACHED[[m]]
  po <- arrow::read_parquet(file.path(a$fit, "pred_oof.parquet"))
  if (!identical(as.integer(po$row), as.integer(rows$row))) stop(m, ": out-of-fold rows out of order", call. = FALSE)
  put(m, "mimic_oof",  po$p, rows$y)
  put(m, "mimic_test", arrow::read_parquet(file.path(a$fit, "pred_test.parquet"))$p, y_test)
  put(m, "eicu",       arrow::read_parquet(file.path(a$fit, "pred_eicu.parquet"))$p, y_eicu)
  ms <- read.csv(file.path(a$score, "tables", "ebm_metrics.csv"))
  mt <- read.csv(file.path(a$test,  "tables", "ebm_metrics.csv"))
  check(paste(m, "OOF AUROC"),  .auroc(po$p, rows$y), ms$auroc[ms$site == "mimic_oof"])
  check(paste(m, "test AUROC"), .auroc(preds[[m]]$mimic_test$p, y_test), mt$auroc[1])
  check(paste(m, "eICU AUROC"), .auroc(preds[[m]]$eicu$p, y_eicu), ms$auroc[ms$site == "eicu"])
}

# --- XGBoost: out of fold from the store, test and eICU from the bundle ---------
y_train <- tar_read(y_train)
for (d in c("xgb_feat", "xgb_raw")) {
  o <- tar_read_raw(paste0("oof_", d))
  put(d, "mimic_oof", o$p_hat, y_train)
  check(paste(d, "OOF slope"), slope_of(o$p_hat, y_train), ref_slope(RUNS$internal, paste0(d, "_oof")))
}
ap_t <- apply_bundle(bundle, tabs_m, cfg_m, test_ids, arms = c("xgb_feat", "xgb_raw"), verbose = FALSE)
ap_e <- apply_bundle(bundle, tabs_e, cfg_e, eicu_ids, arms = c("xgb_feat", "xgb_raw"), verbose = FALSE)
for (d in c("xgb_feat", "xgb_raw")) {
  pb <- bundle$xgb[[d]]$p_bar
  put(d, "mimic_test", stats::plogis(ap_t$scores[[d]] + stats::qlogis(pb)), y_test)
  put(d, "eicu",       stats::plogis(ap_e$scores[[d]] + stats::qlogis(pb)), y_eicu)
  check(paste(d, "test slope"), slope_of(preds[[d]]$mimic_test$p, y_test), ref_slope(RUNS$test, paste0(d, "_frozen")))
  check(paste(d, "eICU slope"), slope_of(preds[[d]]$eicu$p, y_eicu), ref_slope(RUNS$external, paste0(d, "_frozen")))
}

# --- severity scores, on the stays where they can be computed --------------------
sev_tr <- tar_read(severity_train)
sev_oof <- tar_read(severity_oof)
k <- sev_tr$keep
pb_oof <- mean(y_train[k])
sv_t <- apply_severity(bundle, tabs_m, cfg_m, test_ids, verbose = FALSE)
sv_e <- apply_severity(bundle, tabs_e, cfg_e, eicu_ids, verbose = FALSE)
pb_sev <- bundle$severity$p_bar
for (cell in SEVERITY_CELLS) {
  put(cell, "mimic_oof",  stats::plogis(sev_oof[[cell]] + stats::qlogis(pb_oof)), y_train[k])
  put(cell, "mimic_test", stats::plogis(sv_t$scores[[cell]][sv_t$keep] + stats::qlogis(pb_sev)), y_test[sv_t$keep])
  put(cell, "eicu",       stats::plogis(sv_e$scores[[cell]][sv_e$keep] + stats::qlogis(pb_sev)), y_eicu[sv_e$keep])
  check(paste(cell, "test slope"), slope_of(preds[[cell]]$mimic_test$p, preds[[cell]]$mimic_test$y), ref_slope(RUNS$test, paste0(cell, "_sev")))
  check(paste(cell, "eICU slope"), slope_of(preds[[cell]]$eicu$p, preds[[cell]]$eicu$y), ref_slope(RUNS$external, paste0(cell, "_sev")))
}

# --- measures --------------------------------------------------------------------
run <- new_run("compcal", cfg_local, note = "traditional calibration of the comparator models; the score is not included by design")
tab <- list(); cur <- list()
for (m in names(preds)) for (s in SETTINGS) {
  r <- calib(preds[[m]][[s]]$p, preds[[m]][[s]]$y, m, s)
  tab[[length(tab) + 1L]] <- r$table; cur[[length(cur) + 1L]] <- r$curve
}
tab <- do.call(rbind, tab); rownames(tab) <- NULL
save_table(run, tab, "comparator_calibration")
save_table(run, do.call(rbind, cur), "comparator_calibration_curves")
finalize_run(run, extra = list(bundle = bpath, models = paste(names(preds), collapse = ",")))
num <- vapply(tab, is.numeric, logical(1))
tab[num] <- lapply(tab[num], round, 3)
print(tab[, c("model", "setting", "n", "observed", "predicted", "oe_ratio", "citl", "slope", "slope_lo", "slope_hi", "ici")], row.names = FALSE)
cat("run:", run$path, "\n")
