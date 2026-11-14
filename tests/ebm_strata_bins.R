# tests/ebm_strata_bins.R ----------------------------------------------------------
# The EBM arms' 20-bin risk tables, in the format the targets graph exports for
# every other arm, so tests/strata_reporting.R can read them (author, 2026-10-29).
#
# For each EBM arm (`ebm`: default with interactions; `ebmadd`: main effects
# only): the score is logit(p) - logit(p_bar_train), as for the boosters; the
# cut points are score_cutpoints() of the MIMIC-IV out-of-fold score, exactly as
# the graph's `cutpoints` target; and risk_bins(breaks = cut points) is applied
# out of fold, at the MIMIC-IV test set and at eICU, exactly as the apply runners
# do. The test-set table is a post hoc description of each arm's single test
# look, already spent (ebmtest_20261028T224126, ebmaddtest_20261028T234803), as
# for the other arms in tests/strata_reporting.R.
#
# Writes per-bin counts and each arm's raw AUROC/AUPRC (copied from its scoring
# runs) into one run. Aggregates only (hard rule 1).
#
#   Rscript tests/ebm_strata_bins.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(arrow); library(yaml) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

EBM  <- "out/cache/ebm"
# The EBM scoring runs: --ebm --ebmtest --ebmadd --ebmaddtest (options added 2026-11-13; defaults are the
# runs of record).
args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
ARMS <- list(ebm    = list(fit = file.path(EBM, "fit"),          score = .opt("--ebm",    "out/runs/ebm_20261028T224038"),    test = .opt("--ebmtest",    "out/runs/ebmtest_20261028T224126")),
             ebmadd = list(fit = file.path(EBM, "fit_additive"), score = .opt("--ebmadd", "out/runs/ebmadd_20261028T234749"), test = .opt("--ebmaddtest", "out/runs/ebmaddtest_20261028T234803")))
# The joint GAM arms (tests/jointgam.R), named explicitly by their scoring runs,
# never picked as "the latest": --jointgam-gated <score run>,<test run>, and the
# same for --jointgam-ungated. They read the same EBM design rows.
for (v in c("gated", "ungated")) {
  i <- match(paste0("--jointgam-", v), args)
  if (is.na(i)) next
  r <- strsplit(args[i + 1L], ",", fixed = TRUE)[[1]]
  if (length(r) != 2L || !all(dir.exists(r))) stop("--jointgam-", v, " takes two existing run directories: <score>,<test>", call. = FALSE)
  ARMS[[paste0("jgam", v)]] <- list(fit = file.path("out/cache/jointgam", v), score = r[1], test = r[2])
}
cfg    <- tar_read(cfg)
N_BINS <- as.integer(cfg_req(cfg, "metrics", "n_bins"))
rows   <- arrow::read_parquet(file.path(EBM, "rows.parquet"))
p_bar  <- mean(rows$y)
to_score <- function(p) { eps <- 1e-6; logit(pmin(pmax(p, eps), 1 - eps)) - logit(p_bar) }

tabs  <- tar_read(tabs)
t_ids <- readRDS(file.path(EBM, "test_ids.rds"))
y_t   <- as.integer(tabs$cohort$mortality[match(t_ids, as.character(tabs$cohort$stay_id))])
xcfg  <- yaml::read_yaml("config/external.yml")
bundle <- load_bundle(cfg_req(xcfg, "bundle"), verbose = FALSE)
cfg_e <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
tabs_e <- load_tables(cfg_e$paths, cfg_e, site = "eicu", verbose = FALSE)
e_ids <- readRDS(file.path(EBM, "eicu_ids.rds"))
if (!identical(e_ids, as.character(tabs_e$cohort$stay_id))) stop("eICU predictions are not in cohort order", call. = FALSE)
y_e <- as.integer(tabs_e$cohort$mortality)
if (anyNA(y_t) || anyNA(y_e)) stop("an outcome is missing", call. = FALSE)

run <- new_run("ebmbins", cfg, note = "EBM arms' 20-bin risk tables on frozen out-of-fold cut points, for tests/strata_reporting.R")
summ <- list()
for (a in names(ARMS)) {
  fd <- ARMS[[a]]$fit
  po <- arrow::read_parquet(file.path(fd, "pred_oof.parquet"))
  if (!identical(as.integer(po$row), as.integer(rows$row))) stop(a, ": out-of-fold predictions are not in row order", call. = FALSE)
  s_oof <- to_score(po$p)
  s_t   <- to_score(arrow::read_parquet(file.path(fd, "pred_test.parquet"))$p)
  s_e   <- to_score(arrow::read_parquet(file.path(fd, "pred_eicu.parquet"))$p)
  if (length(s_t) != length(y_t) || length(s_e) != length(y_e)) stop(a, ": prediction counts differ from the cohorts", call. = FALSE)
  cp <- score_cutpoints(s_oof, n_bins = N_BINS)
  save_table(run, risk_bins(s_oof, rows$y, n_bins = N_BINS, breaks = cp), sprintf("risk_bins_%s_oof", a))
  save_table(run, risk_bins(s_t, y_t, n_bins = N_BINS, breaks = cp), sprintf("risk_bins_%s_test_frozen", a))
  save_table(run, risk_bins(s_e, y_e, n_bins = N_BINS, breaks = cp), sprintf("risk_bins_%s_eicu_frozen", a))
  m  <- read.csv(file.path(ARMS[[a]]$score, "tables", "ebm_metrics.csv"), stringsAsFactors = FALSE)
  mt <- read.csv(file.path(ARMS[[a]]$test,  "tables", "ebm_metrics.csv"), stringsAsFactors = FALSE)
  summ[[a]] <- data.frame(arm = a, site = c("mimic_train_oof", "mimic_test", "eicu"),
                          auroc = c(m$auroc[m$site == "mimic_oof"], mt$auroc, m$auroc[m$site == "eicu"]),
                          auprc = c(m$auprc[m$site == "mimic_oof"], mt$auprc, m$auprc[m$site == "eicu"]))
}
save_table(run, do.call(rbind, summ), "ebm_raw_summary")
finalize_run(run, extra = list(arms = paste(names(ARMS), collapse = ","), n_bins = N_BINS))
cat("run:", run$path, "\n")
