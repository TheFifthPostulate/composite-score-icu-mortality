# tests/ebm_score.R ---------------------------------------------------------------
# THE EBM COMPARATOR, STEP 3: score the predictions of python/ebm_fit.py through
# the same metric code as every other arm (docs/ebm_comparator_20261028.md).
#
#   default   MIMIC-IV training out of fold, and eICU (frozen final EBM, apply only)
#   --test    the MIMIC-IV held-out test set: the EBM's SINGLE first look (author,
#             2026-10-20). Refuses if an `ebmtest` run already exists.
#
# Each site reports the EBM's AUROC and AUPRC with patient-bootstrap intervals
# (score_metrics()), and the paired difference against the summed evidence
# (llr_sum) and against the booster on the same 99 variables (xgb_feat), on the
# same stays, from the runs of record. The EBM score is put on the arms' scale,
# logit(p) - logit(p_bar_train), as xgb_oof_perfold() does; the rank metrics do
# not depend on it.
#
# Fits nothing. Aggregates only (hard rule 1).
#
#   Rscript tests/ebm_score.R
#   Rscript tests/ebm_score.R --test
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(arrow); library(yaml) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
TEST <- "--test" %in% args
# --additive: the main-effects-only EBM (python/ebm_fit.py --additive), a model
# of its own with its own single test look (author, 2026-10-28).
ADD  <- "--additive" %in% args
# --jointgam <variant>: the joint GAM comparator (tests/jointgam.R), which writes
# its predictions in this script's format and reads the same EBM design rows.
# Comparator arms carry no single-look restriction (author, 2026-10-30).
JG   <- { i <- match("--jointgam", args); if (is.na(i)) NA_character_ else args[i + 1L] }
if (!is.na(JG) && !JG %in% c("gated", "ungated")) stop("--jointgam takes gated or ungated", call. = FALSE)
PFX  <- if (!is.na(JG)) paste0("jgam", JG) else if (ADD) "ebmadd" else "ebm"
EBM  <- "out/cache/ebm"
FIT  <- if (!is.na(JG)) file.path("out/cache/jointgam", JG) else file.path(EBM, if (ADD) "fit_additive" else "fit")
ARM_NOTE <- if (!is.na(JG)) sprintf("joint GAM comparator (%s), fitted by tests/jointgam.R", JG) else
  sprintf("EBM comparator (%s), fitted by python/ebm_fit.py", if (ADD) "main effects only" else "default, with interactions")
# Reference runs whose stored scores fix the row order (options added 2026-11-13; defaults are the runs of
# record). The test-set run is --test-run because --test already selects the test look.
.opt <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
REF_OOF  <- file.path(.opt("--internal", "out/runs/internal_20260909T112643"), "tables/oof_scores.rds")
REF_EICU <- file.path(.opt("--external", "out/runs/external_20261025T211621"), "tables/scores_eicu.rds")
REF_TEST <- file.path(.opt("--test-run", "out/runs/test_20261025T211320"), "tables/scores_test.rds")

cfg    <- tar_read(cfg)
N_BOOT <- as.integer(cfg_req(cfg, "metrics", "n_boot"))
SEED   <- as.integer(cfg_req(cfg, "seed"))
rows   <- arrow::read_parquet(file.path(EBM, "rows.parquet"))
p_bar  <- mean(rows$y)                       # the training event rate, frozen
to_score <- function(p) { eps <- 1e-6; logit(pmin(pmax(p, eps), 1 - eps)) - logit(p_bar) }

#' One site: the EBM's metrics, and its paired difference against two arms.
site_rows <- function(site, ids, y, group, s_ebm, ref) {
  # The graph's `oof_scores` leaves the booster arms unnamed: they are in
  # train_ids order by construction (xgb_oof_perfold() scores in the order of
  # `y_train`), and the target asserts their length only. Name them by the LLR
  # arm's names after asserting those ARE the scored stays in order.
  if (is.null(names(ref$xgb_feat))) {
    if (!identical(names(ref$llr_sum), ids) || length(ref$xgb_feat) != length(ids)) {
      stop(site, ": unnamed booster scores cannot be aligned to the stays", call. = FALSE)
    }
    names(ref$xgb_feat) <- names(ref$llr_sum)
  }
  for (a in c("llr_sum", "xgb_feat")) {
    if (is.null(names(ref[[a]])) || !all(ids %in% names(ref[[a]]))) stop(site, ": reference ", a, " does not cover the stays", call. = FALSE)
  }
  m <- score_metrics(s_ebm, y, label = if (is.na(JG)) "ebm" else PFX, n_boot = N_BOOT, seed = SEED, group = group)
  base <- data.frame(site = site, arm = if (is.na(JG)) "ebm" else PFX, n = length(y), deaths = sum(y),
                     auroc = m$auroc, auroc_lo = m$auroc_lo, auroc_hi = m$auroc_hi,
                     auprc = m$auprc, auprc_lo = m$auprc_lo, auprc_hi = m$auprc_hi,
                     auprc_lift = m$auprc / m$event_rate, boot_unit = m$boot_unit, stringsAsFactors = FALSE)
  pd <- do.call(rbind, lapply(c("llr_sum", "xgb_feat"), function(a) {
    r <- unname(ref[[a]][ids])
    rbind(cbind(site = site, against = a, metric = "auroc",
                paired_boot_diff(s_ebm, r, y, metric = .auroc, n_boot = N_BOOT, seed = SEED, label = if (is.na(JG)) "ebm_minus" else paste0(PFX, "_minus"), group = group)),
          cbind(site = site, against = a, metric = "auprc",
                paired_boot_diff(s_ebm, r, y, metric = .auprc, n_boot = N_BOOT, seed = SEED, label = if (is.na(JG)) "ebm_minus" else paste0(PFX, "_minus"), group = group)))
  }))
  list(metrics = base, paired = pd)
}

out <- list()
if (!TEST) {
  # MIMIC-IV training, out of fold.
  tr_ids <- readRDS(file.path(EBM, "train_ids.rds"))
  po <- arrow::read_parquet(file.path(FIT, "pred_oof.parquet"))
  if (!identical(as.integer(po$row), as.integer(rows$row))) stop("out-of-fold predictions are not in training-row order", call. = FALSE)
  out$mimic_oof <- site_rows("mimic_oof", tr_ids, rows$y, rows$group, to_score(po$p), readRDS(REF_OOF))
  # eICU, the frozen final EBM.
  xcfg  <- yaml::read_yaml("config/external.yml")
  bundle <- load_bundle(cfg_req(xcfg, "bundle"), verbose = FALSE)
  cfg_e <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
  tabs_e <- load_tables(cfg_e$paths, cfg_e, site = "eicu", verbose = FALSE)
  e_ids <- readRDS(file.path(EBM, "eicu_ids.rds"))
  if (!identical(e_ids, as.character(tabs_e$cohort$stay_id))) stop("eICU predictions are not in cohort order", call. = FALSE)
  pe <- arrow::read_parquet(file.path(FIT, "pred_eicu.parquet"))
  if (nrow(pe) != length(e_ids)) stop("eICU prediction count differs from the cohort", call. = FALSE)
  out$eicu <- site_rows("eicu", e_ids, as.integer(tabs_e$cohort$mortality), patient_group_of(tabs_e$cohort, cfg_e, tabs_e$cohort$stay_id),
                        to_score(pe$p), readRDS(REF_EICU))
  run <- new_run(PFX, cfg, note = paste0(ARM_NOTE, ": MIMIC-IV out of fold and eICU; fits nothing here"))
} else {
  prior <- list.dirs("out/runs", recursive = FALSE)
  if (is.na(JG) && any(grepl(paste0("^", PFX, "test_"), basename(prior)))) stop("a ", PFX, "test run exists: this model's single test look is spent", call. = FALSE)
  tabs  <- tar_read(tabs)
  t_ids <- readRDS(file.path(EBM, "test_ids.rds"))
  pt <- arrow::read_parquet(file.path(FIT, "pred_test.parquet"))
  if (nrow(pt) != length(t_ids)) stop("test prediction count differs from the test stays", call. = FALSE)
  y_t <- as.integer(tabs$cohort$mortality[match(t_ids, as.character(tabs$cohort$stay_id))])
  if (anyNA(y_t)) stop("a test stay has no outcome", call. = FALSE)
  out$mimic_test <- site_rows("mimic_test", t_ids, y_t, patient_group_of(tabs$cohort, cfg, t_ids), to_score(pt$p), readRDS(REF_TEST))
  run <- new_run(paste0(PFX, "test"), cfg, note = paste0(ARM_NOTE, ": the MIMIC-IV held-out test set"))
}
met <- do.call(rbind, lapply(out, `[[`, "metrics")); pdf <- do.call(rbind, lapply(out, `[[`, "paired"))
save_table(run, met, "ebm_metrics"); save_table(run, pdf, "ebm_paired_differences")
if (!TEST) { fit_file <- if (is.na(JG)) "ebm_final.pkl" else "fits.rds"; file.copy(file.path(FIT, fit_file), file.path(run$path, fit_file)) }
finalize_run(run, extra = list(cache = EBM, ref_oof = REF_OOF, ref_eicu = REF_EICU, ref_test = if (TEST) REF_TEST else "",
                                fit_dir = FIT, settings = if (is.na(JG)) "interpret-core 0.7.8 defaults, random_state = seed (+ fold)" else "tests/jointgam.R, bam settings of config/config.yml"))
options(width = 200)
r3 <- function(d) { d[] <- lapply(d, function(v) if (is.numeric(v)) round(v, 4) else v); d }
print(r3(met), row.names = FALSE)
print(r3(pdf), row.names = FALSE)
cat("\nrun:", run$path, "\n")
