# tests/hetero_oof_terms.R -----------------------------------------------------
# Heterogeneity stage T1 (docs/heterogeneity_plan_20261027.md, section 8b).
# OUT-OF-FOLD ANCHORED TERM CONTRIBUTIONS for the card-flag analysis.
#
# The patient card bolds a random variable when its anchored contribution
# (R/15_anchors.R: relative to a quiet day on the channel) is at least log 2.
# To count bold values for MIMIC-IV training stays honestly, each stay's terms
# must come from the fold model that held it out. The pipeline never saved fold
# models, so this script refits them: for each of the 19 channels, the card's
# model (`full`, or `meas` for an unpaired channel), in each of the 5 folds,
# with the fold priors from the targets store -- 95 fits, exactly the oof jobs
# of the pipeline (`fit_one()`, `job_ids()`, `priors_for()`). CLAUDE.md rule 7
# was relaxed by the author on 2026-10-28 for this purpose: the stripped fold
# models are saved under the run, and never enter a bundle or a score.
#
# CHECKS, per fold job, before anything is saved as a result:
#   1. the refit's held-out L reproduces the stored `l_oof` (same stays, max
#      absolute difference below TOL);
#   2. the anchored terms plus the reference evidence reproduce that L.
# At eICU the frozen final models of the bundle give the terms directly
# (application, no fitting), checked against the stage A weights of evidence.
# The MIMIC-IV test set is not touched.
#
# RESUMABLE. Each fold job is written to a compute cache (default
# out/cache/oofterms/) keyed by the job and a hash of its inputs; a rerun
# skips cached jobs. The run directory gets the assembled tables, the checks
# and the stripped models.
#
# AGGREGATE CONSOLE OUTPUT ONLY (rule 1): timings, counts, largest differences.
# The per-stay term tables are written to files and never printed.
#
#   Rscript tests/hetero_oof_terms.R --benchmark        # one fit, timed
#   Rscript tests/hetero_oof_terms.R [--skip-eicu] [--signals a,b]
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(mgcv); library(yaml); library(qs2) })
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
BENCH     <- "--benchmark" %in% args
SKIP_EICU <- "--skip-eicu" %in% args
STAYTAB   <- .opt("--staytab", "out/runs/staytab_20261027T215332")
CACHE     <- .opt("--cache", "out/cache/oofterms")
ONLY      <- .opt("--signals", NULL)
ANCHOR_N_PATH <- "paper/figs/coefs/anchor_n.csv"
TOL <- 1e-8

cfg    <- tar_read(cfg)
tabs   <- tar_read(tabs)
folds  <- tar_read(folds)
priors <- tar_read(priors)
l_oof  <- tar_read(l_oof)
jobs   <- layer1_jobs(cfg)
oj     <- jobs[jobs$fit & jobs$role == "oof", , drop = FALSE]
# The model whose L enters the card and the score (the alias rule of l_matrix()).
card_model <- function(sg) if (any(oj$signal == sg & oj$model == "full")) "full" else "meas"
SIGS <- cfg$signals
if (!is.null(ONLY)) { SIGS <- strsplit(ONLY, ",", fixed = TRUE)[[1]]; stopifnot(all(SIGS %in% cfg$signals)) }
FOLDS <- sort(unique(oj$fold))
ANCHOR_N <- read.csv(ANCHOR_N_PATH, stringsAsFactors = FALSE)
n_anc_of <- function(sg) { n <- ANCHOR_N$n_median[match(sg, ANCHOR_N$signal)]; if (is.na(n)) stop("no anchor n for ", sg); n }
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

#' Anchored contribution of every term of model `b` at the rows of `Xp`, the
#' same blocks as the card (run/patient_card.R): one column per smooth, one per
#' parametric column, each relative to its anchor; plus the reference evidence.
anchored_terms <- function(b, Xp, pri, n_anc) {
  cf <- stats::coef(b)
  tm <- list()
  for (s in b$smooth) {
    idx  <- s$first.para:s$last.para
    row0 <- as.numeric(.smooth_rows(s, term_anchor(s$term[1], pri, n_anc)))
    tm[[s$label]] <- as.numeric(sweep(Xp[, idx, drop = FALSE], 2, row0) %*% cf[idx])
  }
  in_sm <- unlist(lapply(b$smooth, function(s) s$first.para:s$last.para))
  for (j in setdiff(seq_along(cf), c(in_sm, match("(Intercept)", names(cf))))) tm[[names(cf)[j]]] <- Xp[, j] * cf[[j]]
  list(terms = as.data.frame(tm, check.names = FALSE), ref = reference_evidence(b, pri, n_anc)$value)
}

#' The intervention covariates the card reads to decide what it displays: the
#' exposure (fraction or hours with a dose) and the admission indicator.
display_cols <- function(sg, nd, cfg) {
  cols <- unlist(lapply(interventions_of(sg, cfg), function(iv)
    intersect(paste0(iv, "__", c("exposure_frac", "n_hours", "present_at_admission")), names(nd))))
  if (length(cols)) stats::setNames(as.data.frame(nd[, cols, drop = FALSE]), paste0("x:", cols)) else NULL
}
# cbind() of data frames fails on a NULL part (an unpaired channel has no display columns).
cbind_parts <- function(...) do.call(cbind, Filter(Negate(is.null), list(...)))

# --- one fold job ---------------------------------------------------------------
run_job <- function(sg, k) {
  md  <- card_model(sg)
  pri <- priors_for(priors, sg, "oof", k)
  ids <- job_ids("oof", k, folds)
  key <- substr(digest::digest(list(sg, md, k, build_formula(sg, md, cfg), pri, ids, n_anc_of(sg)), algo = "md5"), 1, 12)
  f_terms <- file.path(CACHE, sprintf("%s_f%d_%s.rds", sg, k, key))
  f_model <- file.path(CACHE, sprintf("%s_f%d_%s_model.qs2", sg, k, key))
  if (file.exists(f_terms) && file.exists(f_model)) return(readRDS(f_terms))

  t0 <- proc.time()[["elapsed"]]
  r  <- fit_one(sg, md, tabs, cfg, pri, fit_ids = ids$fit_ids, predict_ids = ids$predict_ids,
                keep_model = TRUE, role = "oof", fold = k)
  t_fit <- proc.time()[["elapsed"]] - t0
  b  <- r$model
  nd <- signal_frame(sg, md, tabs, cfg, pri, stay_ids = ids$predict_ids, stage = "predict")
  Xp <- stats::predict(b, newdata = nd, type = "lpmatrix", discrete = FALSE)
  at <- anchored_terms(b, Xp, pri, n_anc_of(sg))

  sid  <- as.character(nd$stay_id)
  L    <- r$l$l[match(sid, as.character(r$l$stay_id))]
  st   <- l_oof[l_oof$signal == sg & l_oof$model == md & l_oof$role == "oof" & l_oof$fold == k, ]
  if (nrow(st) != length(sid) || !setequal(as.character(st$stay_id), sid)) {
    stop(sprintf("%s fold %d: the refit predicts a different stay set from the stored l_oof", sg, k), call. = FALSE)
  }
  d_refit <- max(abs(L - st$l[match(sid, as.character(st$stay_id))]))
  d_terms <- max(abs(rowSums(at$terms) + at$ref - L))
  if (!is.finite(d_refit) || d_refit > TOL) stop(sprintf("%s fold %d: refit does not reproduce l_oof (max diff %.3g)", sg, k, d_refit), call. = FALSE)
  if (!is.finite(d_terms) || d_terms > TOL) stop(sprintf("%s fold %d: terms do not add up to L (max diff %.3g)", sg, k, d_terms), call. = FALSE)

  out <- list(signal = sg, model = md, fold = k,
              table = cbind_parts(data.frame(stay_id = nd$stay_id, fold = k, L = L, ref = at$ref), at$terms, display_cols(sg, nd, cfg)),
              check = data.frame(site = "mimic_oof", signal = sg, model = md, fold = k, n_stays = length(sid),
                                 max_diff_refit_vs_stored = d_refit, max_diff_terms_vs_L = d_terms,
                                 fit_seconds = round(t_fit, 1)))
  qs2::qs_save(strip_gam(b), f_model)
  saveRDS(out, f_terms)
  out
}

if (BENCH) {
  sg <- if (is.null(ONLY)) "heart_rate" else SIGS[1]
  t0 <- proc.time()[["elapsed"]]
  z  <- run_job(sg, FOLDS[1])
  cat(sprintf("benchmark %s fold %d (%s): fit %.1f s, whole job %.1f s; %d held-out stays; max diff refit vs stored %.2g, terms vs L %.2g\n",
              sg, FOLDS[1], z$model, z$check$fit_seconds, proc.time()[["elapsed"]] - t0, z$check$n_stays,
              z$check$max_diff_refit_vs_stored, z$check$max_diff_terms_vs_L))
  cat(sprintf("projected for %d fold jobs at this rate: %.0f min\n", 19L * length(FOLDS), 19 * length(FOLDS) * z$check$fit_seconds / 60))
  quit(save = "no")
}

# --- all fold jobs ----------------------------------------------------------------
run <- new_run("oofterms", cfg, note = "heterogeneity T1: fold refits and out-of-fold anchored terms; eICU terms from the frozen bundle")
dir.create(file.path(run$path, "models"))
checks <- list()
t_all <- proc.time()[["elapsed"]]
for (sg in SIGS) {
  parts <- lapply(FOLDS, function(k) {
    z <- run_job(sg, k)
    cat(sprintf("  %-18s fold %d  %5d stays  fit %6.1f s  refit diff %.1e  terms diff %.1e\n", sg, k,
                z$check$n_stays, z$check$fit_seconds, z$check$max_diff_refit_vs_stored, z$check$max_diff_terms_vs_L))
    z
  })
  checks[[sg]] <- do.call(rbind, lapply(parts, `[[`, "check"))
  tab <- do.call(rbind, lapply(parts, `[[`, "table"))
  save_object(run, tab, paste0("terms_mimic_oof_", sg))
  for (z in parts) file.copy(list.files(CACHE, pattern = sprintf("^%s_f%d_.*_model\\.qs2$", sg, z$fold), full.names = TRUE),
                             file.path(run$path, "models"), overwrite = FALSE)
}
cat(sprintf("fold jobs done in %.1f min\n", (proc.time()[["elapsed"]] - t_all) / 60))

# --- eICU: terms of the frozen final models ---------------------------------------
if (!SKIP_EICU) {
  xcfg   <- yaml::read_yaml("config/external.yml")
  bundle <- load_bundle(cfg_req(xcfg, "bundle"), verbose = FALSE)
  cfg_e  <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
  tabs_e <- load_tables(cfg_e$paths, cfg_e, site = "eicu", verbose = FALSE)
  tab_e  <- readRDS(file.path(STAYTAB, "tables", "stay_table_eicu.rds"))
  W_e    <- readRDS(file.path(STAYTAB, "tables", "woe_eicu.rds"))
  for (sg in SIGS) {
    md <- card_model(sg)
    b  <- bundle$models[[paste(sg, md, sep = "/")]]
    if (is.null(b)) stop("bundle has no ", sg, "/", md, call. = FALSE)
    pri <- priors_for(bundle$priors, sg, "final", NA_integer_)
    nd  <- signal_frame(sg, md, tabs_e, cfg_e, pri, stay_ids = tab_e$stay_id, stage = "predict")
    Xp  <- stats::predict(b, newdata = nd, type = "lpmatrix", discrete = FALSE)
    L   <- as.numeric(Xp %*% stats::coef(b)) - logit(pri$p_bar)
    at  <- anchored_terms(b, Xp, pri, n_anc_of(sg))
    ref_L <- W_e[match(as.character(nd$stay_id), as.character(tab_e$stay_id)), sg]
    d_stored <- max(abs(L - ref_L)); d_terms <- max(abs(rowSums(at$terms) + at$ref - L))
    if (!is.finite(d_stored) || d_stored > TOL) stop(sprintf("eICU %s: L differs from the stage A weights of evidence (%.3g)", sg, d_stored), call. = FALSE)
    if (!is.finite(d_terms) || d_terms > TOL) stop(sprintf("eICU %s: terms do not add up to L (%.3g)", sg, d_terms), call. = FALSE)
    save_object(run, cbind_parts(data.frame(stay_id = nd$stay_id, L = L, ref = at$ref), at$terms, display_cols(sg, nd, cfg_e)),
                paste0("terms_eicu_", sg))
    checks[[paste0("eicu_", sg)]] <- data.frame(site = "eicu", signal = sg, model = md, fold = NA_integer_, n_stays = nrow(nd),
                                                max_diff_refit_vs_stored = d_stored, max_diff_terms_vs_L = d_terms, fit_seconds = NA_real_)
    cat(sprintf("  eICU %-18s %6d stays  L diff %.1e  terms diff %.1e\n", sg, nrow(nd), d_stored, d_terms))
  }
}

ck <- do.call(rbind, checks)
save_table(run, ck, "oof_terms_checks", subdir = "diagnostics")
finalize_run(run, extra = list(staytab = STAYTAB, cache = CACHE, tol = TOL, n_fold_jobs = sum(ck$site == "mimic_oof"),
                                eicu = !SKIP_EICU))
cat(sprintf("\n%d checks, all within %.0e. run: %s\n", nrow(ck), TOL, run$path))
