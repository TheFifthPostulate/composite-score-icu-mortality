# tests/bold_values_replicates.R ------------------------------------------------
# DO THE CARD'S BOLD VALUES COME BACK AFTER RETRAINING?
# (docs/bold_reproducibility_plan_20261028.md, "bold values"; author approved
# 2026-10-28.)
#
# The replicate stores hold one weight of evidence per channel, not its terms,
# so the bold random variables need the card's models refitted in every bag:
#
#   mimic   for bag b, channel g, fold k: the card's model of g (`full`, or
#           `meas` for an unpaired channel) refitted on the fold-k fitting
#           stays inside the bag, with fold k's out-of-fold priors, predicting
#           fold k's held-out stays: EXACTLY the `bootstrap` route of
#           tests/attr_replicates.R (`boot_arm()`: `.bam_fit()`, then
#           `.accept_outcome_fit()`, then predict).
#   eicu    for bag b, channel g: `fit_one()` on all training stays in the bag
#           with the bundle's frozen final priors, applied to the eICU cohort:
#           EXACTLY the LLR route of tests/attr_external_bags.R.
#
# CHECKS before a job's terms are used (CLAUDE.md rule 7, as relaxed 2026-10-28):
#   1. every bag's membership hash equals the internal store's;
#   2. the refit's L equals the stored `llr_full` bootstrap replicate of that
#      bag on every predicted stay (max absolute difference <= TOL);
#   3. the anchored terms plus the reference evidence reproduce that L.
# Tombstoned bags (excluded whole by the generators) are skipped.
#
# Each job's bold items are counted with `channel_items()` (tests/hetero_common.R,
# the code the card-flag analysis uses) and cached as one bitmask per stay
# (bit j = item j of the channel is bold). Jobs are cached per (site, bag,
# channel, fold) under --cache, so the script is resumable and the bags can be
# split across processes (--bags). `--assemble` reads the cache, verifies the
# unperturbed card (bag 0: the pipeline's own out-of-fold terms at MIMIC-IV, the
# frozen bundle's at eICU, from stage T1) against the card-flag analysis's
# per-stay counts, and writes the agreement tables over the pairs that
# `attr_contrast_pairs()` defines.
#
# Fits only reproducibility replicates (CLAUDE.md rule 9); nothing enters a
# bundle or a score. AGGREGATE CONSOLE OUTPUT ONLY (hard rule 1).
#
#   Rscript tests/bold_values_replicates.R --site mimic --benchmark
#   Rscript tests/bold_values_replicates.R --site mimic --bags 1-19
#   Rscript tests/bold_values_replicates.R --site eicu  --bags 1-40
#   Rscript tests/bold_values_replicates.R --site mimic --assemble
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(mgcv); library(yaml); library(qs2) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
source("tests/hetero_common.R")

args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
SITE     <- .opt("--site", "mimic")
if (!SITE %in% c("mimic", "eicu")) stop("--site must be mimic or eicu", call. = FALSE)
BENCH    <- "--benchmark" %in% args
ASSEMBLE <- "--assemble" %in% args
# The job cache holds finished bold masks, so it is keyed by everything that decides them (Codex review,
# finding 7, 2026-10-30): the counting code (tests/hetero_common.R) and the threshold. A cache written under
# other counting code is never reused; the default cache of the saved store (out/cache/boldvals/<site>,
# counted before the singular-tail fix) is therefore not read by a rerun.
CACHE_KEY <- substr(unname(tools::md5sum("tests/hetero_common.R")), 1, 10)
CACHE    <- file.path(.opt("--cache", "out/cache/boldvals"), paste0(SITE, "_", CACHE_KEY))
T1_RUN   <- .opt("--terms", "out/runs/oofterms_20261028T134954")
# Bag 0 is verified against the card-flag recount made with the same counting code: the 2026-10-30 runs, after
# the singular-tail fix. The saved store (boldvals_20261028T230049, boldvalsext_20261028T231235) was counted
# under the earlier rule and verified against heteroflags_20261028T153507 / _153910, as its manifest records.
# --heteroflags <run>: the card-flag run of this site (option added 2026-11-13; defaults are the runs above).
FLAGS    <- file.path(.opt("--heteroflags", c(mimic = "out/runs/heteroflags_20261030T204724",
                                              eicu  = "out/runs/heteroflags_20261030T204938")[[SITE]]),
                      "tables", c(mimic = "card_flags_stay_mimic_oof.rds", eicu = "card_flags_stay_eicu.rds")[[SITE]])
TOL      <- 1e-8
STRONG   <- log(2)                                   # paper/make_card_figure.R
ANCHOR_N <- read.csv("paper/figs/coefs/anchor_n.csv", stringsAsFactors = FALSE)
n_anc_of <- function(sg) { n <- ANCHOR_N$n_median[match(sg, ANCHOR_N$signal)]; if (is.na(n)) stop("no anchor n for ", sg); n }
ecfg     <- yaml::read_yaml("config/attribution_eval.yml")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

# --- the training side and the store ---------------------------------------------
cfg_t  <- tar_read(cfg)
folds  <- tar_read(folds)
tr     <- tar_read(train_ids)
priors <- tar_read(priors)
tabs_m <- tar_read(tabs)
STORE  <- .opt("--store", latest_run(if (SITE == "mimic") "attrgen" else "attrextgen", require_complete = TRUE))
DG     <- qs2::qs_read(file.path(STORE, "design.qs2"))
MAN    <- utils::read.csv(file.path(STORE, "manifest_replicates.csv"), stringsAsFactors = FALSE)
SUB    <- cfg_req(ecfg, "storage", "subdir")
ids_out <- as.character(DG$stay_id)

if (SITE == "mimic") {
  INTERNAL <- STORE
  cfg <- cfg_t
  .fold_vec <- folds$fold[match(tr, folds$stay_id)]
  .want <- attr_design_key(cfg, .fold_vec, DG$k_ti)
  if (nrow(attr_design_diff(DG$design_key, .want))) stop("the store's design key differs from the live design", call. = FALSE)
  if (nrow(attr_fingerprint_diff(DG$fingerprint, attr_design_fingerprint(cfg, .fold_vec, DG$k_ti, priors, ecfg, design_key = .want)))) {
    stop("the store's fingerprint differs from the live design", call. = FALSE)
  }
  if (!identical(ids_out, as.character(tr))) stop("the store's rows are not the training stays", call. = FALSE)
  sf <- tabs_m$signal_features
} else {
  source("tests/attr_external_common.R")
  INTERNAL <- file.path(dirname(STORE), DG$internal_store)
  bundle <- load_bundle(DG$bundle, verbose = FALSE)
  cfg    <- bundle_cfg(bundle, paths = cfg_req(yaml::read_yaml("config/external.yml"), "paths"))
  if (!identical(attr_key_hash(priors_final(priors)), attr_key_hash(bundle$priors))) stop("the targets store's final priors differ from the bundle's", call. = FALSE)
  .id <- eicu_store_identity(qs2::qs_read(file.path(INTERNAL, "design.qs2"))$design_key, bundle, cfg, ecfg, DG$stay_id, DG$methods)
  if (!identical(.id$design_key, DG$design_key) || nrow(attr_fingerprint_diff(DG$fingerprint, .id$fingerprint))) {
    stop("the eICU store's design key or fingerprint differs from the live design", call. = FALSE)
  }
  tabs_e <- load_tables(cfg$paths, cfg, site = "eicu", verbose = FALSE)
  if (!identical(as.character(tabs_e$cohort$stay_id), ids_out)) stop("the eICU cohort is not the store's row set", call. = FALSE)
  sf <- tabs_e$signal_features
}
SIGS <- as.character(unlist(cfg$signals))
card_model <- function(sg) if (is.na(spec_source(sg, "full", cfg))) "full" else "meas"

# --- the bags: verbatim from tests/attr_replicates.R, identity checked -------------
BAG_SEED0 <- as.integer(cfg_req(ecfg, "bootstrap", "seed_base"))
B_BOOT    <- as.integer(cfg_req(ecfg, "bootstrap", "b"))
bag_of <- function(seed) {
  gcol <- resample_cols(tabs_m$cohort, cfg_t)$fold_group
  src  <- if (gcol %in% names(folds)) folds else tabs_m$cohort
  gid  <- as.character(src[[gcol]])[match(tr, src$stay_id)]
  if (anyNA(gid)) stop("bag_of: `", gcol, "` is missing for some training stay", call. = FALSE)
  ug   <- unique(gid)
  with_seed(seed, { keep <- sample(ug, length(ug), replace = TRUE); gid %in% unique(keep) })
}
.ibt <- readRDS(file.path(INTERNAL, "diagnostics", "bootstrap_bags.rds"))
.tp  <- file.path(STORE, "llr_bootstrap_excluded.csv")
TOMB <- if (file.exists(.tp)) {
  .t <- utils::read.csv(.tp, stringsAsFactors = FALSE)
  sort(unique(as.integer(.t$boot_id[.t$design %in% c(attr_key_hash(DG$design_key), DG$migrated_from)])))
} else integer(0)

store_rep <- function(b) {
  z <- MAN[MAN$method == "llr_full" & MAN$route == "bootstrap" & MAN$boot_id == b & MAN$seed_id == 0L & MAN$draw_id == 0L, ]
  if (nrow(z) != 1L) stop(nrow(z), " stored llr_full bootstrap replicates for bag ", b, "; expected one", call. = FALSE)
  M <- qs2::qs_read(file.path(STORE, SUB, paste0(z$key, ".qs2")))
  attr_check_replicate(M, ids_out, SIGS, z$key)
  M
}

#' Anchored terms, as tests/hetero_oof_terms.R and the card compute them.
anchored_terms <- function(b, Xp, pri, n_anc) {
  cf <- stats::coef(b); tm <- list()
  for (s in b$smooth) {
    idx  <- s$first.para:s$last.para
    row0 <- as.numeric(.smooth_rows(s, term_anchor(s$term[1], pri, n_anc)))
    tm[[s$label]] <- as.numeric(sweep(Xp[, idx, drop = FALSE], 2, row0) %*% cf[idx])
  }
  in_sm <- unlist(lapply(b$smooth, function(s) s$first.para:s$last.para))
  for (j in setdiff(seq_along(cf), c(in_sm, match("(Intercept)", names(cf))))) tm[[names(cf)[j]]] <- Xp[, j] * cf[[j]]
  list(terms = as.data.frame(tm, check.names = FALSE), ref = reference_evidence(b, pri, n_anc)$value)
}
display_cols <- function(sg, nd) {
  cols <- unlist(lapply(interventions_of(sg, cfg), function(iv)
    intersect(paste0(iv, "__", c("exposure_frac", "n_hours", "present_at_admission")), names(nd))))
  if (length(cols)) stats::setNames(as.data.frame(nd[, cols, drop = FALSE]), paste0("x:", cols)) else NULL
}

#' Bold-item bitmask per stay from a term table, with the card's item rules.
bold_mask <- function(tt, sg, stay_id) {
  z  <- sf[sf$signal == sg, c("stay_id", "n_obs", "dx_trend_defined", "q05", "q95", "value_min", "value_max")]
  zz <- z[match(as.character(stay_id), as.character(z$stay_id)), ]
  ci <- channel_items(tt, sg, zz$n_obs, as.logical(zz$dx_trend_defined), cfg, strong = STRONG,
                      tails_eq = tails_equal(zz, sg, cfg))
  it <- attr(ci, "items")
  if (ncol(it) > 30L) stop("more than 30 items on channel ", sg, call. = FALSE)
  B  <- !is.na(it) & it >= STRONG
  mask <- as.integer(B %*% (2^(seq_len(ncol(B)) - 1L)))
  if (!identical(as.integer(rowSums(B)), as.integer(ci[, "n_bold"]))) stop("bitmask and bold count disagree on ", sg, call. = FALSE)
  list(mask = mask, n_items = as.integer(ci[, "n_items"]), item_names = colnames(it))
}

# --- one job ------------------------------------------------------------------------
ND_E <- list()
run_job <- function(b, sg, k, bag, M_store) {
  f_out <- file.path(CACHE, sprintf("b%02d_%s_f%s.rds", b, sg, if (is.na(k)) "all" else k))
  if (file.exists(f_out)) return(readRDS(f_out))
  md <- card_model(sg); fo <- build_formula(sg, md, cfg)
  t0 <- proc.time()[["elapsed"]]
  if (SITE == "mimic") {
    pri <- priors_for(priors, sg, "oof", fold = k)
    ji  <- job_ids("oof", k, folds)
    fit_ids <- intersect(ji$fit_ids, tr[bag])
    d_fit <- signal_frame(sg, md, tabs_m, cfg, pri, stay_ids = fit_ids)
    bm <- .bam_fit(fo, d_fit, cfg)
    .accept_outcome_fit(bm, list(signal = sg, model = md, role = "bag", fold = k, converged = .bam_converged(bm)), cfg)
    nd <- signal_frame(sg, md, tabs_m, cfg, pri, stay_ids = ji$predict_ids, stage = "predict")
  } else {
    pri <- priors_for(bundle$priors, sg, "final", NA_integer_)
    bm  <- fit_one(sg, md, tabs_m, cfg, pri, fit_ids = tr[bag], keep_model = TRUE, role = "final")$model
    if (is.null(ND_E[[sg]])) ND_E[[sg]] <<- signal_frame(sg, md, tabs_e, cfg, pri, stay_ids = DG$stay_id, stage = "predict")
    nd <- ND_E[[sg]]
  }
  t_fit <- proc.time()[["elapsed"]] - t0
  Xp   <- stats::predict(bm, newdata = nd, type = "lpmatrix", discrete = FALSE)
  L    <- as.numeric(Xp %*% stats::coef(bm)) - logit(pri$p_bar)
  rows <- match(as.character(nd$stay_id), ids_out)
  if (anyNA(rows)) stop("job stays outside the store's rows", call. = FALSE)
  d_store <- max(abs(L - M_store[rows, sg]))
  at <- anchored_terms(bm, Xp, pri, n_anc_of(sg))
  d_terms <- max(abs(rowSums(at$terms) + at$ref - L))
  if (!is.finite(d_store) || d_store > TOL) stop(sprintf("bag %d %s fold %s: refit does not reproduce the stored replicate (%.3g)", b, sg, k, d_store), call. = FALSE)
  if (!is.finite(d_terms) || d_terms > TOL) stop(sprintf("bag %d %s fold %s: terms do not add up to L (%.3g)", b, sg, k, d_terms), call. = FALSE)
  tt <- do.call(cbind, Filter(Negate(is.null), list(at$terms, display_cols(sg, nd))))
  bmk <- bold_mask(tt, sg, nd$stay_id)
  out <- list(rows = rows, mask = bmk$mask, n_items = bmk$n_items, item_names = bmk$item_names,
              check = data.frame(site = SITE, boot_id = b, signal = sg, fold = k, n_stays = length(rows),
                                 max_diff_vs_store = d_store, max_diff_terms_vs_L = d_terms, fit_seconds = round(t_fit, 1)))
  saveRDS(out, f_out)
  out
}
FOLDS <- if (SITE == "mimic") sort(unique(folds$fold[match(tr, folds$stay_id)])) else NA_integer_

# --- generate -----------------------------------------------------------------------
if (!ASSEMBLE) {
  spec <- .opt("--bags", paste0("1-", B_BOOT))
  BAGS <- sort(unique(unlist(lapply(strsplit(spec, ",")[[1]], function(r) { p <- as.integer(strsplit(r, "-")[[1]]); seq(p[1], p[length(p)]) }))))
  BAGS <- setdiff(BAGS, TOMB)
  if (BENCH) BAGS <- BAGS[1]
  cat(sprintf("%s: store %s; bags %s (tombstoned and skipped: %s)\n", SITE, basename(STORE),
              paste(range(BAGS), collapse = "-"), paste(TOMB, collapse = ",")))
  for (b in BAGS) {
    t0 <- proc.time()[["elapsed"]]
    bag <- bag_of(BAG_SEED0 + b)
    if (!identical(attr_key_hash(which(bag)), .ibt$membership_hash[.ibt$boot_id == b])) stop("bag ", b, " differs from the internal store's", call. = FALSE)
    M_store <- store_rep(b)
    sigs_b <- if (BENCH) SIGS[2] else SIGS
    for (sg in sigs_b) for (k in FOLDS) {
      z <- run_job(b, sg, k, bag, M_store)
      if (BENCH) cat(sprintf("benchmark bag %d %s fold %s: fit %.1f s, %d stays, diff vs store %.1e, terms %.1e\n",
                             b, sg, k, z$check$fit_seconds, z$check$n_stays, z$check$max_diff_vs_store, z$check$max_diff_terms_vs_L))
    }
    cat(sprintf("  bag %2d done in %.1f min\n", b, (proc.time()[["elapsed"]] - t0) / 60))
    rm(M_store); invisible(gc(verbose = FALSE))
  }
  quit(save = "no")
}

# --- assemble -----------------------------------------------------------------------
# Per bag: an n x 19 integer matrix of bold-item bitmasks. Bag 0 is the card itself.
POP10 <- vapply(0:1023, function(v) sum(bitwAnd(v, 2L^(0:9)) > 0), integer(1))
popcount <- function(v) POP10[v %% 1024L + 1L] + POP10[(v %/% 1024L) %% 1024L + 1L] + POP10[v %/% 1048576L + 1L]
stopifnot(popcount(c(0, 1, 1023, 1024, 2^19 - 1)) == c(0, 1, 10, 1, 19))
meas <- qs2::qs_read(file.path(STORE, "measured.qs2"))
BAGS <- setdiff(seq_len(B_BOOT), TOMB)
assemble_bag <- function(b) {
  K <- matrix(0L, length(ids_out), length(SIGS), dimnames = list(NULL, SIGS)); NI <- K; seen <- K
  ck <- list()
  for (sg in SIGS) for (k in FOLDS) {
    f <- file.path(CACHE, sprintf("b%02d_%s_f%s.rds", b, sg, if (is.na(k)) "all" else k))
    if (!file.exists(f)) stop("missing job: bag ", b, " ", sg, " fold ", k, call. = FALSE)
    z <- readRDS(f); K[z$rows, sg] <- z$mask; NI[z$rows, sg] <- z$n_items; seen[z$rows, sg] <- 1L
    ck[[length(ck) + 1L]] <- z$check
  }
  if (any((seen == 1L) != unname(meas[, SIGS]))) stop("bag ", b, ": the jobs do not cover exactly the measured cells", call. = FALSE)
  list(mask = K, n_items = NI, check = do.call(rbind, ck))
}
# Bag 0: the card's own terms, from stage T1, counted with the same function.
card_bag <- function() {
  K <- matrix(0L, length(ids_out), length(SIGS), dimnames = list(NULL, SIGS))
  for (sg in SIGS) {
    tt <- readRDS(file.path(T1_RUN, "tables", sprintf("terms_%s_%s.rds", if (SITE == "mimic") "mimic_oof" else "eicu", sg)))
    keep <- setdiff(names(tt), c("stay_id", "fold", "L", "ref"))
    bm <- bold_mask(tt[, keep, drop = FALSE], sg, tt$stay_id)
    K[match(as.character(tt$stay_id), ids_out), sg] <- bm$mask
  }
  K
}
K0 <- card_bag()
fl <- readRDS(FLAGS)
k0 <- rowSums(matrix(popcount(K0), nrow(K0)))[match(as.character(fl$stay_id), ids_out)]
if (!identical(as.integer(k0), as.integer(fl$k_val))) stop("bag 0 bold-value counts differ from the card-flag analysis", call. = FALSE)
cat(sprintf("bag 0 (the card): bold-value counts equal the card-flag analysis on all %d stays\n", nrow(fl)))

ENC <- list("0" = K0); CK <- list()
for (b in BAGS) { a <- assemble_bag(b); ENC[[as.character(b)]] <- a$mask; CK[[length(CK) + 1L]] <- a$check }
CK <- do.call(rbind, CK)
cat(sprintf("%d bags assembled; %d jobs; max diff vs store %.1e; max terms diff %.1e\n",
            length(BAGS), nrow(CK), max(CK$max_diff_vs_store), max(CK$max_diff_terms_vs_L)))

pair_metrics <- function(A, B) {
  inter <- rowSums(matrix(popcount(bitwAnd(A, B)), nrow(A))); uni <- rowSums(matrix(popcount(bitwOr(A, B)), nrow(A)))
  na <- rowSums(matrix(popcount(A), nrow(A))); nb <- rowSums(matrix(popcount(B), nrow(B))); any_b <- uni > 0
  c(jaccard = mean(inter[any_b] / uni[any_b]), retention = mean(c(sum(inter) / sum(na), sum(inter) / sum(nb))),
    identical = mean(rowSums(A != B) == 0), bold_per_patient = mean(c(na, nb)))
}
# The bag refits only: bag 0 above is the card's check, not a replicate pair.
mp <- data.frame(method = "llr_full", route = "bootstrap", boot_id = BAGS, seed_id = 0L, draw_id = 0L, stringsAsFactors = FALSE)
CP <- attr_contrast_pairs(mp)
rows <- lapply(seq_len(nrow(CP)), function(q) {
  pm <- pair_metrics(ENC[[as.character(CP$boot_i[q])]], ENC[[as.character(CP$boot_j[q])]])
  data.frame(contrast = CP$code[q], held = CP$held[q], boot_i = CP$boot_i[q], boot_j = CP$boot_j[q], t(pm))
})
raw <- do.call(rbind, rows)
summ <- do.call(rbind, lapply(split(raw, paste(raw$contrast, raw$held)), function(d)
  do.call(rbind, lapply(c("jaccard", "retention", "identical", "bold_per_patient"), function(mt) data.frame(
    contrast = d$contrast[1], held = d$held[1], metric = mt, n_pairs = nrow(d), median = stats::median(d[[mt]]),
    p025 = unname(stats::quantile(d[[mt]], 0.025)), p975 = unname(stats::quantile(d[[mt]], 0.975)))))))

run <- new_run(if (SITE == "mimic") "boldvals" else "boldvalsext", cfg_t,
               note = sprintf("bold-value agreement across bag refits (%s); terms refitted per bag, checked against store %s", SITE, basename(STORE)))
save_table(run, summ, "bold_value_agreement")
save_table(run, raw, "bold_value_agreement_pairs")
save_table(run, CK, "job_checks", subdir = "diagnostics")
finalize_run(run, extra = list(site = SITE, store = STORE, t1_run = T1_RUN, flags = FLAGS, tol = TOL, n_bags = length(BAGS),
                                tombstoned = paste(TOMB, collapse = ","), cache = CACHE))
options(width = 200)
s2 <- summ; s2$value <- sprintf("%.3f [%.3f, %.3f]", s2$median, s2$p025, s2$p975)
print(s2[, c("contrast", "held", "metric", "n_pairs", "value")], row.names = FALSE)
cat("\nrun:", run$path, "\n")
