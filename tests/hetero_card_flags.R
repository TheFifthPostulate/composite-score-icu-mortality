# tests/hetero_card_flags.R ----------------------------------------------------
# Heterogeneity stage T2 (docs/heterogeneity_plan_20261027.md, section 8b).
# DO THE PATIENT CARD'S BOLD MARKS SEPARATE DEATHS FROM SURVIVORS WITHIN THE
# REPORTED QUARTILE, BEYOND WHAT THE MEASUREMENT PATTERN EXPLAINS?
#
# Two flags, counted exactly as paper/make_card_figure.R draws them:
#   bold percentile  a channel's weight of evidence at or above the 85th
#                    percentile (rounded, as displayed) of its out-of-fold
#                    distribution among measured MIMIC-IV training stays;
#                    counted per SYSTEM (config/domains.csv): k_sys flagged of
#                    N_sys systems with a measured channel.
#   bold value       a displayed random-variable contribution of at least log 2
#                    nats against a quiet day on the channel (R/15_anchors.R);
#                    EVERY bold item counts: k_val of N_val displayed items.
#                    Items, as on the card: each modeled propensity side; each
#                    modeled tail, or one combined tail for a single observed
#                    hour on two sides ("singular"); the trend when established;
#                    each paired intervention shown (its terms summed).
# Terms come from stage T1 (tests/hetero_oof_terms.R): out of fold at MIMIC-IV,
# the frozen final models at eICU.
#
# THE EXPECTATION MODEL (as for delta). Within each quartile, a binomial model
# of k out of N, logit link, fitted WITHOUT THE OUTCOME on MIMIC-IV out of fold
# and frozen for eICU:
#   bold percentiles  ~ log(m) [+ S]               m = measured channels, S = score
#   bold values       ~ log(m) + log(observed hours over measured channels) [+ S]
# The version without S is primary for claim A (author, 2026-10-28: compare
# within the reported stratum, which the reader sees); the version with S is
# secondary and shows how much of the signal the score position already carries.
# The residual k - N * p_hat is the stay's excess of flags.
#
# COMPARISONS, deaths against survivors within each quartile (primary) and
# decile: AUROC of the residual and of the raw count, minus 0.5, with a
# patient-bootstrap interval (score_metrics()); mortality by raw count (0, 1,
# 2, 3+) within bands of m. Confirmation (section 5): same sign at eICU, with
# the interval excluding zero. Cells under 10 stays show counts only.
#
# CHECK: the counting functions reproduce the bold marks of a rendered card
# (the high-risk demo card): 10 bold values, 5 bold systems.
#
# Fits only the outcome-blind expectation models. AGGREGATE OUTPUT ONLY.
#
#   Rscript tests/hetero_card_flags.R --site mimic_oof --terms out/runs/oofterms_...
#   Rscript tests/hetero_card_flags.R --site eicu --terms out/runs/oofterms_... --expect out/runs/heteroflags_<mimic run>
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(yaml); library(targets) })
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)
source("tests/hetero_common.R")

args   <- commandArgs(trailingOnly = TRUE)
SITE   <- .hopt(args, "--site", "mimic_oof")
STAGE1 <- .hopt(args, "--stage1", "out/runs/heteroasm_20261027T233748")
TERMS  <- .hopt(args, "--terms", NULL)
EXPECT <- .hopt(args, "--expect", NULL)
CARD   <- .hopt(args, "--card", "out/runs/card_20261028T120452")
if (is.null(TERMS)) stop("pass --terms <stage T1 run>", call. = FALSE)
if (SITE == "eicu" && is.null(EXPECT)) stop("eICU needs --expect <the MIMIC-IV run of this script>", call. = FALSE)
PCT_BOLD <- 85            # paper/make_card_figure.R
STRONG   <- log(2)        # paper/make_card_figure.R
M_BANDS  <- c(-Inf, 14, 16, 17, 19)   # the missingness bands: <=14, 15-16, 17, 18-19
M_LAB    <- c("<=14", "15-16", "17", "18-19")

cfg    <- load_config("config/config.yml")
N_BOOT <- as.integer(cfg_req(cfg, "metrics", "n_boot"))
SEED   <- as.integer(cfg_req(cfg, "seed"))
DOM    <- read.csv("config/domains.csv", stringsAsFactors = FALSE)

# --- counting, as the card draws ----------------------------------------------
# channel_items(), the counting of displayed and bold items, lives in
# tests/hetero_common.R (moved 2026-10-28 so tests/bold_values_replicates.R
# counts with the same code); it takes `strong` explicitly.
bold_pct <- function(pct) round(pct) >= PCT_BOLD
sys_count <- function(B, Mm) {          # B: n x 19 bold flags, Mm: n x 19 measured
  doms <- unique(DOM$domain)
  sb <- sapply(doms, function(d) rowSums(B[, DOM$signal[DOM$domain == d], drop = FALSE]) > 0)
  sm <- sapply(doms, function(d) rowSums(Mm[, DOM$signal[DOM$domain == d], drop = FALSE]) > 0)
  cbind(k_sys = rowSums(matrix(sb, nrow(B))), N_sys = rowSums(matrix(sm, nrow(B))))
}

# --- check against a rendered card --------------------------------------------
check_card <- function() {
  ps <- read.csv(file.path(CARD, "tables", "per_signal_high.csv"), stringsAsFactors = FALSE)
  ct <- read.csv(file.path(CARD, "tables", "card_terms_high.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  n_val <- 0L
  for (i in which(ps$measured)) {
    sg <- ps$signal[i]; z <- ct[ct$signal == sg, ]
    tt <- as.data.frame(as.list(stats::setNames(z$contribution, z$term)), check.names = FALSE)
    if (!is.na(ps$exposures[i]) && nzchar(ps$exposures[i])) for (e in strsplit(ps$exposures[i], ";", fixed = TRUE)[[1]]) {
      f <- strsplit(e, "|", fixed = TRUE)[[1]]
      if (f[2] != "absent") tt[[paste0("x:", f[1], "__", f[2])]] <- as.numeric(f[3])
      if (!is.na(f[4]) && f[4] != "NA") tt[[paste0("x:", f[1], "__present_at_admission")]] <- as.numeric(f[4])
    }
    n_val <- n_val + channel_items(tt, sg, ps$n_obs[i], as.logical(ps$trend_defined[i]), cfg, strong = STRONG,
                                   tails_eq = ps$tail_low[i] == ps$tail_high[i])[, "n_bold"]
  }
  B  <- matrix(ps$measured & bold_pct(ps$pct_train), 1, nrow(ps), dimnames = list(NULL, ps$signal))
  Mm <- matrix(ps$measured, 1, nrow(ps), dimnames = list(NULL, ps$signal))
  n_sys <- sys_count(B[, DOM$signal, drop = FALSE], Mm[, DOM$signal, drop = FALSE])[1, "k_sys"]
  data.frame(check = c("bold values on the high-risk demo card", "bold systems on the high-risk demo card"),
             expected = c(10L, 5L), got = c(n_val, n_sys), status = ifelse(c(n_val, n_sys) == c(10L, 5L), "PASS", "FAIL"))
}
chk <- check_card()
print(chk, row.names = FALSE)
if (any(chk$status != "PASS")) stop("the counting does not reproduce the rendered card", call. = FALSE)

# --- per-stay flags ---------------------------------------------------------------
H   <- hetero_read(STAGE1, SITE)
tab <- H$tab; Mm <- H$measured
staytab <- yaml::read_yaml(file.path(STAGE1, "manifest.yml"))$staytab
W   <- readRDS(file.path(staytab, "tables", paste0("woe_", SITE, ".rds")))
W_m <- readRDS(file.path(staytab, "tables", "woe_mimic_oof.rds")); M_m <- readRDS(file.path(staytab, "tables", "measured_mimic_oof.rds"))
sigs <- colnames(W)
if (!identical(sigs, colnames(Mm))) stop("channel order differs between stage A and stage 1", call. = FALSE)

# Bold percentiles: the reference is the MIMIC-IV out-of-fold distribution of
# measured stays, as the card's (run/patient_card.R: 100 * mean(ref <= x)).
B <- sapply(sigs, function(g) { e <- stats::ecdf(W_m[M_m[, g], g]); Mm[, g] & bold_pct(100 * e(W[, g])) })
SY <- sys_count(B[, DOM$signal], Mm[, DOM$signal])

# Bold values, from the stage T1 terms.
sf <- if (SITE == "mimic_oof") tar_read(tabs)$signal_features else {
  xcfg <- yaml::read_yaml("config/external.yml"); bundle <- load_bundle(cfg_req(xcfg, "bundle"), verbose = FALSE)
  cfg_e <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
  load_tables(cfg_e$paths, cfg_e, site = "eicu", verbose = FALSE)$signal_features
}
ids <- as.character(tab$stay_id)
k_val <- N_val <- hours <- numeric(nrow(tab))
for (sg in sigs) {
  tt <- readRDS(file.path(TERMS, "tables", sprintf("terms_%s_%s.rds", SITE, sg)))
  r  <- match(as.character(tt$stay_id), ids)
  if (anyNA(r)) stop("terms of ", sg, " include stays outside the stay table", call. = FALSE)
  if (!all(Mm[r, sg]) || sum(Mm[, sg]) != nrow(tt)) stop("terms of ", sg, " do not cover exactly the measured stays", call. = FALSE)
  if (max(abs(tt$L - W[r, sg])) > 1e-8) stop("terms of ", sg, " carry a different L from stage A", call. = FALSE)
  z  <- sf[sf$signal == sg, c("stay_id", "n_obs", "dx_trend_defined", "q05", "q95", "value_min", "value_max")]
  zz <- z[match(as.character(tt$stay_id), as.character(z$stay_id)), ]
  ci <- channel_items(tt, sg, zz$n_obs, as.logical(zz$dx_trend_defined), cfg, strong = STRONG,
                      tails_eq = tails_equal(zz, sg, cfg))
  k_val[r] <- k_val[r] + ci[, "n_bold"]; N_val[r] <- N_val[r] + ci[, "n_items"]; hours[r] <- hours[r] + zz$n_obs
}
rm(sf)
fl <- data.frame(stay_id = tab$stay_id, group = tab$group, y = tab$y, quartile = tab$quartile, decile = tab$decile,
                 S = tab$S, m = tab$m, hours = hours, k_sys = SY[, "k_sys"], N_sys = SY[, "N_sys"], k_val = k_val, N_val = N_val)
# A stay with no channel measured has a score of 0 and nothing the card could
# mark (20 eICU stays, none at MIMIC-IV). It is left out and counted.
N_NONE <- sum(fl$m == 0); N_NONE_DIED <- sum(fl$y[fl$m == 0])
fl <- fl[fl$m > 0, ]
if (any(fl$hours < 1)) stop("a measured stay with no observed hour", call. = FALSE)
cat(sprintf("stays with no channel measured, left out: %d (%d died)\n", N_NONE, N_NONE_DIED))

# --- the expectation models, within each quartile, without the outcome --------------
# Two versions (author, 2026-10-28). The first (primary, `sys` and `val`)
# adjusts for the measurement pattern only: the comparison is within the
# reported stratum, as the reader sees it. The second (secondary, `*_s`) also
# adjusts for the stay's score S: bold marks raise the score, and deaths sit
# high within their quartile, so it shows how much of the signal the score
# position already carries.
FORM <- list(sys   = cbind(k_sys, N_sys - k_sys) ~ log(m),
             val   = cbind(k_val, N_val - k_val) ~ log(m) + log(hours),
             sys_s = cbind(k_sys, N_sys - k_sys) ~ log(m) + S,
             val_s = cbind(k_val, N_val - k_val) ~ log(m) + log(hours) + S)
if (SITE == "mimic_oof") {
  coefs <- do.call(rbind, lapply(names(FORM), function(fn) do.call(rbind, lapply(sort(unique(fl$quartile)), function(q) {
    g <- stats::glm(FORM[[fn]], family = stats::binomial(), data = fl[fl$quartile == q, ])
    data.frame(flag = fn, quartile = q, term = names(stats::coef(g)), estimate = unname(stats::coef(g)))
  }))))
} else {
  coefs <- read.csv(file.path(EXPECT, "tables", "card_flag_expectation_coefs.csv"), stringsAsFactors = FALSE)
}
expected <- function(fn, d) {
  X <- stats::model.matrix(stats::update(FORM[[fn]], NULL ~ .), d)
  out <- numeric(nrow(d))
  for (q in unique(d$quartile)) {
    j <- d$quartile == q; cf <- coefs[coefs$flag == fn & coefs$quartile == q, ]
    out[j] <- plogis(X[j, cf$term, drop = FALSE] %*% cf$estimate)
  }
  out
}
fl$res_sys   <- fl$k_sys - fl$N_sys * expected("sys", fl)
fl$res_val   <- fl$k_val - fl$N_val * expected("val", fl)
fl$res_sys_s <- fl$k_sys - fl$N_sys * expected("sys_s", fl)
fl$res_val_s <- fl$k_val - fl$N_val * expected("val_s", fl)

# --- comparisons ---------------------------------------------------------------------
auc_rows <- list(); mort_rows <- list(); desc_rows <- list()
for (part in c("quartile", "decile")) {
  strata <- split(seq_len(nrow(fl)), fl[[part]])
  for (st in names(strata)) {
    i <- strata[[st]]; d <- sum(fl$y[i]); s <- length(i) - d
    for (yy in 0:1) {
      j <- i[fl$y[i] == yy]
      desc_rows[[length(desc_rows) + 1L]] <- data.frame(site = SITE, partition = part, stratum = as.integer(st),
        outcome = if (yy == 1) "died" else "survived", n = length(j),
        mean_k_sys = if (length(j) >= HETERO_MIN_CELL) mean(fl$k_sys[j]) else NA, mean_N_sys = if (length(j) >= HETERO_MIN_CELL) mean(fl$N_sys[j]) else NA,
        mean_k_val = if (length(j) >= HETERO_MIN_CELL) mean(fl$k_val[j]) else NA, mean_N_val = if (length(j) >= HETERO_MIN_CELL) mean(fl$N_val[j]) else NA,
        mean_m = if (length(j) >= HETERO_MIN_CELL) mean(fl$m[j]) else NA)
    }
    if (d < HETERO_MIN_CELL || s < HETERO_MIN_CELL) next
    for (meas in c("res_sys_s", "res_sys", "k_sys", "res_val_s", "res_val", "k_val")) {
      sm <- score_metrics(fl[[meas]][i], fl$y[i], label = meas, n_boot = N_BOOT, seed = SEED, group = fl$group[i])
      auc_rows[[length(auc_rows) + 1L]] <- data.frame(site = SITE, partition = part, stratum = as.integer(st), measure = meas,
        deaths = d, survivors = s, auroc = sm$auroc, difference = sm$auroc - 0.5, lo = sm$auroc_lo - 0.5, hi = sm$auroc_hi - 0.5,
        auprc_lift = sm$auprc / sm$event_rate, lift_lo = sm$auprc_lo / sm$event_rate, lift_hi = sm$auprc_hi / sm$event_rate,
        boot_unit = sm$boot_unit, n_boot = sm$n_boot, primary = part == "quartile" && meas %in% c("res_sys", "res_val"))
    }
  }
}
band <- cut(fl$m, M_BANDS, labels = M_LAB)
for (fn in c("k_sys", "k_val")) for (q in sort(unique(fl$quartile))) for (b in M_LAB) for (k in c("0", "1", "2", "3+")) {
  kk <- if (k == "3+") fl[[fn]] >= 3 else fl[[fn]] == as.integer(k)
  j <- fl$quartile == q & band == b & kk
  mort_rows[[length(mort_rows) + 1L]] <- data.frame(site = SITE, flag = fn, quartile = q, m_band = b, count = k, n = sum(j),
    deaths = if (sum(j) >= HETERO_MIN_CELL) sum(fl$y[j]) else NA, mortality = if (sum(j) >= HETERO_MIN_CELL) mean(fl$y[j]) else NA)
}
auc <- do.call(rbind, auc_rows); mort <- do.call(rbind, mort_rows); desc <- do.call(rbind, desc_rows)

run <- new_run("heteroflags", cfg, note = sprintf("heterogeneity T2 (%s): card flags within quartiles, expectation-adjusted", SITE))
save_table(run, chk, "card_flag_checks", subdir = "diagnostics")
save_table(run, auc, paste0("card_flag_contrast_", SITE))
save_table(run, mort, paste0("card_flag_mortality_", SITE))
save_table(run, desc, paste0("card_flag_summary_", SITE))
if (SITE == "mimic_oof") save_table(run, coefs, "card_flag_expectation_coefs")
save_object(run, fl, paste0("card_flags_stay_", SITE))
finalize_run(run, extra = list(site = SITE, stage1 = STAGE1, terms = TERMS, expect = EXPECT %||% "", card = CARD,
                                pct_bold = PCT_BOLD, strong = STRONG, n_no_channel = N_NONE, n_no_channel_died = N_NONE_DIED))

options(width = 200)
r3 <- function(d) { d[] <- lapply(d, function(v) if (is.numeric(v)) round(v, 3) else v); d }
cat("\n== mean flags and opportunities by quartile and outcome\n")
print(r3(desc[desc$partition == "quartile", -2]), row.names = FALSE)
if (SITE == "mimic_oof") { cat("\n== expectation models (outcome-blind)\n"); print(r3(coefs), row.names = FALSE) }
cat("\n== deaths against survivors within quartile: AUROC minus 0.5 (patient bootstrap)\n")
print(r3(auc[auc$partition == "quartile", c("stratum", "measure", "deaths", "auroc", "difference", "lo", "hi", "auprc_lift", "lift_lo", "lift_hi", "primary")]), row.names = FALSE)
cat("\n== same, by decile (residuals only)\n")
print(r3(auc[auc$partition == "decile" & auc$measure %in% c("res_sys_s", "res_val_s", "res_sys", "res_val"), c("stratum", "measure", "deaths", "difference", "lo", "hi")]), row.names = FALSE)
cat("\nrun:", run$path, "\n")
