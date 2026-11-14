# tests/strata_reporting.R -----------------------------------------------------
# Performance of the REPORTED quantity: the observed mortality of the frozen
# reference stratum a stay falls in. NON-FITTING CONSUMER. It reads only the
# per-bin tables (counts of stays and deaths in each of the 20 frozen bins) that
# the internal, test-look and external runs already exported, so it touches no
# stay-level data, fits nothing, and is not a new look at the test set: the
# test-set numbers are a post hoc description of the single look already spent
# (test_20260912T215546), and no design choice is made from them.
#
# Three analyses, for each arm and site, with quartiles (bins 1-5, ..., 16-20)
# and deciles (bins 1-2, ..., 19-20) as the reporting partitions:
#   1. strata_pr          precision and recall when the top k strata are flagged
#                         (the precision-recall curve read at the cut points),
#                         with Wilson intervals and lift over the site's rate
#   2. strata_calibration observed mortality of each stratum at the site against
#                         its reference (MIMIC-IV training, out of fold) rate,
#                         plus one summary row per arm, site and partition
#   3. strata_discrimination AUROC (with Hanley-McNeil interval) and stepwise
#                         average precision of the stratum ASSIGNMENT (ties within
#                         a stratum), beside the raw score's AUROC and AUPRC, so
#                         the discrimination given up by reporting strata is shown
#
#   Rscript tests/strata_reporting.R
#   Rscript tests/strata_reporting.R --internal out/runs/internal_... \
#     --test out/runs/test_... --external out/runs/external_...
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(yaml))
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
INTERNAL <- .opt("--internal", "out/runs/internal_20260909T112643")
TEST     <- .opt("--test", "out/runs/test_20260912T215546")
EXTERNAL <- .opt("--external", "out/runs/external_20260909T112824")
ARMS <- c("llr_sum", "xgb_feat", "xgb_raw")
# --ebm-bins <run>: also the EBM arms, from the bin tables tests/ebm_strata_bins.R
# writes in the graph's format (author, 2026-10-29). Those arms were scored
# outside the graph, so their tables live in that one run, not in the site runs.
EBM_BINS <- .opt("--ebm-bins", NULL)
# The arms come from the bins run itself, so the joint GAM arms (added to that
# run by tests/ebm_strata_bins.R --jointgam-*) are reported when present.
EBM_ARMS <- if (is.null(EBM_BINS)) c("ebm", "ebmadd") else
  unique(read.csv(file.path(EBM_BINS, "tables", "ebm_raw_summary.csv"), stringsAsFactors = FALSE)$arm)
if (!is.null(EBM_BINS)) ARMS <- c(ARMS, EBM_ARMS)
ebm_file <- function(arm, site) file.path(EBM_BINS, "tables", sprintf("risk_bins_%s_%s.csv", arm,
  switch(site, mimic_train_oof = "oof", mimic_test = "test_frozen", eicu = "eicu_frozen")))
# THE BIN COUNT IS DECLARED IN CONFIG, NOT HERE (code review 2026-10-25, C-2).
# `bin = 1:20` and the quartile/decile group sizes were literals, so a different
# `metrics.n_bins` would have grouped the wrong bins into "quartiles" without an
# error. The group sizes follow from the declared count, which must divide.
N_BINS <- as.integer(cfg_req(load_config("config/config.yml"), "metrics", "n_bins"))
if (N_BINS %% 4L != 0L || N_BINS %% 10L != 0L) {
  stop("metrics.n_bins = ", N_BINS, " does not divide into quartiles and deciles", call. = FALSE)
}
PARTITIONS <- list(quartile = N_BINS %/% 4L, decile = N_BINS %/% 10L)

read_bins <- function(run, file) {
  p <- if (is.null(run)) file else file.path(run, "tables", file)
  if (!file.exists(p)) stop("missing bin table: ", p, call. = FALSE)
  b <- read.csv(p, stringsAsFactors = FALSE)
  if (!all(c("bin", "n", "deaths") %in% names(b))) stop("unexpected columns in ", p, call. = FALSE)
  out <- data.frame(bin = seq_len(N_BINS), n = 0L, deaths = 0L)
  out$n[match(b$bin, out$bin)] <- b$n
  out$deaths[match(b$bin, out$bin)] <- b$deaths
  if (sum(out$n) != sum(b$n)) stop("bins outside 1-", N_BINS, " in ", p, call. = FALSE)
  out
}
wilson <- function(k, n, z = 1.96) {
  if (n == 0) return(c(NA_real_, NA_real_))
  p <- k / n; d <- 1 + z^2 / n; c0 <- (p + z^2 / (2 * n)) / d
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(max(0, c0 - h), min(1, c0 + h))
}
group <- function(b, size) {
  g <- (b$bin - 1L) %/% size + 1L
  data.frame(stratum = sort(unique(g)), n = as.integer(tapply(b$n, g, sum)), deaths = as.integer(tapply(b$deaths, g, sum)))
}
# AUROC of an ordinal assignment from counts per level: P(death level > survivor
# level) + 0.5 P(tie); Hanley-McNeil standard error.
auroc_counts <- function(n, d) {
  s <- n - d; D <- sum(d); S <- sum(s)
  below_s <- cumsum(c(0, s))[seq_along(s)]
  a <- sum(d * (below_s + 0.5 * s)) / (D * S)
  q1 <- a / (2 - a); q2 <- 2 * a^2 / (1 + a)
  se <- sqrt((a * (1 - a) + (D - 1) * (q1 - a^2) + (S - 1) * (q2 - a^2)) / (D * S))
  c(auroc = a, lo = a - 1.96 * se, hi = a + 1.96 * se)
}
# Stepwise average precision: flag strata from the top down; AP = sum over k of
# (recall_k - recall_{k-1}) * precision_k.
ap_counts <- function(n, d) {
  o <- rev(seq_along(n)); cn <- cumsum(n[o]); cd <- cumsum(d[o])
  rec <- cd / sum(d); prec <- ifelse(cn > 0, cd / cn, 0)
  sum(diff(c(0, rec)) * prec)
}

SITES <- list(
  list(site = "mimic_train_oof", run = INTERNAL, file = function(a) sprintf("risk_bins_%s_oof.csv", a)),
  list(site = "mimic_test",      run = TEST,     file = function(a) sprintf("risk_bins_%s_frozen.csv", a)),
  list(site = "eicu",            run = EXTERNAL, file = function(a) sprintf("risk_bins_%s_frozen.csv", a)))
raw_summary <- function(run, site, arm) {
  if (arm %in% EBM_ARMS) {
    s <- read.csv(file.path(EBM_BINS, "tables", "ebm_raw_summary.csv"), stringsAsFactors = FALSE)
    r <- s[s$arm == arm & s$site == site, ]
    if (nrow(r) != 1L) stop("no raw summary for ", arm, " at ", site, call. = FALSE)
    return(c(auroc = r$auroc, auprc = r$auprc))
  }
  f <- if (site == "mimic_train_oof") "score_summary_oof.csv" else "score_summary.csv"
  s <- read.csv(file.path(run, "tables", f), stringsAsFactors = FALSE)
  lab <- if (site == "mimic_train_oof") paste0(arm, "_oof") else paste0(arm, "_frozen")
  r <- s[s$label %in% c(lab, arm), , drop = FALSE]
  if (!nrow(r)) return(c(auroc = NA_real_, auprc = NA_real_))
  c(auroc = r$auroc[1], auprc = r$auprc[1])
}

pr_rows <- list(); cal_rows <- list(); cal_sum <- list(); disc_rows <- list()
for (arm in ARMS) {
  is_ebm <- arm %in% EBM_ARMS
  ref20 <- if (is_ebm) read_bins(NULL, ebm_file(arm, "mimic_train_oof")) else read_bins(INTERNAL, sprintf("risk_bins_%s_oof.csv", arm))
  for (st in SITES) {
    b20 <- if (is_ebm) read_bins(NULL, ebm_file(arm, st$site)) else read_bins(st$run, st$file(arm))
    N <- sum(b20$n); D <- sum(b20$deaths); base <- D / N
    rs <- raw_summary(st$run, st$site, arm)
    a20 <- auroc_counts(b20$n, b20$deaths)
    disc_rows[[length(disc_rows) + 1L]] <- data.frame(arm = arm, site = st$site, assignment = "raw score",
      auroc = rs[["auroc"]], auroc_lo = NA, auroc_hi = NA, ap = rs[["auprc"]], base_rate = base)
    disc_rows[[length(disc_rows) + 1L]] <- data.frame(arm = arm, site = st$site, assignment = sprintf("%d bins", N_BINS),
      auroc = a20[["auroc"]], auroc_lo = a20[["lo"]], auroc_hi = a20[["hi"]], ap = ap_counts(b20$n, b20$deaths), base_rate = base)
    for (pn in names(PARTITIONS)) {
      g <- group(b20, PARTITIONS[[pn]]); gr <- group(ref20, PARTITIONS[[pn]])
      ag <- auroc_counts(g$n, g$deaths)
      disc_rows[[length(disc_rows) + 1L]] <- data.frame(arm = arm, site = st$site, assignment = pn,
        auroc = ag[["auroc"]], auroc_lo = ag[["lo"]], auroc_hi = ag[["hi"]], ap = ap_counts(g$n, g$deaths), base_rate = base)
      K <- nrow(g)
      for (k in seq_len(K - 1L)) {          # flag the top k strata
        top <- (K - k + 1L):K
        tn <- sum(g$n[top]); td <- sum(g$deaths[top])
        pci <- wilson(td, tn); rci <- wilson(td, D)
        pr_rows[[length(pr_rows) + 1L]] <- data.frame(arm = arm, site = st$site, partition = pn, top_strata = k,
          share_stays = tn / N, recall = td / D, recall_lo = rci[1], recall_hi = rci[2],
          precision = td / tn, precision_lo = pci[1], precision_hi = pci[2], lift = (td / tn) / base, base_rate = base)
      }
      if (st$site != "mimic_train_oof") {
        obs <- g$deaths / g$n; ref <- gr$deaths / gr$n
        ci <- t(mapply(wilson, g$deaths, g$n))
        inside <- ref >= ci[, 1] & ref <= ci[, 2]
        cal_rows[[length(cal_rows) + 1L]] <- data.frame(arm = arm, site = st$site, partition = pn, stratum = g$stratum,
          n = g$n, share_stays = g$n / N, deaths = g$deaths, observed = obs, observed_lo = ci[, 1], observed_hi = ci[, 2],
          reference = ref, difference = obs - ref, ratio = obs / ref, reference_inside_ci = inside)
        w <- g$n / N
        cal_sum[[length(cal_sum) + 1L]] <- data.frame(arm = arm, site = st$site, partition = pn,
          max_abs_difference = max(abs(obs - ref)), weighted_mean_abs_difference = sum(w * abs(obs - ref)),
          strata_reference_outside_ci = sum(!inside), n_strata = length(obs),
          top_stratum_observed = obs[length(obs)], top_stratum_reference = ref[length(ref)])
      }
    }
  }
}
pr <- do.call(rbind, pr_rows); cal <- do.call(rbind, cal_rows); cals <- do.call(rbind, cal_sum); disc <- do.call(rbind, disc_rows)

cfg <- load_config("config/config.yml")
run <- new_run("strata", cfg, note = "performance of the reported stratum mortality; non-fitting; reads per-bin tables only")
save_table(run, pr, "strata_pr"); save_table(run, cal, "strata_calibration")
save_table(run, cals, "strata_calibration_summary"); save_table(run, disc, "strata_discrimination")
finalize_run(run, extra = list(internal = INTERNAL, test = TEST, external = EXTERNAL, ebm_bins = EBM_BINS %||% "", arms = paste(ARMS, collapse = ",")))

# Console: aggregate tables only.
options(width = 200)
r3 <- function(x) { x[] <- lapply(x, function(v) if (is.numeric(v)) round(v, 3) else v); x }
cat("\n== 1. precision-recall at the cut points (llr_sum)\n")
print(r3(pr[pr$arm == "llr_sum", c("site", "partition", "top_strata", "share_stays", "recall", "precision", "precision_lo", "precision_hi", "lift")]), row.names = FALSE)
cat("\n== 1b. top stratum only, all arms\n")
top1 <- pr[pr$top_strata == 1, c("arm", "site", "partition", "share_stays", "recall", "precision", "lift")]
print(r3(top1), row.names = FALSE)
cat("\n== 2. stratum calibration, quartiles\n")
print(r3(cal[cal$partition == "quartile", c("arm", "site", "stratum", "observed", "observed_lo", "observed_hi", "reference", "difference", "reference_inside_ci")]), row.names = FALSE)
cat("\n== 2b. calibration summary\n"); print(r3(cals), row.names = FALSE)
cat("\n== 3. discrimination of the assignment\n"); print(r3(disc), row.names = FALSE)
cat("\nrun:", run$path, "\n")
