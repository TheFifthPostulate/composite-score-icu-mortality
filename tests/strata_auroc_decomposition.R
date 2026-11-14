# tests/strata_auroc_decomposition.R -------------------------------------------
# HOW MUCH OF THE COMPOSITE SCORE'S DISCRIMINATION THE FOUR REPORTED STRATA CARRY,
# AND HOW MUCH THE ORDERING WITHIN THEM ADDS. DESCRIPTIVE.
#
# Added 2026-11-11 at the author's request, for the main Results. The quartiles
# are intervals of the composite score S, so every death-survivor pair either
# lies within one quartile or across two, and a pair across two is ordered by
# its quartiles alone. Hence, exactly,
#
#   AUROC(S) = AUROC(quartile) + pi_w * (A_w - 0.5),
#
# where AUROC(quartile) ranks the stays by their quartile alone (ties as one
# half), pi_w is the share of death-survivor pairs that lie within a quartile,
# and A_w is the within-quartile AUROC averaged over those pairs (each quartile
# weighted by its deaths x survivors). The second term is what the ordering
# within the strata adds. Reported beside it: the shares of pairs within a
# quartile, across quartiles with the death higher, and across with the death
# lower; the share of the above-chance discrimination the quartiles carry,
# (AUROC(quartile) - 0.5) / (AUROC(S) - 0.5); and the AUPRC lift of the
# quartile alone against that of S, and (added 2026-11-11) the share of the lift
# of S above one that the quartile ranking reaches, (lift_q - 1) / (lift_S - 1),
# with its interval from the same resamples (AUPRC has no exact pairwise split; lift =
# AUPRC / site mortality, interval bounds over the same mortality, as
# tests/hetero_card_flags.R and paper/make_results_figures.R).
# AUROC and AUPRC are .auroc() and .auprc() (R/09_metrics.R, ties as one
# threshold). Intervals: patient bootstrap, config n_boot and seed; the quartile
# of a stay is fixed (MIMIC-IV cut points, frozen for eICU).
#
# CHECKS (no run is written on a failure; CLAUDE.md rule 12):
#   the quartiles are intervals of S (no overlap between neighbours);
#   AUROC(quartile) from the pair counts equals .auroc() of the quartile;
#   the identity above holds to 1e-9;
#   the AUPRC gained within the quartiles, split by quartile (within_gain_by_quartile(), added 2026-11-11),
#   adds up to the whole gain to 1e-9;
#   AUROC(S) and AUPRC(S) reproduce the stored run of record (internal run,
#   llr_sum_oof, at MIMIC-IV; external run, llr_sum_frozen, at eICU).
#
# Reads the stage-1 stay tables only; fits nothing. AGGREGATE OUTPUT ONLY.
#
#   Rscript tests/strata_auroc_decomposition.R
#   Rscript tests/strata_auroc_decomposition.R --stage1 out/runs/heteroasm_<id> --internal out/runs/internal_<id> --external out/runs/external_<id>
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(yaml) })
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)
source("tests/hetero_common.R")

args     <- commandArgs(trailingOnly = TRUE)
STAGE1   <- .hopt(args, "--stage1",   "out/runs/heteroasm_20261027T233748")
INTERNAL <- .hopt(args, "--internal", "out/runs/internal_20260909T112643")
EXTERNAL <- .hopt(args, "--external", "out/runs/external_20260909T112824")
REF <- list(mimic_oof = list(file = file.path(INTERNAL, "tables", "score_summary_oof.csv"), label = "llr_sum_oof"),
            eicu      = list(file = file.path(EXTERNAL, "tables", "score_summary.csv"),     label = "llr_sum_frozen"))

cfg    <- load_config("config/config.yml")
N_BOOT <- as.integer(cfg_req(cfg, "metrics", "n_boot"))
SEED   <- as.integer(cfg_req(cfg, "seed"))

# One evaluation of the decomposition on rows (S, y, q).
decomp <- function(S, y, q) {
  y <- as.integer(y); lv <- sort(unique(q))
  d <- vapply(lv, function(k) sum(y[q == k]), 0); s <- vapply(lv, function(k) sum(q == k) - sum(y[q == k]), 0)
  A <- vapply(lv, function(k) { j <- q == k; if (sum(y[j]) == 0 || all(y[j] == 1)) NA_real_ else .auroc(S[j], y[j]) }, 0)
  pairs <- sum(d) * sum(s); w <- d * s
  above <- sum(vapply(seq_along(lv), function(i) d[i] * sum(s[seq_len(i - 1L)]), 0))
  below <- sum(vapply(seq_along(lv), function(i) d[i] * sum(s[-seq_len(i)]), 0))
  pi_w <- sum(w) / pairs
  A_w  <- sum((w * A)[w > 0]) / sum(w[w > 0])
  a_q  <- (above + 0.5 * sum(w)) / pairs
  a_s  <- .auroc(S, y)
  rate <- mean(y)
  c(pairs_within = pi_w, pairs_across_death_higher = above / pairs, pairs_across_death_lower = below / pairs,
    within_auroc = A_w, auroc_quartile = a_q, within_increment = pi_w * (A_w - 0.5), auroc_composite = a_s,
    quartile_share = (a_q - 0.5) / (a_s - 0.5),
    auprc_quartile = .auprc(q, y), auprc_composite = .auprc(S, y), mortality = rate,
    # added 2026-11-11: the share of the composite score's AUPRC lift above one that the quartile ranking reaches,
    # (lift_quartile - 1) / (lift_composite - 1), the AUPRC counterpart of quartile_share (not an exact split)
    quartile_share_auprc = (.auprc(q, y) - rate) / (.auprc(S, y) - rate),
    identity_gap = a_s - (a_q + pi_w * (A_w - 0.5)), quartile_auroc_direct = .auroc(q, y),
    within_gain_by_quartile(S, y, q, lv))
}

# Added 2026-11-11: the AUPRC gained by the ordering within the quartiles, AUPRC(S) - AUPRC(quartile), split by
# quartile. Erasing the order inside quartile k (its stays set to one tied value inside its own interval) changes
# only the precision at the deaths of quartile k: a stay of another quartile has the same stays above it whatever
# the order inside k. Hence the contributions add up to the whole gain exactly (checked below). Reported as the
# share of the gain that each quartile's ordering carries.
within_gain_by_quartile <- function(S, y, q, lv) {
  a_s <- .auprc(S, y)
  contrib <- vapply(lv, function(k) { S2 <- S; j <- q == k; S2[j] <- mean(S[j]); a_s - .auprc(S2, y) }, 0)
  gain <- a_s - .auprc(q, y)
  c(stats::setNames(contrib / gain, paste0("within_gain_share_q", lv)), within_gain_sum_gap = sum(contrib) - gain)
}

rows <- list(); chk <- list()
add_chk <- function(name, ok, detail) chk[[length(chk) + 1L]] <<- data.frame(check = name, detail = detail, status = if (isTRUE(ok)) "PASS" else "FAIL")
for (site in HETERO_SITES) {
  tab <- readRDS(file.path(STAGE1, "tables", paste0("hetero_stay_", site, ".rds")))
  need <- c("S", "y", "quartile", "group")
  if (!all(need %in% names(tab))) stop(site, ": stage-1 table lacks ", paste(setdiff(need, names(tab)), collapse = ", "), call. = FALSE)
  S <- tab$S; y <- as.integer(tab$y); q <- tab$quartile; g <- tab$group
  lv <- sort(unique(q))
  gaps <- vapply(seq_len(length(lv) - 1L), function(i) min(S[q == lv[i + 1L]]) - max(S[q == lv[i]]), 0)
  add_chk(sprintf("%s quartiles are intervals of the composite score", site), all(gaps > 0),
          sprintf("smallest gap between neighbouring quartiles %.3g", min(gaps)))
  est <- decomp(S, y, q)
  add_chk(sprintf("%s quartile AUROC from pair counts", site), abs(est[["auroc_quartile"]] - est[["quartile_auroc_direct"]]) < 1e-12,
          sprintf("%.6f vs %.6f", est[["auroc_quartile"]], est[["quartile_auroc_direct"]]))
  add_chk(sprintf("%s decomposition identity", site), abs(est[["identity_gap"]]) < 1e-9, sprintf("gap %.2e", est[["identity_gap"]]))
  add_chk(sprintf("%s AUPRC gain within the quartiles adds up over the quartiles", site), abs(est[["within_gain_sum_gap"]]) < 1e-9,
          sprintf("gap %.2e", est[["within_gain_sum_gap"]]))
  ref <- read.csv(REF[[site]]$file, stringsAsFactors = FALSE); ref <- ref[ref$label == REF[[site]]$label, ]
  add_chk(sprintf("%s stays and deaths match the run of record", site), nrow(ref) == 1L && ref$n == length(y) && ref$n_events == sum(y),
          sprintf("%d stays, %d deaths", length(y), sum(y)))
  add_chk(sprintf("%s AUROC matches the run of record", site), nrow(ref) == 1L && abs(ref$auroc - est[["auroc_composite"]]) < 1e-4,
          sprintf("%.5f vs %.5f", est[["auroc_composite"]], ref$auroc))
  add_chk(sprintf("%s AUPRC matches the run of record", site), nrow(ref) == 1L && abs(ref$auprc - est[["auprc_composite"]]) < 1e-4,
          sprintf("%.5f vs %.5f", est[["auprc_composite"]], ref$auprc))

  bu <- .boot_units(length(y), g)
  B <- with_seed(SEED, t(vapply(seq_len(N_BOOT), function(b) { j <- bu$draw(); decomp(S[j], y[j], q[j]) }, numeric(length(est)))))
  colnames(B) <- names(est)
  ci <- function(v) stats::quantile(v, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  keep <- c("pairs_within", "pairs_across_death_higher", "pairs_across_death_lower", "within_auroc", "auroc_quartile", "quartile_share_auprc",
            "within_increment", "auroc_composite", "quartile_share", grep("^within_gain_share_q", names(est), value = TRUE))
  r <- data.frame(site = site, stays = length(y), deaths = sum(y), mortality = est[["mortality"]])
  for (k in keep) { r[[k]] <- est[[k]]; r[[paste0(k, "_lo")]] <- ci(B[, k])[1]; r[[paste0(k, "_hi")]] <- ci(B[, k])[2] }
  for (k in c("quartile", "composite")) {
    a <- paste0("auprc_", k); r[[paste0("lift_", k)]] <- est[[a]] / est[["mortality"]]
    r[[paste0("lift_", k, "_lo")]] <- ci(B[, a])[1] / est[["mortality"]]; r[[paste0("lift_", k, "_hi")]] <- ci(B[, a])[2] / est[["mortality"]]
  }
  r$boot_unit <- bu$unit; r$n_boot <- N_BOOT
  rows[[site]] <- r
}
tab_out <- do.call(rbind, rows); chk <- do.call(rbind, chk)

options(width = 220)
show <- c("site", "stays", "deaths", "pairs_within", "pairs_across_death_higher", "pairs_across_death_lower", "within_auroc",
          "auroc_quartile", "auroc_quartile_lo", "auroc_quartile_hi", "within_increment", "within_increment_lo", "within_increment_hi",
          "auroc_composite", "quartile_share", "quartile_share_lo", "quartile_share_hi",
          "lift_quartile", "lift_quartile_lo", "lift_quartile_hi", "lift_composite", "lift_composite_lo", "lift_composite_hi",
          "quartile_share_auprc", grep("^within_gain_share_q", names(tab_out), value = TRUE))
r4 <- function(d) { d[] <- lapply(d, function(v) if (is.numeric(v) && !is.integer(v)) round(v, 4) else v); d }
cat("== decomposition of the composite score's discrimination into the strata and the ordering within them\n")
print(t(r4(tab_out[, show])))
cat("\n== checks\n"); print(chk, row.names = FALSE)
if (any(chk$status != "PASS")) stop("a check failed; no run written", call. = FALSE)

run <- new_run("stratadecomp", cfg, note = "AUROC of the composite score split into the quartiles and the ordering within them; descriptive, fits nothing")
save_table(run, chk, "strata_decomposition_checks", subdir = "diagnostics")
save_table(run, tab_out, "strata_decomposition")
finalize_run(run, extra = list(stage1_run = STAGE1, internal_run = INTERNAL, external_run = EXTERNAL))
cat("\nrun:", run$path, "\n")
