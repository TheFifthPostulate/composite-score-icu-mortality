# tests/hetero_bold_vs_score.R ----------------------------------------------------
# THE CARD'S BOLD COUNT AGAINST THE SCORE WITHIN EACH REPORTED QUARTILE
# (docs/heterogeneity_plan_20261027.md, section 8b.3). DESCRIPTIVE.
#
# Section 8b.3 was a scratch check with no script; the Results and Methods cite
# it ("the bold count reaches 50 to 90% of the score's own discrimination above
# chance, and it adds no detection beyond the score"). This script reproduces it
# from the per-stay tables of the claim A runs. Within each quartile and site:
#   AUROC of the raw bold-value count k_val, and of the score S
#   retained share   (AUROC_bold - 0.5) / (AUROC_score - 0.5): how much of the
#                    score's own within-quartile ordering the bold count keeps
#   tie-broken rank  AUROC of a ranking by bold count with ties broken by the
#                    score (bold count first, then S)
#   same burden      for k = 1..6, the deaths caught by flagging the stays with
#                    at least k bold values, against the deaths caught by
#                    flagging the same number of stays with the highest S
#   within bands     (added 2026-11-01) the model-free check of claim A: the
#                    AUROC of the raw bold count within each declared band of
#                    measured channels (<=14, 15-16, 17, 18-19, as
#                    tests/hetero_card_flags.R), bands with at least 10 deaths
#                    and 10 survivors, averaged with stay weights; beside it
#                    the AUROC of the measurement-adjusted residual. Section
#                    8b.1 of the plan quotes a scratch version with other bands
#                    (MIMIC-IV quartile 1: 0.553); it is not a check here.
# AUROC is .auroc() (R/09_metrics.R), deaths against survivors, ties as one half.
#
# CHECK: the figures of section 8b.3 (quartiles 2 and 3, both sites, the
# printed example k; and the retained shares of all quartiles) are reproduced
# at their printed precision. On a failure the script stops before writing a
# run: then identify whether the per-stay tables, this code or the scratch
# figures differ (CLAUDE.md rule 12).
#
# Reads the per-stay tables only; fits nothing. Cells under 10 stays or 10
# deaths show counts only. AGGREGATE OUTPUT ONLY.
#
#   Rscript tests/hetero_bold_vs_score.R
#   Rscript tests/hetero_bold_vs_score.R --mimic out/runs/heteroflags_<id> --eicu out/runs/heteroflags_<id>
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(yaml) })
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)
source("tests/hetero_common.R")

args <- commandArgs(trailingOnly = TRUE)
RUNS <- c(mimic_oof = .hopt(args, "--mimic", "out/runs/heteroflags_20261030T204724"),
          eicu      = .hopt(args, "--eicu",  "out/runs/heteroflags_20261030T204938"))
K <- 1:6

# The eICU run must have used the MIMIC-IV run's frozen models.
.ex <- yaml::read_yaml(file.path(RUNS[["eicu"]], "manifest.yml"))$expect
if (normalizePath(.ex, mustWork = FALSE) != normalizePath(RUNS[["mimic_oof"]], mustWork = FALSE)) {
  stop("the eICU run was not built on the named MIMIC-IV run", call. = FALSE)
}
cfg <- load_config("config/config.yml")

M_BANDS <- c(-Inf, 14, 16, 17, 19)   # tests/hetero_card_flags.R

auc_rows <- list(); burden_rows <- list(); band_rows <- list()
for (site in names(RUNS)) {
  fl <- readRDS(file.path(RUNS[[site]], "tables", paste0("card_flags_stay_", site, ".rds")))
  band <- cut(fl$m, M_BANDS)
  for (q in sort(unique(fl$quartile))) {
    a <- w <- numeric(0); used <- 0L
    for (lv in levels(band)) {
      j <- fl$quartile == q & band == lv; dj <- sum(fl$y[j])
      if (dj >= HETERO_MIN_CELL && sum(j) - dj >= HETERO_MIN_CELL) { a <- c(a, .auroc(fl$k_val[j], fl$y[j])); w <- c(w, sum(j)); used <- used + 1L }
    }
    iq <- fl$quartile == q
    band_rows[[length(band_rows) + 1L]] <- data.frame(site = site, quartile = q, bands_used = used,
      stays_in_used_bands = sum(w), auroc_within_bands = if (length(a)) sum(a * w) / sum(w) else NA,
      auroc_residual = .auroc(fl$res_val[iq], fl$y[iq]))
  }
  for (q in sort(unique(fl$quartile))) {
    d <- fl[fl$quartile == q, ]; y <- d$y; n <- nrow(d); nd <- sum(y)
    if (nd < HETERO_MIN_CELL || n - nd < HETERO_MIN_CELL) {
      auc_rows[[length(auc_rows) + 1L]] <- data.frame(site = site, quartile = q, stays = n, deaths = nd,
        auroc_bold = NA, auroc_score = NA, retained_share = NA, auroc_bold_then_score = NA)
      next
    }
    a_b <- .auroc(d$k_val, y); a_s <- .auroc(d$S, y)
    key <- d$k_val + rank(d$S) / (n + 1)                    # bold count first, the score breaks ties
    auc_rows[[length(auc_rows) + 1L]] <- data.frame(site = site, quartile = q, stays = n, deaths = nd,
      auroc_bold = a_b, auroc_score = a_s, retained_share = (a_b - 0.5) / (a_s - 0.5),
      auroc_bold_then_score = .auroc(key, y))
    by_s <- order(-d$S)
    for (k in K) {
      fk <- d$k_val >= k; n_fl <- sum(fk)
      ok <- n_fl >= HETERO_MIN_CELL && n - n_fl >= HETERO_MIN_CELL
      top <- if (n_fl > 0) by_s[seq_len(n_fl)] else integer(0)
      burden_rows[[length(burden_rows) + 1L]] <- data.frame(site = site, quartile = q, k = k, stays = n, deaths = nd,
        n_flagged = n_fl, share_flagged = n_fl / n,
        caught_bold  = if (ok) sum(y[fk]) / nd else NA,
        caught_score = if (ok) sum(y[top]) / nd else NA)
    }
  }
}
auc <- do.call(rbind, auc_rows); bur <- do.call(rbind, burden_rows); bands <- do.call(rbind, band_rows)
bands$difference <- bands$auroc_within_bands - bands$auroc_residual
bur$score_minus_bold <- bur$caught_score - bur$caught_bold

# --- check against the scratch figures of section 8b.3 ---------------------------
REF <- data.frame(site = c("mimic_oof", "eicu", "mimic_oof", "eicu"), quartile = c(2L, 2L, 3L, 3L), k = c(2L, 2L, 3L, 3L),
                  auroc_bold = c(0.602, 0.592, 0.550, 0.600), auroc_score = c(0.644, 0.628, 0.596, 0.610),
                  auroc_bold_then_score = c(0.628, 0.616, 0.565, 0.612), share_flagged = c(0.57, 0.45, 0.75, 0.56),
                  caught_bold = c(0.753, 0.591, 0.805, 0.701), caught_score = c(0.760, 0.614, 0.844, 0.699))
SHARE <- data.frame(site = rep(c("mimic_oof", "eicu"), each = 4), quartile = rep(1:4, 2),
                    retained_share = c(0.31, 0.71, 0.52, 0.77, 0.69, 0.72, 0.91, 0.78))
near <- function(got, want, digits) isTRUE(abs(round(got, digits) - want) < 1e-9)
chk <- list()
for (i in seq_len(nrow(REF))) {
  r <- REF[i, ]; a <- auc[auc$site == r$site & auc$quartile == r$quartile, ]
  b <- bur[bur$site == r$site & bur$quartile == r$quartile & bur$k == r$k, ]
  for (v in c("auroc_bold", "auroc_score", "auroc_bold_then_score"))
    chk[[length(chk) + 1L]] <- data.frame(check = sprintf("%s Q%d %s", r$site, r$quartile, v), expected = r[[v]], got = round(a[[v]], 3),
                                          status = if (near(a[[v]], r[[v]], 3)) "PASS" else "FAIL")
  for (v in c("share_flagged", "caught_bold", "caught_score")) {
    dg <- if (v == "share_flagged") 2 else 3
    chk[[length(chk) + 1L]] <- data.frame(check = sprintf("%s Q%d k=%d %s", r$site, r$quartile, r$k, v), expected = r[[v]],
                                          got = round(b[[v]], dg), status = if (near(b[[v]], r[[v]], dg)) "PASS" else "FAIL")
  }
}
# The scratch shares were computed from AUROCs already rounded to three decimals
# (found 2026-11-01: only eICU quartile 3 differs at two decimals, 0.905 unrounded
# against 0.909 from rounded AUROCs). The check reproduces that computation; the
# saved table carries the unrounded share.
for (i in seq_len(nrow(SHARE))) {
  r <- SHARE[i, ]; a <- auc[auc$site == r$site & auc$quartile == r$quartile, ]
  sh3 <- (round(a$auroc_bold, 3) - 0.5) / (round(a$auroc_score, 3) - 0.5)
  chk[[length(chk) + 1L]] <- data.frame(check = sprintf("%s Q%d retained_share (from 3-decimal AUROCs)", r$site, r$quartile),
                                        expected = r$retained_share, got = round(sh3, 2), status = if (near(sh3, r$retained_share, 2)) "PASS" else "FAIL")
}
chk <- do.call(rbind, chk)

options(width = 200)
r3 <- function(d) { d[] <- lapply(d, function(v) if (is.numeric(v) && !is.integer(v)) round(v, 3) else v); d }
cat("== bold-value count against the score within each quartile\n"); print(r3(auc), row.names = FALSE)
cat("\n== deaths caught at the same burden: at least k bold values, or the same number of highest scores\n")
print(r3(bur[, c("site", "quartile", "k", "share_flagged", "caught_bold", "caught_score", "score_minus_bold")]), row.names = FALSE)
cat("\n== model-free check: raw bold count within bands of measured channels, against the adjusted residual\n")
print(r3(bands), row.names = FALSE)
cat("\n== check against docs/heterogeneity_plan_20261027.md section 8b.3\n"); print(chk, row.names = FALSE)
if (any(chk$status != "PASS")) stop("the scratch figures of section 8b.3 are not reproduced; no run written", call. = FALSE)

run <- new_run("heteroboldscore", cfg, note = "bold-value count against the score within quartiles (8b.3); descriptive, fits nothing")
save_table(run, chk, "bold_vs_score_checks", subdir = "diagnostics")
save_table(run, auc, "bold_vs_score_auroc")
save_table(run, bur, "bold_vs_score_same_burden")
save_table(run, bands, "bold_count_within_bands")
finalize_run(run, extra = list(mimic_run = RUNS[["mimic_oof"]], eicu_run = RUNS[["eicu"]]))
cat("\nrun:", run$path, "\n")
