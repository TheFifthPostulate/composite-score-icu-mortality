# tests/rv_assumptions_summary.R --------------------------------------------------
# ONE TABLE OF THE LIKELIHOOD-MODEL CHECKS, AT BOTH SITES, FOR THE PAPER'S
# APPENDIX ON THE ASSUMPTIONS OF THE RANDOM VARIABLES (sa:rv_assumptions).
#
# The code review of 2026-10-25 (docs/code_review_20261025.md, items 10 and 12)
# found that the reported counts were formed by hand: the 57 deviation cells
# include the four high tails that config/signal_tails declares empty, and the
# tail-extremity and fold-stability tables include the 12 tail extremities
# that enter no model. This script applies those two rules in code:
#   deviation categories  drop side "high" of a signal whose signal_tails entry
#                         omits it (19 x 3 - 4 = 53 categories)
#   tail extremities      keep rows with modeled == TRUE (26 of 38)
#   fold stability        keep rows with modeled == TRUE
# and summarizes each check per site. MIMIC-IV: the run-of-record internal
# diagnostics (final priors on the training stays). eICU: the external run's
# diagnostics of the same frozen priors applied to the eICU stays.
#
# Reads aggregate diagnostic tables only; fits nothing. AGGREGATE OUTPUT ONLY.
#
#   Rscript tests/rv_assumptions_summary.R [--internal out/runs/internal_...] [--external out/runs/external_...]
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(yaml) })
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)
args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
INT <- .opt("--internal", "out/runs/internal_20261025T211111")
EXT <- .opt("--external", "out/runs/external_20261025T211621")
cfg <- load_config("config/config.yml")
rd  <- function(dir, f) read.csv(file.path(dir, "diagnostics", f), stringsAsFactors = FALSE)

# --- declared rules ------------------------------------------------------------
declared_empty <- function(p) {
  st <- cfg$signal_tails
  vapply(seq_len(nrow(p)), function(i) {
    s <- p$signal[i]; c <- p$category[i]
    !is.null(st[[s]]) && c %in% c("low", "high") && !(c %in% unlist(st[[s]]))
  }, logical(1))
}
ppc_i <- rd(INT, "dm_fit_ppc.csv");     ppc_i <- ppc_i[!declared_empty(ppc_i), ]
ppc_e <- rd(EXT, "prior_dm_ppc.csv");   ppc_e <- ppc_e[!declared_empty(ppc_e), ]
if (nrow(ppc_i) != 53L || nrow(ppc_e) != 53L) stop("expected 53 modeled deviation categories at each site", call. = FALSE)
del_i <- rd(INT, "delta_fit.csv")
del_e <- rd(EXT, "prior_delta_fit.csv")
mod   <- del_i[del_i$modeled, c("signal", "variable")]
del_i <- merge(mod, del_i); del_e <- merge(mod, del_e)
if (nrow(del_i) != 26L || nrow(del_e) != 26L) stop("expected 26 modeled tail extremities", call. = FALSE)
gcs_i <- rd(INT, "delta_ordinal_fit.csv"); gcs_e <- rd(EXT, "prior_delta_ordinal_fit.csv")
gmod  <- mod[mod$signal %in% gcs_i$signal, ]
gcs_i <- merge(gmod, gcs_i); gcs_e <- merge(gmod, gcs_e)
lam_i <- rd(INT, "lambda_fit.csv"); lam_e <- rd(EXT, "prior_lambda_fit.csv")
shr_i <- rd(INT, "dm_shrinkage.csv"); shr_e <- rd(EXT, "prior_dm_shrinkage.csv")
dr_i  <- rd(INT, "dm_drift.csv"); dr_e <- rd(EXT, "prior_dm_drift.csv")
stab  <- rd(INT, "prior_stability.csv"); stab <- stab[stab$modeled, ]

# --- summaries -------------------------------------------------------------------
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
mr  <- function(x, d = 2) sprintf("%s [%s, %s]", fmt(median(x, na.rm = TRUE), d), fmt(min(x, na.rm = TRUE), d), fmt(max(x, na.rm = TRUE), d))
ks_max_mimic <- max(ppc_i$pit_ks)
sat <- function(p) { j <- which.max(p$pfull_obs - p$pfull_exp); sprintf("%s %s: %s%% against %s%%", p$signal[j], p$category[j], fmt(100 * p$pfull_obs[j], 1), fmt(100 * p$pfull_exp[j], 1)) }
amt <- function(l) { l <- l[l$family != "binomial", ]; paste(sprintf("%s %s", l$intervention, fmt(l$resid_mean, 2)), collapse = "; ") }
rows <- list(
  c("deviation counts", "dispersion of the counts, median [range] over the 53 categories (1 under the model)", mr(ppc_i$dispersion), mr(ppc_e$dispersion)),
  c("deviation counts", "Kolmogorov-Smirnov distance of the quantile residuals, median [range]", mr(ppc_i$pit_ks, 3), mr(ppc_e$pit_ks, 3)),
  c("deviation counts", sprintf("categories with a distance above the largest at MIMIC-IV (%s)", fmt(ks_max_mimic, 4)),
    as.character(sum(ppc_i$pit_ks > ks_max_mimic)), as.character(sum(ppc_e$pit_ks > ks_max_mimic))),
  c("deviation counts", "categories with a distance above 0.10", as.character(sum(ppc_i$pit_ks > 0.10)), as.character(sum(ppc_e$pit_ks > 0.10))),
  c("deviation counts", "largest excess of stays deviating in every observed hour, observed against expected", sat(ppc_i), sat(ppc_e)),
  c("deviation counts", "dispersion at high coverage over low coverage, median [range] over tails with a ratio",
    mr(dr_i$ratio), mr(dr_e$ratio)),
  c("deviation counts", "tails without a coverage ratio (coverage bins collapse)", as.character(sum(is.na(dr_i$ratio))), as.character(sum(is.na(dr_e$ratio)))),
  c("prior strength", "weight on the stay's own hours at the median observed hours, median [range] over the 19 measurements",
    mr(shr_i$w_med), mr(shr_e$w_med)),
  c("prior strength", "weight on the stay's own hours with one observed hour, median [range]", mr(shr_i$w_at_n1), mr(shr_e$w_at_n1)),
  c("tail extremity", "R2 of a smooth of the residual on the predictors, median [range] over the 26 modeled tails",
    mr(del_i$resid_smooth_r2, 3), mr(del_e$resid_smooth_r2, 3)),
  c("tail extremity", "modeled tails without the R2 check (fewer stays than the check's floor)", as.character(sum(is.na(del_i$resid_smooth_r2))), as.character(sum(is.na(del_e$resid_smooth_r2)))),
  c("Glasgow tail", "largest calibration difference of the proportional-odds model over the three components", fmt(max(gcs_i$cal_max_abs_diff), 3), fmt(max(gcs_e$cal_max_abs_diff), 3)),
  c("Glasgow tail", "Kolmogorov-Smirnov distance of the mid-probability transform, largest", fmt(max(gcs_i$pit_ks), 3), fmt(max(gcs_e$pit_ks), 3)),
  c("intensity, amounts", "mean residual (0 under the model)", amt(lam_i), amt(lam_e)),
  c("intensity, agent counts", "residual check", "not computed by the diagnostics", "not computed by the diagnostics"),
  c("fold stability", "coefficient of variation over the five folds, median [range]: Dirichlet concentration", mr(abs(stab$cv[stab$block == "alpha"]), 3), "-"),
  c("fold stability", "the same: tail-extremity coefficients", mr(abs(stab$cv[stab$block == "delta"]), 3), "-"),
  c("fold stability", "the same: intensity coefficients", mr(abs(stab$cv[stab$block == "lambda"]), 3), "-"),
  c("fold stability", "parameters whose sign differs between folds", sprintf("%d of %d", sum(stab$sign_flip, na.rm = TRUE), nrow(stab)), "-"))
tab <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(tab) <- c("model", "check", "mimic", "eicu")

cfg_run <- cfg
run <- new_run("rvassume", cfg_run, note = "likelihood-model checks at both sites, with the declared-empty tails and unmodeled extremities removed in code; fits nothing")
save_table(run, tab, "rv_assumption_checks")
save_table(run, ppc_i, "dm_ppc_modeled_mimic", subdir = "diagnostics")
save_table(run, ppc_e, "dm_ppc_modeled_eicu", subdir = "diagnostics")
finalize_run(run, extra = list(internal = INT, external = EXT))
options(width = 250)
print(tab, row.names = FALSE)
cat("\neICU categories above the MIMIC-IV maximum distance:\n")
z <- ppc_e[ppc_e$pit_ks > ks_max_mimic, c("signal", "category", "pit_ks")]; print(z[order(-z$pit_ks), ], row.names = FALSE)
cat("\nrun:", run$path, "\n")
