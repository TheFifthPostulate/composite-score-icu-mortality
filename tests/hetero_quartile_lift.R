# tests/hetero_quartile_lift.R -----------------------------------------------------
# THE SCORE AND THE RESIDUAL BOLD MARKS WITHIN EACH REPORTED QUARTILE: AUROC AND
# AUPRC LIFT, DEATHS AGAINST SURVIVORS. DESCRIPTIVE.
#
# Added 2026-11-11 at the author's request, elevated from a scratch check, for the
# Results table of the bold marks within a stratum (tab:quartile_lift). The run of
# tests/hetero_card_flags.R stores the AUROC and AUPRC lift of the residual bold
# marks (res_val, the excess of bold values over the count expected from the
# channels measured and the hours observed); tests/hetero_bold_vs_score.R stores
# the AUROC of the score S but not its AUPRC. This script computes both for S and
# res_val on the same stays, through the same call:
#   score_metrics() (R/09_metrics.R), patient bootstrap with config n_boot and seed,
#   AUPRC lift = AUPRC / mortality of the quartile, interval bounds divided by the
#   same mortality (as tests/hetero_card_flags.R).
#
# CHECK: res_val reproduces the stored contrast table of each heteroflags run
# (AUROC, lift and its interval), and the AUROC of S reproduces the stored
# bold_vs_score_auroc table. On a failure the script stops before writing a run:
# then identify whether the per-stay tables, this code or the stored runs differ
# (CLAUDE.md rule 12).
#
# Reads the per-stay tables only; fits nothing. Quartiles with fewer than 10
# deaths or 10 survivors are not computed. AGGREGATE OUTPUT ONLY.
#
#   Rscript tests/hetero_quartile_lift.R
#   Rscript tests/hetero_quartile_lift.R --mimic out/runs/heteroflags_<id> --eicu out/runs/heteroflags_<id> --boldscore out/runs/heteroboldscore_<id>
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(yaml) })
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)
source("tests/hetero_common.R")

args <- commandArgs(trailingOnly = TRUE)
RUNS <- c(mimic_oof = .hopt(args, "--mimic", "out/runs/heteroflags_20261030T204724"),
          eicu      = .hopt(args, "--eicu",  "out/runs/heteroflags_20261030T204938"))
BOLDSCORE <- .hopt(args, "--boldscore", "out/runs/heteroboldscore_20261101T115051")

# The eICU run must have used the MIMIC-IV run's frozen models.
.ex <- yaml::read_yaml(file.path(RUNS[["eicu"]], "manifest.yml"))$expect
if (normalizePath(.ex, mustWork = FALSE) != normalizePath(RUNS[["mimic_oof"]], mustWork = FALSE)) {
  stop("the eICU run was not built on the named MIMIC-IV run", call. = FALSE)
}
# The score AUROC check reads the bold-versus-score run, so it must have been
# built on the same two runs (guard added 2026-11-11).
.bs <- yaml::read_yaml(file.path(BOLDSCORE, "manifest.yml"))
for (.s in c(mimic_oof = "mimic_run", eicu = "eicu_run")) {
  .site <- names(which(c(mimic_oof = "mimic_run", eicu = "eicu_run") == .s))
  if (is.null(.bs[[.s]]) || normalizePath(.bs[[.s]], mustWork = FALSE) != normalizePath(RUNS[[.site]], mustWork = FALSE)) {
    stop("the bold-versus-score run ", BOLDSCORE, " was not built on ", RUNS[[.site]], call. = FALSE)
  }
}
cfg    <- load_config("config/config.yml")
N_BOOT <- as.integer(cfg_req(cfg, "metrics", "n_boot"))
SEED   <- as.integer(cfg_req(cfg, "seed"))
REF_S  <- read.csv(file.path(BOLDSCORE, "tables", "bold_vs_score_auroc.csv"), stringsAsFactors = FALSE)

rows <- list(); chk <- list()
# A reference value that is missing or duplicated in the stored table is a named
# FAIL, not an error from data.frame() (guard added 2026-11-11).
add_chk <- function(name, expected, got, tol = 1e-4) {
  one <- length(expected) == 1L && !is.na(expected)
  chk[[length(chk) + 1L]] <<- data.frame(check = name, expected = if (one) expected else NA_real_, got = got,
    status = if (!one) "FAIL (reference row missing or duplicated)" else if (isTRUE(abs(got - expected) < tol)) "PASS" else "FAIL")
}
for (site in names(RUNS)) {
  fl  <- readRDS(file.path(RUNS[[site]], "tables", paste0("card_flags_stay_", site, ".rds")))
  ref <- read.csv(file.path(RUNS[[site]], "tables", paste0("card_flag_contrast_", site, ".csv")), stringsAsFactors = FALSE)
  for (q in sort(unique(fl$quartile))) {
    i <- which(fl$quartile == q); y <- fl$y[i]; nd <- sum(y)
    if (nd < HETERO_MIN_CELL || length(i) - nd < HETERO_MIN_CELL) next
    s <- score_metrics(fl$S[i], y, label = "S", n_boot = N_BOOT, seed = SEED, group = fl$group[i])
    r <- score_metrics(fl$res_val[i], y, label = "res_val", n_boot = N_BOOT, seed = SEED, group = fl$group[i])
    rows[[length(rows) + 1L]] <- data.frame(site = site, quartile = q, stays = length(i), deaths = nd, mortality = s$event_rate,
      score_auroc = s$auroc, score_auroc_lo = s$auroc_lo, score_auroc_hi = s$auroc_hi,
      score_lift = s$auprc / s$event_rate, score_lift_lo = s$auprc_lo / s$event_rate, score_lift_hi = s$auprc_hi / s$event_rate,
      residual_auroc = r$auroc, residual_auroc_lo = r$auroc_lo, residual_auroc_hi = r$auroc_hi,
      residual_lift = r$auprc / r$event_rate, residual_lift_lo = r$auprc_lo / r$event_rate, residual_lift_hi = r$auprc_hi / r$event_rate,
      boot_unit = s$boot_unit, n_boot = s$n_boot)
    r0 <- ref[ref$partition == "quartile" & ref$stratum == q & ref$measure == "res_val", ]
    add_chk(sprintf("%s Q%d residual AUROC", site, q), r0$auroc, r$auroc)
    add_chk(sprintf("%s Q%d residual lift", site, q), r0$auprc_lift, r$auprc / r$event_rate)
    add_chk(sprintf("%s Q%d residual lift lower", site, q), r0$lift_lo, r$auprc_lo / r$event_rate)
    add_chk(sprintf("%s Q%d residual lift upper", site, q), r0$lift_hi, r$auprc_hi / r$event_rate)
    s0 <- REF_S[REF_S$site == site & REF_S$quartile == q, ]
    add_chk(sprintf("%s Q%d score AUROC", site, q), s0$auroc_score, s$auroc)
  }
}
tab <- do.call(rbind, rows); chk <- do.call(rbind, chk)

options(width = 200)
r3 <- function(d) { d[] <- lapply(d, function(v) if (is.numeric(v) && !is.integer(v)) round(v, 3) else v); d }
cat("== within each quartile: the score and the residual bold marks, deaths against survivors\n")
print(r3(tab[, setdiff(names(tab), c("boot_unit", "n_boot"))]), row.names = FALSE)
cat("\n== checks against the stored runs\n"); print(chk, row.names = FALSE)
if (any(chk$status != "PASS")) stop("the stored runs are not reproduced; no run written", call. = FALSE)

run <- new_run("heteroquartlift", cfg, note = "score and residual bold marks within quartiles: AUROC and AUPRC lift; descriptive, fits nothing")
save_table(run, chk, "quartile_lift_checks", subdir = "diagnostics")
save_table(run, tab, "quartile_lift")
finalize_run(run, extra = list(mimic_run = RUNS[["mimic_oof"]], eicu_run = RUNS[["eicu"]], boldscore_run = BOLDSCORE))
cat("\nrun:", run$path, "\n")
