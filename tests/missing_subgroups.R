# tests/missing_subgroups.R ----------------------------------------------------
# Analysis queue item 1a: does the score work as well for sparsely measured
# stays? NON-FITTING SUBGROUP CHECK of the reported quantity. Reads only the
# stage A per-stay objects (tests/stay_table.R), whose checks tie them to the
# run of record; fits nothing and does not touch the MIMIC-IV test set.
#
# The score gives an unmeasured channel zero evidence, so a sparsely measured
# stay has a score pulled toward zero and lands in the middle strata more often.
# A stratum then mixes densely and sparsely measured stays, and its reported
# mortality is an average over both. This script asks, by subgroup and site:
#   1. discrimination: AUROC of the score inside the subgroup, with a patient
#      bootstrap interval (score_metrics(), the canonical owner)
#   2. strata: where the subgroup's stays fall among the quartiles, and the
#      observed mortality of each quartile against
#        - the reference rate (MIMIC-IV training, out of fold, all stays), the
#          number the patient card reports, and
#        - the site's own pooled rate for that quartile, which separates the
#          subgroup effect from the site's transport gap at eICU
#
# SUBGROUPS, fixed 2026-10-27 from the OUTCOME-BLIND distribution of measured
# channels in staytab_20261027T215332, before any outcome was read by subgroup
# (author, on the agent's recommendation):
#   count     measured channels <= 14, 15-16, 17, 18-19
#   gcs       any Glasgow component measured, or none
#   lactate   lactate measured, or not
#   urine     urine output rate measured, or not
# The channel splits separate charting practice from sickness more cleanly than
# the count: at eICU whole channels go uncharted (Glasgow 59% measured against
# 89% at MIMIC-IV, urine output 66% against 96%).
#
# AUROC inside a subgroup also moves with the subgroup's case mix, so a lower
# value can mean a narrower range of severity rather than a weaker score; the
# event rate and stratum spread are printed beside it for that reading.
#
#   Rscript tests/missing_subgroups.R [--staytab out/runs/staytab_...]
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(yaml))
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
STAYTAB <- .opt("--staytab", "out/runs/staytab_20261027T215332")
SITES <- c("mimic_oof", "eicu")
COUNT_BREAKS <- c(-Inf, 14, 16, 17, 19)           # <= 14, 15-16, 17, 18-19
COUNT_LABELS <- c("<=14", "15-16", "17", "18-19")
# Cumulative thresholds below the lowest band (author, 2026-10-29): each is its
# own two-way grouping, "<= t" against "> t", so the sparsely measured tail can
# be read at several depths. Fixed from the outcome-blind count, like the bands.
CUM_THRESHOLDS <- c(5L, 7L, 10L)
GCS <- c("gcs_motor", "gcs_eyes", "gcs_verbal")
MIN_CELL <- 10L                                    # no rate from fewer stays (DUA floor)

cfg <- load_config("config/config.yml")
N_BOOT <- as.integer(cfg_req(cfg, "metrics", "n_boot"))
SEED   <- as.integer(cfg_req(cfg, "seed"))

rd <- function(nm) readRDS(file.path(STAYTAB, "tables", paste0(nm, ".rds")))
chk <- read.csv(file.path(STAYTAB, "diagnostics", "stay_table_checks.csv"), stringsAsFactors = FALSE)
if (any(chk$status != "PASS")) stop("stage A run ", STAYTAB, " has failed checks; rerun tests/stay_table.R", call. = FALSE)

#' The four groupings for one site, as one factor each, aligned to the table.
groupings <- function(tab, meas) {
  stopifnot(nrow(meas) == nrow(tab), all(c(GCS, "lactate", "urine_output_rate") %in% colnames(meas)))
  cum <- lapply(CUM_THRESHOLDS, function(t)
    factor(ifelse(tab$m <= t, paste0("<=", t), paste0(">", t)), levels = c(paste0("<=", t), paste0(">", t))))
  names(cum) <- paste0("count_le", CUM_THRESHOLDS)
  c(list(
    all     = factor(rep("all stays", nrow(tab))),
    count   = factor(COUNT_LABELS[findInterval(tab$m, COUNT_BREAKS, left.open = TRUE)], levels = COUNT_LABELS),
    gcs     = factor(ifelse(rowSums(meas[, GCS, drop = FALSE]) > 0, "measured", "not measured"),
                     levels = c("measured", "not measured")),
    lactate = factor(ifelse(meas[, "lactate"], "measured", "not measured"), levels = c("measured", "not measured")),
    urine   = factor(ifelse(meas[, "urine_output_rate"], "measured", "not measured"), levels = c("measured", "not measured"))),
    cum)
}

quartile_rates <- function(tab) {
  q <- sort(unique(tab$quartile))
  data.frame(quartile = q, n = as.integer(tabulate(tab$quartile, max(q))[q]),
             deaths = as.integer(tapply(tab$y, factor(tab$quartile, levels = q), sum)))
}

ref <- quartile_rates(rd("stay_table_mimic_oof"))
ref$rate <- ref$deaths / ref$n

disc_rows <- list(); strata_rows <- list(); size_rows <- list()
for (site in SITES) {
  tab <- rd(paste0("stay_table_", site)); meas <- rd(paste0("measured_", site))
  pooled <- quartile_rates(tab); pooled$rate <- pooled$deaths / pooled$n
  G <- groupings(tab, meas)
  for (gn in names(G)) {
    g <- G[[gn]]
    if (anyNA(g)) stop("grouping `", gn, "` left stays unassigned at ", site, call. = FALSE)
    for (lv in levels(g)) {
      s <- g == lv
      n <- sum(s); k <- sum(tab$y[s])
      size_rows[[length(size_rows) + 1L]] <- data.frame(site = site, grouping = gn, group = lv,
        n_stays = n, share_of_site = n / nrow(tab), deaths = k, mortality = if (n >= MIN_CELL) k / n else NA_real_,
        median_measured = if (n >= MIN_CELL) stats::median(tab$m[s]) else NA_real_)
      if (n >= MIN_CELL && k > 0 && k < n) {
        sm <- score_metrics(tab$S[s], tab$y[s], label = paste(gn, lv), n_boot = N_BOOT, seed = SEED, group = tab$group[s])
      } else {
        sm <- data.frame(label = paste(gn, lv), n = n, n_events = k, event_rate = NA, auroc = NA, auprc = NA,
                         auprc_lift = NA, auroc_lo = NA, auroc_hi = NA, auprc_lo = NA, auprc_hi = NA,
                         n_boot = 0L, boot_unit = NA, n_units = NA)
      }
      disc_rows[[length(disc_rows) + 1L]] <- cbind(data.frame(site = site, grouping = gn, group = lv), sm[, -1])
      for (qq in ref$quartile) {
        sq <- s & tab$quartile == qq
        nq <- sum(sq); kq <- sum(tab$y[sq])
        ok <- nq >= MIN_CELL
        ci <- if (ok) .wilson(kq, nq) else c(NA_real_, NA_real_)
        obs <- if (ok) kq / nq else NA_real_
        strata_rows[[length(strata_rows) + 1L]] <- data.frame(site = site, grouping = gn, group = lv, quartile = qq,
          n = nq, share_of_group = nq / n, deaths = kq, observed = obs, observed_lo = ci[1], observed_hi = ci[2],
          reference = ref$rate[ref$quartile == qq], site_pooled = pooled$rate[pooled$quartile == qq],
          diff_vs_reference = obs - ref$rate[ref$quartile == qq], diff_vs_site = obs - pooled$rate[pooled$quartile == qq],
          site_inside_ci = if (ok) pooled$rate[pooled$quartile == qq] >= ci[1] & pooled$rate[pooled$quartile == qq] <= ci[2] else NA,
          ci_method = "wilson_stay", below_floor = !ok)
      }
    }
  }
}
sizes <- do.call(rbind, size_rows); disc <- do.call(rbind, disc_rows); strata <- do.call(rbind, strata_rows)
rownames(disc) <- NULL

run <- new_run("missing", cfg, note = "analysis queue 1a: score discrimination and stratum mortality by measurement subgroup; non-fitting")
save_table(run, sizes, "missing_subgroup_sizes")
save_table(run, disc, "missing_subgroup_discrimination")
save_table(run, strata, "missing_subgroup_strata")
finalize_run(run, extra = list(staytab = STAYTAB, n_boot = N_BOOT, seed = SEED, min_cell = MIN_CELL,
                               count_groups = paste(COUNT_LABELS, collapse = ","), cumulative_thresholds = paste(CUM_THRESHOLDS, collapse = ",")))

# --- console: aggregate only ------------------------------------------------
options(width = 200)
r3 <- function(x) { x[] <- lapply(x, function(v) if (is.numeric(v)) round(v, 3) else v); x }
cat("\n== subgroup sizes\n")
print(r3(sizes), row.names = FALSE)
cat("\n== 1. discrimination inside each subgroup (patient bootstrap)\n")
print(r3(disc[, c("site", "grouping", "group", "n", "n_events", "event_rate", "auroc", "auroc_lo", "auroc_hi", "boot_unit")]), row.names = FALSE)
cat("\n== 2. quartile mortality by subgroup (observed vs MIMIC-IV reference and vs the site's pooled rate)\n")
print(r3(strata[strata$grouping != "all", c("site", "grouping", "group", "quartile", "n", "share_of_group", "observed",
                                             "observed_lo", "observed_hi", "reference", "site_pooled", "diff_vs_site", "site_inside_ci")]),
      row.names = FALSE)
cat("\nrun:", run$path, "\n")
