# tests/stay_table.R -----------------------------------------------------------
# STAGE A of the missingness analysis (analysis queue item 1) and the channel
# heterogeneity analysis (docs/heterogeneity_plan_20261027.md). NON-FITTING.
#
# Builds, once, the per-stay table both analyses read, at the two sites:
#   MIMIC-IV training, out of fold   the targets store (l_mats_zero / l_mats_na,
#                                    oof_scores, cutpoints, y_train, group_of_train)
#   eICU                             the external run's saved `l_eicu` and
#                                    `scores_eicu`, the bundle's frozen cut points
# Per stay: outcome, patient group, score S (llr_sum), frozen bin / quartile /
# decile, number of measured channels m, evidence toward death P, toward
# survival M, and total evidence T = P + M. Beside it, the n x 19 paired
# weight-of-evidence matrix (zero when unmeasured) and the measured indicator.
#
# Nothing is analysed against outcome here. The checks that must pass before any
# later stage reads these objects:
#   1. P - M, and the row sum of the matrix, reproduce the stored llr_sum
#   2. the zero-filled and NA-filled matrices agree, so "measured" is well defined
#   3. the frozen bins reproduce the stored per-bin counts of the run of record
#   4. the store's cut points equal the bundle's
#   5. the eICU weights of evidence are identical under the 2026-09-09 and
#      2026-10-25 external runs (the bundles differ only in severity cut points)
#   6. the channel columns match config/domains.csv
#
# Row-level objects are written to the run directory only, as run/external.R
# writes `scores_eicu`; the console gets checks and outcome-blind counts.
#
#   Rscript tests/stay_table.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(yaml); library(targets); library(qs2) })
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
INTERNAL  <- .opt("--internal", "out/runs/internal_20260909T112643")
EXTERNAL  <- .opt("--external", "out/runs/external_20260909T112824")
EXTERNAL2 <- .opt("--external2", "out/runs/external_20261025T211621")
STORE     <- .opt("--store", "_targets")
ARM <- "llr_sum"
MT  <- LLR_ARM_MATRIX[[ARM]]

cfg <- load_config("config/config.yml")
N_BINS <- as.integer(cfg_req(cfg, "metrics", "n_bins"))
if (N_BINS %% 4L != 0L || N_BINS %% 10L != 0L) {
  stop("metrics.n_bins = ", N_BINS, " does not divide into quartiles and deciles", call. = FALSE)
}

checks <- list()
check <- function(site, name, ok, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(site = site, check = name,
    status = if (isTRUE(ok)) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE)
  cat(sprintf("  [%s] %-14s %-44s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", site, name, detail))
}

# The SAME bin rule as risk_bins(breaks = ): sorted unique breaks, ends opened.
assign_bin <- function(score, breaks) {
  br <- unique(sort(as.numeric(breaks))); br[1L] <- -Inf; br[length(br)] <- Inf
  if (length(br) - 1L != N_BINS) stop("frozen breaks give ", length(br) - 1L, " bins, not ", N_BINS, call. = FALSE)
  cut(score, breaks = br, include.lowest = TRUE, labels = FALSE)
}
bin_counts_match <- function(bin, y, run, file) {
  b <- read.csv(file.path(run, "tables", file), stringsAsFactors = FALSE)
  n <- tabulate(bin, N_BINS); d <- tabulate(bin[y == 1L], N_BINS)
  n_ref <- integer(N_BINS); d_ref <- integer(N_BINS)
  n_ref[b$bin] <- b$n; d_ref[b$bin] <- b$deaths
  c(ok = identical(n, n_ref) && identical(d, d_ref),
    max_n = max(abs(n - n_ref)), max_d = max(abs(d - d_ref)))
}

#' One site: the per-stay table and its matrices, with checks 1 and 2.
build_site <- function(site, Lz, Ln, S_ref, y, grp, ids, breaks) {
  if (!is.null(rownames(Lz))) check(site, "matrix rows are the scored stays", identical(rownames(Lz), ids))
  check(site, "zero/NA matrices same shape", identical(dim(Lz), dim(Ln)) && identical(colnames(Lz), colnames(Ln)))
  meas <- !is.na(Ln)
  check(site, "zero fill only where unmeasured",
        all(Lz[!meas] == 0) && isTRUE(all(Lz[meas] == Ln[meas])),
        sprintf("%d measured cells", sum(meas)))
  S <- rowSums(Lz)
  dS <- max(abs(S - S_ref))
  check(site, "row sum reproduces stored llr_sum", dS <= 1e-10, sprintf("max |diff| %.2e", dS))
  P <- rowSums(pmax(Lz, 0)); M <- rowSums(pmax(-Lz, 0))
  dPM <- max(abs((P - M) - S_ref))
  check(site, "P - M reproduces stored llr_sum", dPM <= 1e-10, sprintf("max |diff| %.2e", dPM))
  bin <- assign_bin(S_ref, breaks)
  tab <- data.frame(site = site, stay_id = ids, group = grp, y = as.integer(y), S = as.numeric(S_ref),
                    bin = bin, quartile = (bin - 1L) %/% (N_BINS %/% 4L) + 1L,
                    decile = (bin - 1L) %/% (N_BINS %/% 10L) + 1L,
                    m = as.integer(rowSums(meas)), P = P, M = M, T = P + M,
                    stringsAsFactors = FALSE)
  list(tab = tab, L = Lz, measured = meas)
}

# --- MIMIC-IV, out of fold --------------------------------------------------
cat("\n== MIMIC-IV training, out of fold (store: ", STORE, ")\n", sep = "")
rd <- function(nm) targets::tar_read_raw(nm, store = STORE)
train_ids <- rd("train_ids")
Lz_m <- rd("l_mats_zero")[[MT]]; Ln_m <- rd("l_mats_na")[[MT]]
S_m  <- rd("oof_scores")[[ARM]]
cut_store <- rd("cutpoints")[[ARM]]
mi <- build_site("mimic_oof", Lz_m, Ln_m, S_m, rd("y_train"), rd("group_of_train"),
                 as.character(train_ids), cut_store)
r <- bin_counts_match(mi$tab$bin, mi$tab$y, INTERNAL, sprintf("risk_bins_%s_oof.csv", ARM))
check("mimic_oof", "frozen bins reproduce run-of-record counts", r[["ok"]] == 1,
      sprintf("max |n diff| %d, max |deaths diff| %d (%s)", r[["max_n"]], r[["max_d"]], basename(INTERNAL)))

# --- eICU -------------------------------------------------------------------
cat("\n== eICU (", basename(EXTERNAL), ")\n", sep = "")
xcfg <- yaml::read_yaml("config/external.yml")
bundle <- load_bundle(cfg_req(xcfg, "bundle"), verbose = FALSE)
cfg_b <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
dcut <- max(abs(as.numeric(bundle$cutpoints[[ARM]]) - as.numeric(cut_store)))
check("both", "store cut points equal the bundle's", length(bundle$cutpoints[[ARM]]) == length(cut_store) && dcut == 0,
      sprintf("max |diff| %.2e", dcut))

tabs_e <- load_tables(cfg_b$paths, cfg_b, site = "eicu", verbose = FALSE)
ids_e <- sort(unique(tabs_e$cohort$stay_id))              # as run/external.R
y_e <- as.integer(tabs_e$cohort$mortality[match(as.character(ids_e), as.character(tabs_e$cohort$stay_id))])
if (anyNA(y_e)) stop("stay_table: an eICU stay has no mortality value", call. = FALSE)
grp_e <- patient_group_of(tabs_e$cohort, cfg_b, ids_e)
l_e <- readRDS(file.path(EXTERNAL, "tables", "l_eicu.rds"))
sc_e <- readRDS(file.path(EXTERNAL, "tables", "scores_eicu.rds"))[[ARM]]
check("eicu", "stored scores cover the scored stays",
      !is.null(names(sc_e)) && setequal(names(sc_e), as.character(ids_e)), sprintf("%d stays", length(ids_e)))
S_e <- as.numeric(sc_e[as.character(ids_e)])
Lm_e <- l_matrices(l_e, tabs_e, cfg_b, ids_e, fill = "zero")[[MT]]
Ln_e <- l_matrices(l_e, tabs_e, cfg_b, ids_e, fill = "na")[[MT]]
ei <- build_site("eicu", Lm_e, Ln_e, S_e, y_e, grp_e, as.character(ids_e), bundle$cutpoints[[ARM]])
r <- bin_counts_match(ei$tab$bin, ei$tab$y, EXTERNAL, sprintf("risk_bins_%s_frozen.csv", ARM))
check("eicu", "frozen bins reproduce run-of-record counts", r[["ok"]] == 1,
      sprintf("max |n diff| %d, max |deaths diff| %d (%s)", r[["max_n"]], r[["max_d"]], basename(EXTERNAL)))
l_e2 <- readRDS(file.path(EXTERNAL2, "tables", "l_eicu.rds"))
check("eicu", "weights identical under both bundles", identical(l_e, l_e2),
      sprintf("%s vs %s", basename(EXTERNAL), basename(EXTERNAL2)))
rm(l_e, l_e2, tabs_e)

# --- channel columns ----------------------------------------------------------
dom <- read.csv("config/domains.csv", stringsAsFactors = FALSE)
check("both", "channel columns match domains.csv",
      identical(colnames(mi$L), colnames(ei$L)) && setequal(colnames(mi$L), dom$signal) && ncol(mi$L) == nrow(dom),
      sprintf("%d channels", ncol(mi$L)))

chk <- do.call(rbind, checks)

# --- outcome-blind summaries ------------------------------------------------
m_dist <- do.call(rbind, lapply(list(mi$tab, ei$tab), function(t) {
  n <- tabulate(t$m + 1L, ncol(mi$L) + 1L)
  data.frame(site = t$site[1], m = 0:ncol(mi$L), n_stays = n, share = n / nrow(t), cum_share = cumsum(n) / nrow(t))
}))
ch_share <- data.frame(channel = colnames(mi$L), domain = dom$domain[match(colnames(mi$L), dom$signal)],
                       measured_share_mimic_oof = colMeans(mi$measured), measured_share_eicu = colMeans(ei$measured))

# --- export -------------------------------------------------------------------
run <- new_run("staytab", cfg, note = "stage A: per-stay table for the missingness and heterogeneity analyses; non-fitting")
for (s in list(mi, ei)) {
  nm <- s$tab$site[1]
  save_object(run, s$tab, paste0("stay_table_", nm))
  save_object(run, s$L, paste0("woe_", nm))
  save_object(run, s$measured, paste0("measured_", nm))
}
save_table(run, chk, "stay_table_checks", subdir = "diagnostics")
save_table(run, m_dist, "measured_count_distribution")
save_table(run, ch_share, "channel_measured_share")
finalize_run(run, extra = list(internal = INTERNAL, external = EXTERNAL, external2 = EXTERNAL2,
                               store = STORE, arm = ARM, matrix = MT,
                               checks_failed = sum(chk$status == "FAIL")))

# --- console: aggregate only ------------------------------------------------
options(width = 160)
cat(sprintf("\nstays: MIMIC-IV out of fold %d, eICU %d\n", nrow(mi$tab), nrow(ei$tab)))
cat("\n== measured channels per stay (outcome-blind)\n")
w <- reshape(m_dist[, c("site", "m", "n_stays", "cum_share")], idvar = "m", timevar = "site", direction = "wide")
w[, grepl("cum_share", names(w))] <- round(w[, grepl("cum_share", names(w))], 3)
print(w, row.names = FALSE)
cat("\n== share of stays with each channel measured\n")
ch_share[, 3:4] <- round(ch_share[, 3:4], 3)
print(ch_share, row.names = FALSE)
cat(sprintf("\nchecks: %d PASS, %d FAIL\nrun: %s\n", sum(chk$status == "PASS"), sum(chk$status == "FAIL"), run$path))
if (any(chk$status == "FAIL")) quit(status = 1L)
