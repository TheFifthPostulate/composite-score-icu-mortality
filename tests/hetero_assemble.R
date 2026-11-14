# tests/hetero_assemble.R ------------------------------------------------------
# STAGE 1 of the channel heterogeneity analysis
# (docs/heterogeneity_plan_20261027.md). NON-FITTING.
#
# Reads the stage A per-stay objects (tests/stay_table.R) and adds, at both sites:
#   - the frozen extreme cut-offs: for each channel, the 5th, 10th, 80th and 90th
#     percentiles of its out-of-fold weights of evidence among MEASURED MIMIC-IV
#     training stays, computed once without reference to outcome and applied
#     unchanged at eICU
#   - per stay, the number of channels extreme toward death at each high level
#     (L at or above the cut-off AND above zero) and toward survival at each low
#     level (at or below the cut-off AND below zero), by channel and by domain
#     (config/domains.csv; a domain counts once)
#   - the strongest death-ward weight of evidence and its channel
#   - flags for channels unmeasured BY DESIGN: every Glasgow hour removed by the
#     sedation/paralytic rule, or urine output / creatinine removed by the
#     dialysis rule (`dx_n_masked`, as in tests/missing_masking_audit.R)
#
# Prints only the cut-offs (quantiles over thousands of stays) and outcome-blind
# shares of measured stays flagged extreme, so a cut-off that the sign condition
# makes unreachable is visible before any outcome is read.
#
#   Rscript tests/hetero_assemble.R [--staytab out/runs/staytab_...]
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(yaml))
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)
source("tests/hetero_common.R")

args <- commandArgs(trailingOnly = TRUE)
STAYTAB <- .hopt(args, "--staytab", "out/runs/staytab_20261027T215332")
GCS <- c("gcs_motor", "gcs_eyes", "gcs_verbal")
RENAL <- c("urine_output_rate", "creatinine")

st_chk <- read.csv(file.path(STAYTAB, "diagnostics", "stay_table_checks.csv"), stringsAsFactors = FALSE)
if (any(st_chk$status != "PASS")) stop("stage A run ", STAYTAB, " has failed checks", call. = FALSE)
rd <- function(nm) readRDS(file.path(STAYTAB, "tables", paste0(nm, ".rds")))

checks <- list()
check <- function(name, ok, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(check = name, status = if (isTRUE(ok)) "PASS" else "FAIL",
                                               detail = detail, stringsAsFactors = FALSE)
  cat(sprintf("  [%s] %-52s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, detail))
}

# --- cut-offs, MIMIC-IV out of fold, measured stays only ---------------------
W_m <- rd("woe_mimic_oof"); M_m <- rd("measured_mimic_oof")
ch <- colnames(W_m)
probs <- c(HETERO_LOW, HETERO_HIGH) / 100
cut <- t(vapply(ch, function(g) stats::quantile(W_m[M_m[, g], g], probs, names = FALSE), numeric(length(probs))))
colnames(cut) <- paste0("p", c(HETERO_LOW, HETERO_HIGH))
cutoffs <- data.frame(channel = ch, n_measured = colSums(M_m), cut, stringsAsFactors = FALSE, row.names = NULL)
for (h in HETERO_HIGH) cutoffs[[paste0("reachable_high", h)]] <- cutoffs[[paste0("p", h)]] > 0
for (l in HETERO_LOW)  cutoffs[[paste0("reachable_low", l)]]  <- cutoffs[[paste0("p", l)]] < 0

dom <- read.csv("config/domains.csv", stringsAsFactors = FALSE)
check("channels match domains.csv", setequal(ch, dom$signal) && length(ch) == nrow(dom), sprintf("%d channels", length(ch)))
domains <- split(dom$signal, dom$domain)

# --- channels removed by design ------------------------------------------------
cfg <- load_config("config/config.yml")
xcfg <- yaml::read_yaml("config/external.yml")
bundle <- load_bundle(cfg_req(xcfg, "bundle"), verbose = FALSE)
cfg_e <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
qbreaks <- {
  br <- unique(sort(as.numeric(bundle$cutpoints[["llr_sum"]])))
  nb <- as.integer(cfg_req(cfg, "metrics", "n_bins"))
  br[1L + (nb %/% 4L) * 1:3]                     # the boundaries between quartiles
}
sf_of <- list(mimic_oof = function() load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)$signal_features,
              eicu      = function() load_tables(cfg_e$paths, cfg_e, site = "eicu", verbose = FALSE)$signal_features)
removed_by_design <- function(sf, ids, meas, signals) {
  hit <- sf[sf$signal %in% signals & sf$dx_n_masked > 0, c("stay_id", "signal")]
  vapply(signals, function(g) !meas[, g] & as.character(ids) %in% as.character(hit$stay_id[hit$signal == g]),
         logical(length(ids)))
}

# --- per site ---------------------------------------------------------------
build <- function(site) {
  tab <- rd(paste0("stay_table_", site)); W <- rd(paste0("woe_", site)); M <- rd(paste0("measured_", site))
  stopifnot(identical(colnames(W), ch), identical(colnames(M), ch), nrow(W) == nrow(tab))
  ext <- list()
  for (h in HETERO_HIGH) ext[[paste0("high", h)]] <- M & sweep(W, 2, cut[, paste0("p", h)], ">=") & W > 0
  for (l in HETERO_LOW)  ext[[paste0("low", l)]]  <- M & sweep(W, 2, cut[, paste0("p", l)], "<=") & W < 0
  for (nm in names(ext)) {
    E <- ext[[nm]]
    tab[[paste0("n_", nm)]] <- as.integer(rowSums(E))
    tab[[paste0("ndom_", nm)]] <- as.integer(rowSums(vapply(domains, function(d) rowSums(E[, d, drop = FALSE]) > 0,
                                                             logical(nrow(E)))))
  }
  mx <- apply(W, 1, max)
  tab$strongest_deathward <- pmax(mx, 0)
  tab$strongest_channel <- ifelse(mx > 0, ch[max.col(W, ties.method = "first")], NA_character_)
  sf <- sf_of[[site]]()
  if (!"dx_n_masked" %in% names(sf)) stop("signal_features has no dx_n_masked at ", site, call. = FALSE)
  R_g <- removed_by_design(sf, tab$stay_id, M, GCS)
  R_r <- removed_by_design(sf, tab$stay_id, M, RENAL)
  tab$glasgow_removed <- rowSums(M[, GCS, drop = FALSE]) == 0 & rowSums(R_g) > 0
  tab$renal_removed <- rowSums(R_r) > 0
  tab$any_removed <- tab$glasgow_removed | tab$renal_removed
  check(sprintf("%s: P - M still equals S", site), max(abs(tab$P - tab$M - tab$S)) <= 1e-10)
  check(sprintf("%s: extremes only on measured channels", site),
        all(vapply(ext, function(E) !any(E & !M), logical(1))))
  check(sprintf("%s: sign condition holds", site),
        all(W[ext$high80 | ext$high90] > 0) && all(W[ext$low5 | ext$low10] < 0))
  list(tab = tab, ext = ext, M = M)
}
S <- list(mimic_oof = build("mimic_oof"), eicu = build("eicu"))
check("90th-level extremes are a subset of 80th", all(vapply(S, function(s) !any(s$ext$high90 & !s$ext$high80), logical(1))))
check("5th-level extremes are a subset of 10th", all(vapply(S, function(s) !any(s$ext$low5 & !s$ext$low10), logical(1))))

# Outcome-blind: share of MEASURED stays flagged per channel and level.
flag_share <- do.call(rbind, lapply(names(S), function(site) {
  s <- S[[site]]
  data.frame(site = site, channel = ch, n_measured = colSums(s$M),
             vapply(s$ext, function(E) colSums(E) / pmax(colSums(s$M), 1), numeric(length(ch))),
             row.names = NULL, check.names = FALSE)
}))
chk <- do.call(rbind, checks)

# --- export -------------------------------------------------------------------
run <- new_run("heteroasm", cfg, note = "heterogeneity stage 1: frozen extreme cut-offs, extreme counts, removed-by-design flags; non-fitting")
for (site in names(S)) {
  save_object(run, S[[site]]$tab, paste0("hetero_stay_", site))
  save_object(run, S[[site]]$ext, paste0("hetero_extreme_", site))
}
save_object(run, qbreaks, "quartile_breaks")
save_table(run, cutoffs, "hetero_cutoffs")
save_table(run, flag_share, "hetero_flag_share")
save_table(run, chk, "hetero_assemble_checks", subdir = "diagnostics")
finalize_run(run, extra = list(staytab = STAYTAB, high = paste(HETERO_HIGH, collapse = ","),
                               low = paste(HETERO_LOW, collapse = ","), checks_failed = sum(chk$status == "FAIL")))

# --- console: aggregate only ------------------------------------------------
options(width = 200)
r3 <- function(x) { x[] <- lapply(x, function(v) if (is.numeric(v)) round(v, 3) else v); x }
cat("\n== frozen cut-offs (weights of evidence, nats; measured MIMIC-IV training stays)\n")
print(r3(cutoffs), row.names = FALSE)
cat("\n== share of measured stays flagged extreme (outcome-blind)\n")
print(r3(flag_share), row.names = FALSE)
for (site in names(S)) {
  t <- S[[site]]$tab
  cat(sprintf("\n%s: Glasgow removed by design %d stays, renal removed by design %d, either %d\n",
              site, sum(t$glasgow_removed), sum(t$renal_removed), sum(t$any_removed)))
}
cat(sprintf("\nchecks: %d PASS, %d FAIL\nrun: %s\n", sum(chk$status == "PASS"), sum(chk$status == "FAIL"), run$path))
if (any(chk$status == "FAIL")) quit(status = 1L)
