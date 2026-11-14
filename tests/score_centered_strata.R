# tests/score_centered_strata.R ---------------------------------------------------
# CONTROL (i) FOR THE FIRST LADDER STEP (author, 2026-10-30). Fits nothing.
#
# QUESTION. The gated joint GAM differs from the score in three ways: joint
# estimation, one intervention function in place of a copy per channel, and no
# channel level. A channel's weight of evidence, averaged over the stays that
# measured it, is negative for every channel (L is calibrated against the
# measured mortality, and the mean of logits lies below the logit of the mean),
# so in the score each unmeasured channel raises a stay by that channel's level.
# The gated joint GAM centres every measurement term on the measured stays and
# has no such level. Does the level explain why the score's strata transport and
# the gated joint GAM's do not?
#
# THE CENTRED SCORE. c_g = mean out-of-fold L_g over the measured training stays,
# estimated once on training and frozen. Centred score = sum_g M_g (L_g - c_g):
# an average measured channel contributes 0, as in the gated joint GAM. Reference
# strata: 20 equal-count bins of the out-of-fold centred training score
# (score_cutpoints), quartiles = 5 bins, deciles = 2, cut points frozen and
# applied to the MIMIC-IV test set and eICU, exactly as for the reported score.
#
# CHECK. The uncentred score built here must reproduce the reported per-bin
# counts of llr_sum (stays and deaths in each of the 20 bins) at all three sites,
# exactly; otherwise nothing is reported.
#
# Inputs: the runs behind the reported strata (strata_20261030T151638). Test
# outcomes are read under the 2026-10-30 rule (descriptive; no design choice).
# Aggregates only (hard rule 1).
#
#   Rscript tests/score_centered_strata.R [--internal <run>] [--test <run>] [--external <run>]
# The options (2026-11-13) default to the runs of record; a fresh chain passes its own runs.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(yaml) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
INTERNAL <- .opt("--internal", "out/runs/internal_20260909T112643")
TEST     <- .opt("--test",     "out/runs/test_20260912T215546")
EXTERNAL <- .opt("--external", "out/runs/external_20260909T112824")

xcfg   <- yaml::read_yaml("config/external.yml")
bundle <- load_bundle(as.character(cfg_req(xcfg, "bundle")), verbose = FALSE)
cfg    <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
sigs   <- as.character(unlist(cfg$signals))
N_BINS <- as.integer(cfg_req(load_config("config/config.yml"), "metrics", "n_bins"))
folds  <- tar_read(folds); tr <- tar_read(train_ids); y_tr <- as.integer(tar_read(y_train))
tabs_m <- tar_read(tabs)
test_ids <- folds$stay_id[folds$split == "test"]
tabs_e <- load_tables(cfg$paths, cfg, site = "eicu", verbose = FALSE)
eicu_ids <- tabs_e$cohort$stay_id
y_of <- function(tabs, ids) { v <- as.integer(tabs$cohort$mortality[match(as.character(ids), as.character(tabs$cohort$stay_id))]); if (anyNA(v)) stop("missing outcomes", call. = FALSE); v }
SITES <- list(
  train = list(L = file.path(INTERNAL, "tables/l_oof.rds"),  role = "oof",   ids = tr,       y = y_tr,                     bins = file.path(INTERNAL, "tables/risk_bins_llr_sum_oof.csv")),
  test  = list(L = file.path(TEST,     "tables/l_test.rds"), role = NULL,    ids = test_ids, y = y_of(tabs_m, test_ids),   bins = file.path(TEST, "tables/risk_bins_llr_sum_frozen.csv")),
  eicu  = list(L = file.path(EXTERNAL, "tables/l_eicu.rds"), role = NULL,    ids = eicu_ids, y = y_of(tabs_e, eicu_ids),   bins = file.path(EXTERNAL, "tables/risk_bins_llr_sum_frozen.csv")))

# --- per site: the stay x channel matrix of the score's L, and the measured mask --------
wide_L <- function(S) {
  L <- readRDS(S$L)
  if (!is.null(S$role)) L <- L[L$role == S$role, ]
  paired <- unique(L$signal[L$model == "full"])
  L <- L[(L$model == "full") | (L$model == "meas" & !(L$signal %in% paired)), ]   # the score's model per channel
  if (!setequal(unique(L$signal), sigs)) stop("channels differ from the design", call. = FALSE)
  if (anyDuplicated(paste(L$signal, L$stay_id))) stop("a stay has two rows for one channel", call. = FALSE)
  i <- match(as.character(L$stay_id), as.character(S$ids))
  if (anyNA(i)) stop("L holds stays outside the site's cohort", call. = FALSE)
  W <- matrix(0, length(S$ids), length(sigs), dimnames = list(NULL, sigs)); M <- W > 0
  j <- match(L$signal, sigs)
  W[cbind(i, j)] <- L$l; M[cbind(i, j)] <- TRUE
  list(W = W, M = M)
}
D <- lapply(SITES, wide_L)

# --- check: the uncentred score reproduces the reported bins -------------------------
score <- lapply(D, function(d) rowSums(d$W))
br <- score_cutpoints(score$train, N_BINS)
for (st in names(SITES)) {
  rb  <- risk_bins(score[[st]], SITES[[st]]$y, N_BINS, breaks = br)
  ref <- utils::read.csv(SITES[[st]]$bins, stringsAsFactors = FALSE)
  if (!identical(as.integer(rb$n), as.integer(ref$n)) || !identical(as.integer(rb$deaths), as.integer(ref$deaths))) {
    stop(st, ": the score built here does not reproduce the reported bins", call. = FALSE)
  }
}
cat("check passed: the uncentred score reproduces the reported bins at all three sites\n")

# --- the channel levels and the centred score ------------------------------------------
c_g <- vapply(sigs, function(g) mean(D$train$W[D$train$M[, g], g]), 0)
centred <- lapply(D, function(d) rowSums(d$W) - as.numeric(d$M %*% c_g))

strata <- function(s, br, size) (cut(s, br, include.lowest = TRUE, labels = FALSE) - 1L) %/% size + 1L
ROWS <- list(); MOVE <- list()
for (arm in c("score", "centred")) {
  s <- if (arm == "score") score else centred
  b <- score_cutpoints(s$train, N_BINS); b[1] <- -Inf; b[length(b)] <- Inf
  for (part in c("quartile", "decile")) {
    size <- if (part == "quartile") 5L else 2L
    ref <- tapply(SITES$train$y, strata(s$train, b, size), mean)
    for (st in c("test", "eicu")) {
      k <- strata(s[[st]], b, size); y <- SITES[[st]]$y
      obs <- tapply(y, factor(k, seq_along(ref)), mean)
      ROWS[[length(ROWS) + 1L]] <- data.frame(arm = arm, site = st, partition = part, stratum = seq_along(ref),
        share = as.numeric(tabulate(k, length(ref)) / length(k)), observed = as.numeric(obs), reference = as.numeric(ref),
        difference = as.numeric(obs - ref), auroc_score = .auroc(s[[st]], y))
      if (part == "quartile") MOVE[[paste(arm, st)]] <- k
    }
  }
}
res <- do.call(rbind, ROWS)
summ <- do.call(rbind, lapply(split(res, paste(res$arm, res$site, res$partition)), function(d) data.frame(
  arm = d$arm[1], site = d$site[1], partition = d$partition[1], max_abs_difference = max(abs(d$difference)),
  top_share = d$share[nrow(d)], top_observed = d$observed[nrow(d)], top_reference = d$reference[nrow(d)],
  auroc = d$auroc_score[1])))
moved <- data.frame(site = c("test", "eicu"),
  same_quartile = c(mean(MOVE[["score test"]] == MOVE[["centred test"]]), mean(MOVE[["score eicu"]] == MOVE[["centred eicu"]])),
  mean_level_shift = c(mean(as.numeric(D$test$M %*% c_g)), mean(as.numeric(D$eicu$M %*% c_g))),
  mean_measured = c(mean(rowSums(D$test$M)), mean(rowSums(D$eicu$M))))
levels_tab <- data.frame(signal = sigs, level = unname(c_g), n_measured_train = colSums(D$train$M))

run <- new_run("centstrata", load_config("config/config.yml"), note = "control (i) for the first ladder step: the score with each channel's measured-stay level removed; fits nothing")
save_table(run, levels_tab, "channel_levels"); save_table(run, res, "strata_by_stratum")
save_table(run, summ, "strata_summary"); save_table(run, moved, "quartile_movement")
finalize_run(run, extra = list(internal = INTERNAL, test = TEST, external = EXTERNAL))
options(width = 200)
r2 <- function(d) { d[] <- lapply(d, function(v) if (is.numeric(v)) round(v, 3) else v); d }
cat("\n== channel levels (mean out-of-fold L over measured training stays, nats)\n"); print(r2(levels_tab), row.names = FALSE)
cat(sprintf("sum of levels %.2f nats\n", sum(c_g)))
cat("\n== strata, score against centred score\n"); print(r2(summ[order(summ$partition, summ$site, summ$arm), ]), row.names = FALSE)
cat("\n== quartile assignment, score against centred score\n"); print(r2(moved), row.names = FALSE)
cat("\nrun:", run$path, "\n")
