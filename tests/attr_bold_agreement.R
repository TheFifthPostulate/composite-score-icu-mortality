# tests/attr_bold_agreement.R ---------------------------------------------------
# DO THE CARD'S BOLD SYSTEMS COME BACK AFTER RETRAINING?
# (docs/bold_reproducibility_plan_20261028.md)
#
# A non-fitting consumer of the attribution replicate store, like
# tests/attr_metrics.R, and guarded the same way: the live design key and the
# store fingerprint are recomputed and compared before any replicate is read.
#
# THE RULE, AS THE CARD DRAWS IT (paper/make_card_figure.R). A measured channel
# is bold when its percentile among measured stays, round(100 * mean(ref <= x)),
# is 85 or above. A system is bold when any of its channels is
# (config/domains.csv). Unmeasured channels are never bold, for either method.
#   reference "own":    each replicate's own measured values, as a retrained
#                       model would recompute the card's reference (primary)
#   reference "frozen": the card's cut-offs, the method's unperturbed replicate
#                       (boot 0, draw 0, seed 0) of the INTERNAL store (MIMIC-IV
#                       training, out of fold), held fixed for every replicate
#                       (sensitivity). At eICU the card judges every stay
#                       against these MIMIC-IV cut-offs.
#
# PER PAIR OF REPLICATES, OVER PATIENTS:
#   jaccard    mean Jaccard index of the two bold-system sets, over patients with
#              at least one bold system in either replicate
#   retention  of all bold system marks in one replicate, the share also bold in
#              the other, averaged over the two directions
#   identical  share of patients whose two sets are identical, both empty included
#   *_channel  the same at the channel level, a control: every channel bolds
#              the same share of its measured stays under both methods, so a
#              method whose bold channels cluster within one system cannot
#              gain agreement from the grouping here
#   top1       the existing leading-channel agreement on the same pair
#              (attr_topk_agree, signal level, k = 1), reported beside them.
#              NOT MASKED: it ranks |contribution| over all 19 columns. An
#              unmeasured channel is exactly 0 under LLR but carries SHAP's
#              attribution to the missing value, so SHAP's leading channel can
#              be an unmeasured one (found 2026-10-30). Kept as it was for
#              continuity with tests/attr_metrics.R, and reported beside:
#   top1_measured  the same with unmeasured channels set to 0 for both
#              methods, the masking the bold rule uses; and, per method, the
#              share of patients whose unmasked leading channel is unmeasured
#              (unperturbed replicate).
#
# BY MISSINGNESS (author, 2026-10-30): jaccard and identical within patients
# grouped by how many channels they have unmeasured (0, 1-2, 3 or more). The
# bold rule already masks unmeasured channels; this asks whether SHAP's lower
# agreement concentrates in sparse stays, where its booster leans on
# missingness, or holds in fully measured ones.
# Each metric is then a distribution over pairs. The pairs and their contrasts
# come from attr_contrast_pairs(); the population is the common bag population,
# bounded and seeded as in tests/attr_metrics.R.
#
# NEAR THE EDGE (expectation 2 of the plan): the share of identical sets by how
# close a patient's nearest measured channel sits to the 85th-percentile line in
# the pair's first replicate.
#
# Reads the store only. Aggregates only (hard rule 1).
#
#   Rscript tests/attr_bold_agreement.R [out/runs/attrgen_<id>]
#   Rscript tests/attr_bold_agreement.R --site eicu [out/runs/attrextgen_<id>]
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(arrow); library(yaml) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args  <- commandArgs(trailingOnly = TRUE)
# --site eicu: consume the eICU whole-bag store (tests/attr_external_bags.R),
# guarded as tests/attr_metrics.R --site eicu. The own reference is each
# replicate's eICU values; the frozen reference is the card's, the MIMIC-IV
# training out-of-fold distribution, read from the internal store's
# unperturbed replicate (the plan: "the frozen cut-offs of the card").
.si   <- which(args == "--site")
SITE  <- if (length(.si) && .si[1] < length(args)) args[.si[1] + 1L] else "mimic"
if (!SITE %in% c("mimic", "eicu")) stop("--site must be mimic or eicu", call. = FALSE)
pos   <- args[!grepl("^--", args) & !(seq_along(args) %in% (.si + 1L))]
# The default store is the latest complete generator run, as tests/attr_metrics.R.
gen_d <- if (length(pos)) pos[1] else latest_run(if (SITE == "mimic") "attrgen" else "attrextgen", require_complete = TRUE)
ALLOW_PARTIAL <- "--allow-partial" %in% args
ecfg  <- yaml::read_yaml("config/attribution_eval.yml")

METHODS  <- c("llr_full", "shap_xgb_feat")   # the card's channels, and SHAP grouped to channels
PCT_BOLD <- 85                               # paper/make_card_figure.R
EDGE_BANDS <- c(0, 2, 5, Inf)                # percentile points from the 85th line
EDGE_LAB   <- c("<2", "2-5", ">=5")
MISS_BANDS <- c(0, 1, 3, Inf)                # number of unmeasured channels per patient
MISS_LAB   <- c("0", "1-2", ">=3")

# --- the store and its guards (as tests/attr_metrics.R) ------------------------
.gm <- read_manifest(gen_d)
if (!identical(.gm$status, "complete")) stop("the generator run ", basename(gen_d), " is not marked complete", call. = FALSE)
man    <- utils::read.csv(file.path(gen_d, "manifest_replicates.csv"), stringsAsFactors = FALSE)
meas   <- qs2::qs_read(file.path(gen_d, "measured.qs2"))
sub    <- cfg_req(ecfg, "storage", "subdir")
.dg    <- qs2::qs_read(file.path(gen_d, "design.qs2"))

if (SITE == "eicu") {
  source("tests/attr_external_common.R")
  if (!identical(.dg$site, "eicu") || !identical(.dg$generator, EICU_GEN_VERSION)) {
    stop(basename(gen_d), " is not an eICU whole-bag store (", EICU_GEN_VERSION, ")", call. = FALSE)
  }
  bundle <- load_bundle(.dg$bundle, verbose = FALSE)
  cfg    <- bundle_cfg(bundle, paths = cfg_req(yaml::read_yaml("config/external.yml"), "paths"))
  ids_ch <- as.character(.dg$stay_id)
  sigs   <- as.character(unlist(cfg$signals))
  .tabs_e <- load_tables(cfg$paths, cfg, site = "eicu", verbose = FALSE)
  if (!identical(as.character(.tabs_e$cohort$stay_id), ids_ch)) stop("the eICU cohort on disk is not the store's row set", call. = FALSE)
  if (!identical(measured_matrix(.tabs_e, cfg, .dg$stay_id), meas)) stop("the store's measured mask differs from the eICU tables on disk", call. = FALSE)
  rm(.tabs_e)
  INT_D  <- file.path(dirname(gen_d), .dg$internal_store)
  .int_dg <- qs2::qs_read(file.path(INT_D, "design.qs2"))
  .id <- eicu_store_identity(.int_dg$design_key, bundle, cfg, ecfg, .dg$stay_id, .dg$methods)
  if (!identical(.id$design_key, .dg$design_key)) stop("the eICU store's design key differs from the bundle and internal store today", call. = FALSE)
  if (nrow(attr_fingerprint_diff(.dg$fingerprint, .id$fingerprint))) stop("the eICU store's fingerprint differs from the live design", call. = FALSE)
  .want <- .dg$design_key; .fp <- .id$fingerprint
  .plan <- eicu_store_plan(ecfg, n_folds = as.integer(cfg_req(cfg, "n_folds")), .dg$methods)
} else {
cfg    <- tar_read(cfg)
folds  <- tar_read(folds)
tr     <- tar_read(train_ids)
priors <- tar_read(priors)
ids_ch <- as.character(tr)
sigs   <- as.character(unlist(cfg$signals))
INT_D  <- gen_d
if (is.null(.dg$k_ti) || is.null(.dg$fingerprint)) stop("the store carries no k_ti or no fingerprint", call. = FALSE)
.fold_vec <- folds$fold[match(tr, folds$stay_id)]
.want <- attr_design_key(cfg, .fold_vec, .dg$k_ti)
if (nrow(attr_design_diff(.dg$design_key, .want))) stop("replicate store design key differs from the live design", call. = FALSE)
stopifnot(identical(as.character(.dg$stay_id), ids_ch))
.fp <- attr_design_fingerprint(cfg, .fold_vec, .dg$k_ti, priors, ecfg, design_key = .want)
.fd <- attr_fingerprint_diff(.dg$fingerprint, .fp)
if (nrow(.fd)) {
  # Route-selective, as tests/attr_metrics.R: refuse only if a replicate this
  # script reads is on a route the differing field(s) invalidate.
  .st <- attr_stale_routes(.fd)
  .stale <- if (.st$all) rep(TRUE, nrow(man)) else attr_stale_mask(man, .st$rules)
  if (any(.stale & man$method %in% METHODS)) {
    print(.fd[, c("field", "routes_affected")], row.names = FALSE)
    stop("the store fingerprint differs from the live design on routes this script reads", call. = FALSE)
  }
  cat("  fingerprint differs only on routes this script does not read; proceeding.\n")
}
.plan <- attr_replicate_plan(ecfg, n_folds = as.integer(cfg_req(cfg, "n_folds")), n_pass_fits = 0L)
}
cat(sprintf("%s store: design key %s and fingerprint %s match the live design\n", SITE, attr_key_hash(.want), attr_key_hash(.fp)))
attr_validate_manifest(man, "tests/attr_bold_agreement.R")
if (!identical(rownames(meas), ids_ch) || !identical(colnames(meas), sigs)) stop("measured mask is not in store order", call. = FALSE)

.tomb_p <- file.path(gen_d, "llr_bootstrap_excluded.csv")
TOMB <- if (file.exists(.tomb_p)) {
  .t <- utils::read.csv(.tomb_p, stringsAsFactors = FALSE)
  sort(unique(as.integer(.t$boot_id[.t$design %in% c(attr_key_hash(.dg$design_key), .dg$migrated_from)])))
} else integer(0)
COV <- attr_replicate_coverage(.plan, man, tombstoned = TOMB)
COV <- COV[COV$method %in% METHODS, ]
if (any(COV$missing > 0L) && !ALLOW_PARTIAL) stop("the store is partial for the methods read here", call. = FALSE)
STORE_STATUS <- if (any(COV$missing > 0L)) "partial" else "complete"

man <- man[man$method %in% METHODS, ]
boot_man <- man[man$route == "bootstrap", ]
COMMON_BAGS <- Reduce(intersect, lapply(METHODS, function(m) sort(unique(boot_man$boot_id[boot_man$method == m]))))
in_common <- function(b) b == 0L | b %in% COMMON_BAGS
cat(sprintf("common bag population: %d bags; tombstoned LLR bags: %s\n", length(COMMON_BAGS), paste(TOMB, collapse = ", ")))

# --- per replicate: bold-system sets as bitmasks --------------------------------
DOM  <- utils::read.csv("config/domains.csv", stringsAsFactors = FALSE)
DOMS <- unique(DOM$domain)
dom_of <- DOM$domain[match(sigs, DOM$signal)]
if (anyNA(dom_of) || length(DOMS) > 30L) stop("domains.csv does not cover the signals", call. = FALSE)
BIT  <- 2L^(match(dom_of, DOMS) - 1L)          # one bit per system, per channel
CBIT <- 2L^(seq_along(sigs) - 1L)              # one bit per channel (the channel-level control)
if (length(sigs) > 30L) stop("more channels than an integer bitmask holds", call. = FALSE)
POP10 <- vapply(0:1023, function(v) sum(bitwAnd(v, 2L^(0:9)) > 0), integer(1))
popcount <- function(v) POP10[v %% 1024L + 1L] + POP10[(v %/% 1024L) %% 1024L + 1L] + POP10[v %/% 1048576L + 1L]
stopifnot(popcount(c(0, 1, 1023, 1024, 2^19 - 1, 2^30 - 1, 2^18 + 2^10 + 1)) == c(0, 1, 10, 1, 19, 30, 3))

#' Percentile of each measured value against a reference, as the card computes
#' it: 100 * mean(ref <= x), with ref the sorted measured values of the channel.
pct_against <- function(x, ref_sorted) 100 * findInterval(x, ref_sorted) / length(ref_sorted)

encode <- function(M, ref) {
  # A patient with no measured channel has no edge: set far from the line.
  mask <- integer(nrow(M)); cmask <- integer(nrow(M)); edge <- rep(1e9, nrow(M))
  for (j in seq_along(sigs)) {
    w <- which(meas[, j]); if (!length(w)) next
    p <- pct_against(M[w, j], ref[[j]])
    b <- round(p) >= PCT_BOLD
    mask[w[b]]  <- bitwOr(mask[w[b]], BIT[j])
    cmask[w[b]] <- bitwOr(cmask[w[b]], CBIT[j])
    edge[w] <- pmin(edge[w], abs(p - (PCT_BOLD - 0.5)))
  }
  list(mask = mask, cmask = cmask, edge = edge)
}
own_ref <- function(M) lapply(seq_along(sigs), function(j) sort(M[meas[, j], j]))

load_rep <- function(key) {
  M <- qs2::qs_read(file.path(gen_d, sub, paste0(key, ".qs2")))
  attr_check_replicate(M, ids_ch, sigs, key)
  M
}
# THE FROZEN REFERENCE is the card's: the unperturbed replicate of the INTERNAL
# store (MIMIC-IV training, out of fold), on that store's own measured mask.
# At MIMIC-IV this is the store being read; at eICU it is the internal store
# the eICU store mirrors.
int_man  <- utils::read.csv(file.path(INT_D, "manifest_replicates.csv"), stringsAsFactors = FALSE)
int_meas <- qs2::qs_read(file.path(INT_D, "measured.qs2"))
int_ids  <- as.character(qs2::qs_read(file.path(INT_D, "design.qs2"))$stay_id)
if (!identical(colnames(int_meas), sigs) || !identical(rownames(int_meas), int_ids)) stop("internal store mask is not in store order", call. = FALSE)
anchor_ref <- function(m) {
  z <- int_man[int_man$method == m & int_man$boot_id == 0L & int_man$seed_id == 0L & int_man$draw_id == 0L, ]
  if (nrow(z) != 1L) stop(nrow(z), " unperturbed replicates for ", m, " in the internal store; expected exactly one", call. = FALSE)
  M <- qs2::qs_read(file.path(INT_D, sub, paste0(z$key, ".qs2")))
  attr_check_replicate(M, int_ids, sigs, z$key)
  lapply(seq_along(sigs), function(j) sort(M[int_meas[, j], j]))
}
FROZEN <- lapply(stats::setNames(METHODS, METHODS), anchor_ref)

MISS_BAND <- cut(rowSums(!meas), MISS_BANDS, MISS_LAB, right = FALSE)
if (anyNA(MISS_BAND)) stop("a patient falls outside the missingness bands", call. = FALSE)
ENC <- list(); LEAD <- list()
for (i in seq_len(nrow(man))) {
  if (!in_common(man$boot_id[i])) next
  M <- load_rep(man$key[i]); m <- man$method[i]
  own <- encode(M, own_ref(M))
  Mm <- M; Mm[!meas] <- 0                      # unmeasured channels set to 0, as the bold rule masks them
  ENC[[man$key[i]]] <- list(own = own$mask, chan = own$cmask, edge = own$edge, frozen = encode(M, FROZEN[[m]])$mask,
                            top = attr_topk(M, 1L), top_meas = attr_topk(Mm, 1L))
  if (man$boot_id[i] == 0L && man$seed_id[i] == 0L && man$draw_id[i] == 0L) {
    lead <- max.col(abs(M), ties.method = "first")
    LEAD[[m]] <- data.frame(method = m, share_leading_channel_unmeasured = mean(!meas[cbind(seq_len(nrow(M)), lead)]),
                            share_patients_any_unmeasured = mean(rowSums(!meas) > 0), stringsAsFactors = FALSE)
  }
}
lead_s <- do.call(rbind, LEAD)
rm(M); gc(verbose = FALSE)
cat(sprintf("encoded %d replicates\n", length(ENC)))

# --- pair metrics ---------------------------------------------------------------
pair_metrics <- function(a, b) {
  inter <- popcount(bitwAnd(a, b)); uni <- popcount(bitwOr(a, b))
  na <- popcount(a); nb <- popcount(b); any_bold <- uni > 0
  c(jaccard = mean(inter[any_bold] / uni[any_bold]),
    retention = mean(c(sum(inter) / sum(na), sum(inter) / sum(nb))),
    identical = mean(a == b),
    bold_per_patient = mean(c(na, nb)))
}

D_MAXP <- as.integer(cfg_req(ecfg, "distribution", "max_pairs"))
D_SEED <- as.integer(cfg_req(ecfg, "distribution", "pair_seed"))
STRAT  <- as.character(cfg_req(ecfg, "selection", "stability_stratum"))
CP <- attr_contrast_pairs(man)
CP$common <- in_common(CP$boot_i) & in_common(CP$boot_j)
CP <- CP[CP$common, ]
cells <- unique(CP[, c("method", "route_family", "code", "held")])

RAW <- list(); EDGE <- list(); MISSR <- list()
for (ci in seq_len(nrow(cells))) {
  z <- CP[CP$method == cells$method[ci] & CP$route_family == cells$route_family[ci] &
          CP$code == cells$code[ci] & CP$held == cells$held[ci], ]
  if (nrow(z) > D_MAXP) z <- z[with_seed(D_SEED, sort(sample.int(nrow(z), D_MAXP))), ]
  for (q in seq_len(nrow(z))) {
    A <- ENC[[man$key[z$i[q]]]]; B <- ENC[[man$key[z$j[q]]]]
    own <- pair_metrics(A$own, B$own); frz <- pair_metrics(A$frozen, B$frozen); chn <- pair_metrics(A$chan, B$chan)
    RAW[[length(RAW) + 1L]] <- data.frame(method = cells$method[ci], route = cells$route_family[ci],
      contrast = cells$code[ci], held = cells$held[ci], boot_i = z$boot_i[q], boot_j = z$boot_j[q],
      jaccard = own[["jaccard"]], retention = own[["retention"]], identical = own[["identical"]],
      bold_per_patient = own[["bold_per_patient"]],
      jaccard_frozen = frz[["jaccard"]], retention_frozen = frz[["retention"]], identical_frozen = frz[["identical"]],
      jaccard_channel = chn[["jaccard"]], identical_channel = chn[["identical"]], bold_channels_per_patient = chn[["bold_per_patient"]],
      top1_agreement = attr_topk_agree(A$top, B$top, 1L),
      top1_measured_agreement = attr_topk_agree(A$top_meas, B$top_meas, 1L), stringsAsFactors = FALSE)
    inter <- popcount(bitwAnd(A$own, B$own)); uni <- popcount(bitwOr(A$own, B$own))
    MISSR[[length(MISSR) + 1L]] <- data.frame(method = cells$method[ci], contrast = cells$code[ci], held = cells$held[ci],
      band = MISS_LAB, share_patients = as.numeric(table(MISS_BAND)) / length(MISS_BAND),
      jaccard = vapply(MISS_LAB, function(l) { w <- MISS_BAND == l & uni > 0; if (any(w)) mean(inter[w] / uni[w]) else NA_real_ }, numeric(1)),
      identical = vapply(MISS_LAB, function(l) mean(A$own[MISS_BAND == l] == B$own[MISS_BAND == l]), numeric(1)),
      stringsAsFactors = FALSE)
    band <- cut(A$edge, EDGE_BANDS, EDGE_LAB, right = FALSE)
    same <- A$own == B$own
    EDGE[[length(EDGE) + 1L]] <- data.frame(method = cells$method[ci], contrast = cells$code[ci], held = cells$held[ci],
      band = EDGE_LAB, share_patients = as.numeric(table(band)) / length(band),
      identical = vapply(EDGE_LAB, function(l) mean(same[band == l]), numeric(1)), stringsAsFactors = FALSE)
  }
}
raw <- do.call(rbind, RAW); edg <- do.call(rbind, EDGE); mis <- do.call(rbind, MISSR)

MET <- c("jaccard", "retention", "identical", "bold_per_patient", "jaccard_frozen", "retention_frozen",
         "identical_frozen", "jaccard_channel", "identical_channel", "bold_channels_per_patient", "top1_agreement",
         "top1_measured_agreement")
summ <- do.call(rbind, lapply(split(raw, paste(raw$method, raw$route, raw$contrast, raw$held)), function(d) {
  do.call(rbind, lapply(MET, function(mt) {
    v <- d[[mt]]
    data.frame(method = d$method[1], route = d$route[1], contrast = d$contrast[1], held = d$held[1],
               primary_stratum = d$held[1] == STRAT, metric = mt, n_pairs = length(v),
               median = stats::median(v), p025 = unname(stats::quantile(v, 0.025)),
               p975 = unname(stats::quantile(v, 0.975)), min = min(v), max = max(v), stringsAsFactors = FALSE)
  }))
}))
# THE PRIMARY COMPARISON, PAIRED. The bags are shared, so the two methods' L3
# pairs in the stability stratum are the same bag pairs: the difference is taken
# pair by pair, with the share of pairs where LLR agrees more.
pk <- function(d) paste(pmin(d$boot_i, d$boot_j), pmax(d$boot_i, d$boot_j))
la <- raw[raw$method == METHODS[1] & raw$contrast == "L3" & raw$held == STRAT, ]
sa <- raw[raw$method == METHODS[2] & raw$contrast == "L3" & raw$held == STRAT, ]
pm <- merge(cbind(la, key = pk(la)), cbind(sa, key = pk(sa)), by = "key", suffixes = c("_a", "_b"))
paired <- do.call(rbind, lapply(c("jaccard", "retention", "identical", "jaccard_frozen", "jaccard_channel",
                                  "identical_channel", "top1_agreement", "top1_measured_agreement"), function(mt) {
  d <- pm[[paste0(mt, "_a")]] - pm[[paste0(mt, "_b")]]
  data.frame(metric = mt, method_a = METHODS[1], method_b = METHODS[2], contrast = "L3", held = STRAT,
             n_bag_pairs = length(d), median_difference = stats::median(d),
             p025 = unname(stats::quantile(d, 0.025)), p975 = unname(stats::quantile(d, 0.975)),
             share_a_higher = mean(d > 0), stringsAsFactors = FALSE)
}))
edge_s <- stats::aggregate(cbind(share_patients, identical) ~ method + contrast + held + band, data = edg, FUN = stats::median)
miss_s <- stats::aggregate(cbind(share_patients, jaccard, identical) ~ method + contrast + held + band, data = mis,
                           FUN = function(v) stats::median(v, na.rm = TRUE), na.action = stats::na.pass)

run <- new_run(if (SITE == "mimic") "attrbold" else "attrboldext", cfg,
               note = sprintf("bold-system agreement across replicates (%s), store %s; fits nothing", SITE, basename(gen_d)))
save_table(run, summ, "bold_agreement")
save_table(run, raw, "bold_agreement_pairs")
save_table(run, edge_s, "bold_agreement_edge")
save_table(run, paired, "bold_agreement_paired")
save_table(run, miss_s, "bold_agreement_by_missingness")
save_table(run, lead_s, "leading_channel_unmeasured")
save_table(run, COV, "replicate_coverage", subdir = "diagnostics")
finalize_run(run, extra = list(site = SITE, store = gen_d, frozen_reference_store = INT_D, store_status = STORE_STATUS, design_key = attr_key_hash(.want),
                                store_fingerprint = attr_key_hash(.fp), methods = paste(METHODS, collapse = ","),
                                pct_bold = PCT_BOLD, n_common_bags = length(COMMON_BAGS), stability_stratum = STRAT))

options(width = 200)
r3 <- function(d) { d[] <- lapply(d, function(v) if (is.numeric(v)) round(v, 3) else v); d }
cat("\n== bold-system agreement per pair of replicates: median [2.5th, 97.5th percentile] over pairs\n")
w <- summ[summ$metric %in% c("jaccard", "identical", "bold_per_patient", "jaccard_channel", "identical_channel",
                             "bold_channels_per_patient", "top1_agreement", "top1_measured_agreement"), ]
w$value <- sprintf("%.3f [%.3f, %.3f]", w$median, w$p025, w$p975)
wide <- stats::reshape(w[, c("method", "contrast", "held", "n_pairs", "metric", "value")],
                       idvar = c("method", "contrast", "held", "n_pairs"), timevar = "metric", direction = "wide")
names(wide) <- sub("^value\\.", "", names(wide))
print(wide[order(wide$contrast, wide$method), ], row.names = FALSE)
cat(sprintf("\n== %s minus %s on the same bag pairs (L3, %s)\n", METHODS[1], METHODS[2], STRAT))
print(r3(paired[, c("metric", "n_bag_pairs", "median_difference", "p025", "p975", "share_a_higher")]), row.names = FALSE)
cat("\n== identical bold-system sets by distance of the nearest measured channel to the 85th-percentile line\n")
cat("   (median over pairs; distance in percentile points, in the pair's first replicate)\n")
print(r3(edge_s[order(edge_s$contrast, edge_s$method, edge_s$held, edge_s$band), ]), row.names = FALSE)
cat("\n== bold-system agreement by number of unmeasured channels (L3, stability stratum; median over pairs)\n")
ms <- miss_s[miss_s$contrast == "L3" & miss_s$held == STRAT, ]
print(r3(ms[order(ms$band, ms$method), ]), row.names = FALSE)
cat("\n== share of patients whose leading channel (unmasked |contribution|) is unmeasured, unperturbed replicate\n")
print(r3(lead_s), row.names = FALSE)
cat("\nrun:", run$path, "\n")
