# tests/stab_card_marks.R ------------------------------------------------------------
# THE CARD'S BOLD MARKS UNDER RETRAINING, ON THE MIMIC-IV TEST SET (author, 2026-10-31).
#
# Results 4.2 reported the reproducibility of the card's bold marks from the retired out-of-fold bag
# protocol at MIMIC-IV (attrbold_20261028T200313, boldvals_20261028T230049). The author retired that
# protocol: this script recomputes the MIMIC-IV half on the test set from the whole-subsample refits of
# the stability ladder (tests/stab_ladder_gen.R, store stabgen_20261030T153227), the protocol the eICU
# half already used (attrboldext_20261028T215001, boldvalsext_20261028T231235: the same fits, checked
# bit for bit by the generator).
#
# WHAT IS NEW. Only the score's per-channel weight of evidence and terms at the test set: the ladder
# store holds 31 columns, in which an intervention shared by several channel models (vasopressor,
# the sedatives, renal replacement) is summed over its copies, so the card's 19 channel values cannot
# be read back from it. The saved channel models of each subsample (out/cache/stabladder/<store>/score)
# are applied to the test stays with the final priors. SHAP needs nothing new: the store's 19 channel
# columns are the card's SHAP columns (the generator checked them against the eICU attribution store).
#
# THE RULES, as the earlier analyses (no design choice is made here):
#   bold channel   tests/attr_bold_agreement.R: measured channel, round(100 * mean(ref <= x)) >= 85.
#                  reference "own" (primary): the refit's own measured values at the site;
#                  reference "frozen" (sensitivity): the card's cut-offs, the unperturbed out-of-fold
#                  replicate of the internal attribution store.
#                  Reported per channel (4.2) and per system (config/domains.csv).
#   bold value     tests/bold_values_replicates.R: channel_items() of tests/hetero_common.R with the
#                  card's threshold log(2), on terms anchored at the card's reference values.
#   pairs          all pairs of the 38 subsamples (attr_contrast_pairs on one seed): 703.
#
# CHECKS (stop on failure; identify the branch before going on, CLAUDE.md rule 12):
#   1. the store identity (bags, skipped bags, test stays) equals the live design;
#   2. per subsample, the 19 channel values of a stay sum to the sum of its 31 stored columns (1e-8);
#   3. per channel model, the anchored terms plus the reference evidence reproduce L (1e-8);
#   4. subsample 0 (the bundle) reproduces the test look's stored L (1e-8);
#   5. the pair metrics, run on the store's eICU SHAP columns, reproduce attrboldext's SHAP channel
#      Jaccard (median over the 703 pairs, own reference) to 1e-12.
#
# No test outcome is read (the outcomes are blanked as in the generator). Fits nothing.
# Aggregates only (hard rule 1).
#
#   Rscript tests/stab_card_marks.R --store out/runs/stabgen_20261030T153227 [--bags 0,1,2] [--benchmark]
#   Rscript tests/stab_card_marks.R --store out/runs/stabgen_20261030T153227 --assemble
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(mgcv); library(yaml); library(qs2) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
source("tests/hetero_common.R")

args <- commandArgs(trailingOnly = TRUE)
.opt <- function(nm, default = NULL) { i <- which(args == nm); if (!length(i) || i[1] == length(args)) default else args[i[1] + 1L] }
STORE    <- .opt("--store")
ASSEMBLE <- "--assemble" %in% args
BENCH    <- "--benchmark" %in% args
if (is.null(STORE) || !dir.exists(STORE)) stop("give --store out/runs/stabgen_<id>", call. = FALSE)
# Input runs (options added 2026-11-13; defaults are the runs of record). --attrgen and --test come from the
# store's inputs.qs2 when tests/stab_ladder_gen.R recorded one, and an option naming a different run stops the
# script; --attrboldext is the eICU bold-agreement run that check 5 reproduces.
INPUTS <- local({
  saved <- if (file.exists(file.path(STORE, "inputs.qs2"))) qs2::qs_read(file.path(STORE, "inputs.qs2")) else NULL
  dflt  <- list(attrgen = "out/runs/attrgen_20260909T115047", test = "out/runs/test_20261025T211320")
  lapply(stats::setNames(names(dflt), names(dflt)), function(k) {
    given <- .opt(paste0("--", k))
    if (!is.null(saved) && !is.null(given) && !identical(given, saved[[k]]))
      stop("--", k, " ", given, " differs from the store's recorded input ", saved[[k]], call. = FALSE)
    given %||% saved[[k]] %||% dflt[[k]]
  })
})
INTERNAL  <- INPUTS$attrgen
TEST_REF  <- file.path(INPUTS$test, "tables/l_test.rds")
OLD_EICU  <- file.path(.opt("--attrboldext", "out/runs/attrboldext_20261028T215001"), "tables/bold_agreement.csv")
TOL       <- 1e-8
STRONG    <- log(2)                 # paper/make_card_figure.R
PCT_BOLD  <- 85                     # paper/make_card_figure.R
MISS_BANDS <- c(0, 1, 3, Inf); MISS_LAB <- c("0", "1-2", ">=3")
ANCHOR_N  <- read.csv("paper/figs/coefs/anchor_n.csv", stringsAsFactors = FALSE)
n_anc_of  <- function(sg) { n <- ANCHOR_N$n_median[match(sg, ANCHOR_N$signal)]; if (is.na(n)) stop("no anchor n for ", sg); n }
t_all <- start_timer()

# --- design, as the generator ----------------------------------------------------------
ecfg   <- yaml::read_yaml("config/attribution_eval.yml")
xcfg   <- yaml::read_yaml("config/external.yml")
bpath  <- as.character(cfg_req(xcfg, "bundle"))
bundle <- load_bundle(bpath, verbose = FALSE)
cfg    <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
cfg_t  <- tar_read(cfg)
sigs   <- as.character(unlist(cfg$signals))
ID     <- qs2::qs_read(file.path(STORE, "identity.qs2"))
if (!identical(ID$bundle, bpath) || !identical(ID$design, attr_key_hash(bundle$cfg)) || !identical(ID$priors, attr_key_hash(bundle$priors))) {
  stop("the store was built on another bundle, design or priors", call. = FALSE)
}
folds  <- tar_read(folds); tr <- tar_read(train_ids)
tabs_m <- tar_read(tabs)
test_ids <- folds$stay_id[folds$split == "test"]
tabs_m$cohort$mortality[match(as.character(test_ids), as.character(tabs_m$cohort$stay_id))] <- 0L   # no test outcome
if (!identical(attr_key_hash(as.character(test_ids)), ID$test_ids)) stop("the test stays differ from the store's", call. = FALSE)
meas <- measured_matrix(tabs_m, cfg, test_ids)
ids  <- as.character(test_ids)
if (!identical(colnames(meas), sigs)) stop("measured mask columns are not the channels", call. = FALSE)
sf <- tabs_m$signal_features

BAG_SEED0 <- as.integer(cfg_req(ecfg, "bootstrap", "seed_base"))
B_BOOT    <- as.integer(cfg_req(ecfg, "bootstrap", "b"))
bag_of <- function(seed) {                    # verbatim, tests/stab_ladder_gen.R
  gcol <- resample_cols(tabs_m$cohort, cfg)$fold_group
  src  <- if (gcol %in% names(folds)) folds else tabs_m$cohort
  gid  <- as.character(src[[gcol]])[match(tr, src$stay_id)]
  ug   <- unique(gid)
  with_seed(seed, { keep <- sample(ug, length(ug), replace = TRUE); gid %in% unique(keep) })
}
BAG_HASH <- vapply(seq_len(B_BOOT), function(b) attr_key_hash(which(bag_of(BAG_SEED0 + b))), "")
if (!identical(BAG_HASH, ID$bags)) stop("the bags differ from the store's", call. = FALSE)
BAGS <- setdiff(seq_len(B_BOOT), ID$skipped)
cat(sprintf("store %s: %d subsamples (skipped %s), %d test stays\n", basename(STORE), length(BAGS), paste(ID$skipped, collapse = ","), length(ids)))

FITS  <- file.path("out/cache/stabladder", basename(STORE))
CACHE <- file.path(FITS, "card_test"); dir.create(CACHE, showWarnings = FALSE)
cpath <- function(arm, site, b) file.path(STORE, "contrib", arm, site, sprintf("b%02d_s0.qs2", b))
PRI <- lapply(stats::setNames(sigs, sigs), function(sg) priors_for(bundle$priors, sg, "final", NA_integer_))
MD  <- vapply(sigs, function(sg) if (is.na(spec_source(sg, "full", cfg))) "full" else "meas", "")

#' Anchored terms, verbatim from tests/bold_values_replicates.R.
anchored_terms <- function(b, Xp, pri, n_anc) {
  cf <- stats::coef(b); tm <- list()
  for (s in b$smooth) {
    idx  <- s$first.para:s$last.para
    row0 <- as.numeric(.smooth_rows(s, term_anchor(s$term[1], pri, n_anc)))
    tm[[s$label]] <- as.numeric(sweep(Xp[, idx, drop = FALSE], 2, row0) %*% cf[idx])
  }
  in_sm <- unlist(lapply(b$smooth, function(s) s$first.para:s$last.para))
  for (j in setdiff(seq_along(cf), c(in_sm, match("(Intercept)", names(cf))))) tm[[names(cf)[j]]] <- Xp[, j] * cf[[j]]
  list(terms = as.data.frame(tm, check.names = FALSE), ref = reference_evidence(b, pri, n_anc)$value)
}
display_cols <- function(sg, nd) {
  cols <- unlist(lapply(interventions_of(sg, cfg), function(iv)
    intersect(paste0(iv, "__", c("exposure_frac", "n_hours", "present_at_admission")), names(nd))))
  if (length(cols)) stats::setNames(as.data.frame(nd[, cols, drop = FALSE]), paste0("x:", cols)) else NULL
}
bold_mask <- function(tt, sg, stay_id) {     # verbatim, tests/bold_values_replicates.R
  z  <- sf[sf$signal == sg, c("stay_id", "n_obs", "dx_trend_defined", "q05", "q95", "value_min", "value_max")]
  zz <- z[match(as.character(stay_id), as.character(z$stay_id)), ]
  ci <- channel_items(tt, sg, zz$n_obs, as.logical(zz$dx_trend_defined), cfg, strong = STRONG,
                      tails_eq = tails_equal(zz, sg, cfg))
  it <- attr(ci, "items")
  if (ncol(it) > 30L) stop("more than 30 items on channel ", sg, call. = FALSE)
  B  <- !is.na(it) & it >= STRONG
  mask <- as.integer(B %*% (2^(seq_len(ncol(B)) - 1L)))
  if (!identical(as.integer(rowSums(B)), as.integer(ci[, "n_bold"]))) stop("bitmask and bold count disagree on ", sg, call. = FALSE)
  mask
}

# =========================================================================================
# GENERATE: per subsample, the 19 channel values and the bold-value masks at the test set
# =========================================================================================
if (!ASSEMBLE) {
  ND <- lapply(stats::setNames(sigs, sigs), function(sg) signal_frame(sg, MD[[sg]], tabs_m, cfg, PRI[[sg]], stay_ids = test_ids, stage = "predict"))
  RW <- lapply(ND, function(nd) { r <- match(as.character(nd$stay_id), ids); if (anyNA(r)) stop("frame stays outside the test set", call. = FALSE); r })
  for (sg in sigs) if (!identical(sort(RW[[sg]]), unname(which(meas[, sg])))) stop(sg, ": the frame's stays are not the measured stays", call. = FALSE)
  todo <- { a <- .opt("--bags"); if (is.null(a)) c(0L, BAGS) else as.integer(strsplit(a, ",", fixed = TRUE)[[1]]) }
  if (BENCH) todo <- todo[1]
  if (!all(todo %in% c(0L, BAGS))) stop("--bags outside 0 and the store's subsamples", call. = FALSE)
  for (b in todo) {
    f_out <- file.path(CACHE, sprintf("b%02d.qs2", b))
    if (file.exists(f_out)) next
    t0 <- start_timer()
    models <- if (b == 0L) lapply(stats::setNames(sigs, sigs), function(sg) bundle$models[[paste0(sg, "/", MD[[sg]])]]) else qs2::qs_read(file.path(FITS, "score", sprintf("b%02d_s0.qs2", b)))
    L19 <- matrix(0, length(ids), length(sigs), dimnames = list(ids, sigs)); K <- matrix(0L, length(ids), length(sigs), dimnames = list(ids, sigs))
    d_terms <- 0
    for (sg in sigs) {
      bm <- models[[sg]]; nd <- ND[[sg]]; pri <- PRI[[sg]]
      Xp <- stats::predict(bm, newdata = nd, type = "lpmatrix", discrete = FALSE)
      L  <- as.numeric(Xp %*% stats::coef(bm)) - logit(pri$p_bar)
      at <- anchored_terms(bm, Xp, pri, n_anc_of(sg))
      d_terms <- max(d_terms, abs(rowSums(at$terms) + at$ref - L))
      tt <- do.call(cbind, Filter(Negate(is.null), list(at$terms, display_cols(sg, nd))))
      L19[RW[[sg]], sg] <- L; K[RW[[sg]], sg] <- bold_mask(tt, sg, nd$stay_id)
    }
    if (!is.finite(d_terms) || d_terms > TOL) stop("subsample ", b, ": terms do not add up to L (", format(d_terms), ")", call. = FALSE)
    C <- qs2::qs_read(cpath("score", "test", b))
    if (!identical(rownames(C), ids)) stop("the store's test rows differ", call. = FALSE)
    d_sum <- max(abs(rowSums(L19) - rowSums(C)))
    if (!is.finite(d_sum) || d_sum > TOL) stop("subsample ", b, ": channel values do not sum to the stored 31 columns (", format(d_sum), ")", call. = FALSE)
    d_look <- NA_real_
    if (b == 0L) {
      lt <- readRDS(TEST_REF); d_look <- 0
      for (sg in sigs) {
        z <- lt[lt$signal == sg & lt$model == MD[[sg]], ]
        if (!nrow(z) && MD[[sg]] == "full") z <- lt[lt$signal == sg & lt$model == "meas", ]   # l_long alias trap
        i <- match(as.character(z$stay_id), ids)
        if (anyNA(i) || nrow(z) != sum(meas[, sg])) stop("test look L does not cover ", sg, call. = FALSE)
        d_look <- max(d_look, abs(L19[i, sg] - z$l))
      }
      if (d_look > TOL) stop("subsample 0 does not reproduce the test look's L (", format(d_look), ")", call. = FALSE)
    }
    qs2::qs_save(list(L19 = L19, K = K, check = data.frame(boot_id = b, terms_vs_L = d_terms, sum_vs_store = d_sum, vs_test_look = d_look)), f_out)
    cat(sprintf("  subsample %2d: %.1f min; terms %.1e; sum vs store %.1e%s\n", b, t0()$elapsed_sec / 60, d_terms, d_sum,
                if (b == 0L) sprintf("; vs test look %.1e", d_look) else ""))
    rm(models); invisible(gc(verbose = FALSE))
  }
  cat(sprintf("done in %.1f min\n", t_all()$elapsed_sec / 60))
  quit(save = "no")
}

# =========================================================================================
# ASSEMBLE: pair metrics
# =========================================================================================
POP10 <- vapply(0:1023, function(v) sum(bitwAnd(v, 2L^(0:9)) > 0), integer(1))
popcount <- function(v) POP10[v %% 1024L + 1L] + POP10[(v %/% 1024L) %% 1024L + 1L] + POP10[v %/% 1048576L + 1L]
stopifnot(popcount(c(0, 1, 1023, 1024, 2^19 - 1, 2^30 - 1)) == c(0, 1, 10, 1, 19, 30))
DOM  <- utils::read.csv("config/domains.csv", stringsAsFactors = FALSE)
DOMS <- unique(DOM$domain); dom_of <- DOM$domain[match(sigs, DOM$signal)]
if (anyNA(dom_of)) stop("domains.csv does not cover the signals", call. = FALSE)
BIT  <- 2L^(match(dom_of, DOMS) - 1L); CBIT <- 2L^(seq_along(sigs) - 1L)
pct_against <- function(x, ref_sorted) 100 * findInterval(x, ref_sorted) / length(ref_sorted)
#' Bold channels (cmask) and systems (mask), as tests/attr_bold_agreement.R encode().
encode <- function(M, ref, me) {
  mask <- integer(nrow(M)); cmask <- integer(nrow(M))
  for (j in seq_along(sigs)) {
    w <- which(me[, j]); if (!length(w)) next
    b <- round(pct_against(M[w, j], ref[[j]])) >= PCT_BOLD
    mask[w[b]] <- bitwOr(mask[w[b]], BIT[j]); cmask[w[b]] <- bitwOr(cmask[w[b]], CBIT[j])
  }
  list(mask = mask, cmask = cmask)
}
own_ref <- function(M, me) lapply(seq_along(sigs), function(j) sort(M[me[, j], j]))
pm1 <- function(a, b) {                       # one bitmask per stay, as attr_bold_agreement.R pair_metrics()
  inter <- popcount(bitwAnd(a, b)); uni <- popcount(bitwOr(a, b)); any_b <- uni > 0
  c(jaccard = mean(inter[any_b] / uni[any_b]), identical = mean(a == b),
    retention = mean(c(sum(inter) / sum(popcount(a)), sum(inter) / sum(popcount(b)))))
}
pmK <- function(A, B) {                       # n x 19 bitmasks, as bold_values_replicates.R pair_metrics()
  inter <- rowSums(matrix(popcount(bitwAnd(A, B)), nrow(A))); uni <- rowSums(matrix(popcount(bitwOr(A, B)), nrow(A))); any_b <- uni > 0
  c(jaccard = mean(inter[any_b] / uni[any_b]), identical = mean(rowSums(A != B) == 0))
}
CP <- attr_contrast_pairs(data.frame(method = "x", route = "bootstrap", boot_id = BAGS, seed_id = 0L, draw_id = 0L, stringsAsFactors = FALSE))
if (nrow(CP) != choose(length(BAGS), 2)) stop("the contrast pairs are not all pairs of subsamples", call. = FALSE)

# The card's frozen cut-offs: the unperturbed out-of-fold replicate of the internal store.
int_man  <- utils::read.csv(file.path(INTERNAL, "manifest_replicates.csv"), stringsAsFactors = FALSE)
int_meas <- qs2::qs_read(file.path(INTERNAL, "measured.qs2"))
SUB <- cfg_req(ecfg, "storage", "subdir")
frozen_ref <- function(m) {
  z <- int_man[int_man$method == m & int_man$boot_id == 0L & int_man$seed_id == 0L & int_man$draw_id == 0L, ]
  if (nrow(z) != 1L) stop("no unique unperturbed replicate for ", m, call. = FALSE)
  M <- qs2::qs_read(file.path(INTERNAL, SUB, paste0(z$key, ".qs2")))
  lapply(seq_along(sigs), function(j) sort(M[int_meas[, j], j]))
}
FROZEN <- list(llr_full = frozen_ref("llr_full"), shap_xgb_feat = frozen_ref("shap_xgb_feat"))

# Check 5 first: the metric code on the store's eICU SHAP columns against attrboldext.
me_e <- qs2::qs_read(file.path(STORE, "masks", "shap_eicu.qs2"))[, sigs]
enc_e <- lapply(stats::setNames(BAGS, BAGS), function(b) { M <- qs2::qs_read(cpath("shap", "eicu", b))[, sigs]; encode(M, own_ref(M, me_e), me_e)$cmask })
jE <- vapply(seq_len(nrow(CP)), function(q) pm1(enc_e[[as.character(CP$boot_i[q])]], enc_e[[as.character(CP$boot_j[q])]])[["jaccard"]], 0)
old <- utils::read.csv(OLD_EICU, stringsAsFactors = FALSE)
ref5 <- old$median[old$method == "shap_xgb_feat" & old$route == "refit" & old$contrast == "L3" & old$primary_stratum & old$metric == "jaccard_channel"]
if (length(ref5) != 1L || abs(stats::median(jE) - ref5) > 1e-12) stop(sprintf("check 5: eICU SHAP channel Jaccard %.6f, attrboldext %.6f", stats::median(jE), ref5), call. = FALSE)
cat(sprintf("check 5: eICU SHAP channel Jaccard %.4f equals attrboldext (%d pairs)\n", stats::median(jE), length(jE)))
rm(enc_e)

# Test set: the score from the cache, SHAP from the store.
me_t <- meas
G <- list(); CK <- list()
for (b in BAGS) {
  z <- qs2::qs_read(file.path(CACHE, sprintf("b%02d.qs2", b))); CK[[length(CK) + 1L]] <- z$check
  S <- qs2::qs_read(cpath("shap", "test", b))[, sigs]
  G[[as.character(b)]] <- list(
    llr_own  = encode(z$L19, own_ref(z$L19, me_t), me_t), llr_frz = encode(z$L19, FROZEN$llr_full, me_t),
    shap_own = encode(S, own_ref(S, me_t), me_t), shap_frz = encode(S, FROZEN$shap_xgb_feat, me_t), K = z$K)
}
CK <- do.call(rbind, CK)
z0 <- qs2::qs_read(file.path(CACHE, "b00.qs2"))$check
cat(sprintf("%d subsamples: max terms vs L %.1e, max sum vs store %.1e; subsample 0 vs test look %.1e\n",
            nrow(CK), max(CK$terms_vs_L), max(CK$sum_vs_store), z0$vs_test_look))
band <- cut(rowSums(!me_t), MISS_BANDS, MISS_LAB, right = FALSE)

raw <- do.call(rbind, lapply(seq_len(nrow(CP)), function(q) {
  A <- G[[as.character(CP$boot_i[q])]]; B <- G[[as.character(CP$boot_j[q])]]
  row <- function(method, ref, level, v) data.frame(boot_i = CP$boot_i[q], boot_j = CP$boot_j[q], method = method, reference = ref,
                                                   level = level, metric = names(v), value = unname(v), stringsAsFactors = FALSE)
  out <- list(
    row("score", "own", "channel", pm1(A$llr_own$cmask, B$llr_own$cmask)), row("score", "own", "system", pm1(A$llr_own$mask, B$llr_own$mask)),
    row("score", "frozen", "channel", pm1(A$llr_frz$cmask, B$llr_frz$cmask)), row("score", "frozen", "system", pm1(A$llr_frz$mask, B$llr_frz$mask)),
    row("shap", "own", "channel", pm1(A$shap_own$cmask, B$shap_own$cmask)), row("shap", "own", "system", pm1(A$shap_own$mask, B$shap_own$mask)),
    row("shap", "frozen", "channel", pm1(A$shap_frz$cmask, B$shap_frz$cmask)), row("shap", "frozen", "system", pm1(A$shap_frz$mask, B$shap_frz$mask)),
    row("score", "card", "value", pmK(A$K, B$K)))
  for (lv in MISS_LAB) {
    w <- band == lv
    out[[length(out) + 1L]] <- row("score", "own", paste0("channel_unmeasured_", lv), pm1(A$llr_own$cmask[w], B$llr_own$cmask[w])[c("jaccard", "identical")])
    out[[length(out) + 1L]] <- row("shap", "own", paste0("channel_unmeasured_", lv), pm1(A$shap_own$cmask[w], B$shap_own$cmask[w])[c("jaccard", "identical")])
  }
  do.call(rbind, out)
}))
key <- paste(raw$method, raw$reference, raw$level, raw$metric, sep = "|")
summ <- do.call(rbind, lapply(split(raw, key), function(d) data.frame(
  method = d$method[1], reference = d$reference[1], level = d$level[1], metric = d$metric[1], n_pairs = nrow(d),
  median = stats::median(d$value), p025 = unname(stats::quantile(d$value, 0.025)), p975 = unname(stats::quantile(d$value, 0.975)))))
# Paired: the share of pairs in which the score agrees more than SHAP, channel level, own reference.
jj <- raw[raw$level == "channel" & raw$reference == "own" & raw$metric == "jaccard", ]
w  <- merge(jj[jj$method == "score", c("boot_i", "boot_j", "value")], jj[jj$method == "shap", c("boot_i", "boot_j", "value")], by = c("boot_i", "boot_j"))
paired <- data.frame(contrast = "score vs SHAP, bold channels, own reference", n_pairs = nrow(w), share_score_higher = mean(w$value.x > w$value.y),
                     median_difference = stats::median(w$value.x - w$value.y))

run <- new_run("stabcard", cfg_t, note = sprintf("card bold marks under whole-subsample retraining, MIMIC-IV test set; store %s", basename(STORE)))
save_table(run, summ, "card_mark_agreement")
save_table(run, raw, "card_mark_agreement_pairs")
save_table(run, paired, "card_mark_paired")
save_table(run, CK, "job_checks", subdir = "diagnostics")
finalize_run(run, extra = list(store = STORE, internal = INTERNAL, test_ref = TEST_REF, check5_eicu_shap = ref5, tol = TOL, n_pairs = nrow(CP)))
options(width = 200)
s2 <- summ[summ$metric %in% c("jaccard", "identical"), ]; s2$value <- sprintf("%.3f [%.3f, %.3f]", s2$median, s2$p025, s2$p975)
print(s2[, c("method", "reference", "level", "metric", "n_pairs", "value")], row.names = FALSE)
print(paired, row.names = FALSE)
cat(sprintf("\nrun: %s (%.1f min)\n", run$path, t_all()$elapsed_sec / 60))
