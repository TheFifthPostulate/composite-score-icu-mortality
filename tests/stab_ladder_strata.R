# tests/stab_ladder_strata.R ---------------------------------------------------------
# STRATA STABILITY ON THE LADDER'S BAGS (author, 2026-10-30). Fits nothing: it
# applies the fits saved by tests/stab_ladder_gen.R.
#
# QUESTION. The ladder's transport column rests on one fit per arm. Retrain each
# arm on the 38 bootstrap bags: (1) how widely does the largest stratum gap at
# eICU (and at the MIMIC-IV test set) vary, and is the ladder's ordering of the
# arms the same in every bag; (2) does a stay land in the same quartile under two
# retrained models?
#
# REFERENCE STRATA PER BAG. Each bag is fitted on the patients it drew; its
# out-of-bag training stays (patients not drawn, about 37%) are scored by the
# bag's fit and give the bag's reference strata, as the out-of-fold training
# scores give the reported ones: 20 equal-count bins of the reference score
# (stats::quantile, as R/09_metrics.R risk_bins), quartiles = 5 bins, deciles =
# 2 bins, cut points frozen and opened to +/-Inf at the ends. The reference
# mortality of a stratum is the mortality of its out-of-bag stays. The random
# variables of every stay come from the FINAL priors, estimated on all training
# stays, out-of-bag ones included, as in the generator; the same holds for every
# arm.
#
# A score per arm, on the arm's own scale (any constant offset cancels in the
# quantiles):
#   score         the sum of the 19 channel weights of evidence
#   jgam_*        the joint GAM's linear predictor
#   shap          the booster's margin
# Bag 0 (the whole-training-set fit) has no out-of-bag stays and is left out;
# its strata are the reported ones (tests/strata_reporting.R).
#
# TEST OUTCOMES. Read here, under the 2026-10-30 rule: test-set outcomes may not
# inform any design choice of the score; this describes the frozen design.
#
# Aggregates only (hard rule 1).
#
#   Rscript tests/stab_ladder_strata.R --store out/runs/stabgen_<id> --arm <arm>   (worker: scores per bag)
#   Rscript tests/stab_ladder_strata.R --store out/runs/stabgen_<id> --summarise
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(targets); library(mgcv); library(xgboost); library(arrow); library(yaml); library(qs2)
})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
args <- commandArgs(trailingOnly = TRUE)
.opt <- function(nm) { i <- which(args == nm); if (!length(i) || i[1] == length(args)) NULL else args[i[1] + 1L] }
STORE <- .opt("--store"); ARM <- .opt("--arm"); SUMM <- "--summarise" %in% args
ARMS <- c("score", "jgam_gated", "jgam_ungated", "shap")
if (is.null(STORE) || !dir.exists(STORE)) stop("give --store out/runs/stabgen_<id>", call. = FALSE)
if (!SUMM && !isTRUE(ARM %in% ARMS)) stop("give --arm (", paste(ARMS, collapse = ", "), ") or --summarise", call. = FALSE)
FITS <- file.path("out/cache/stabladder", basename(STORE))
SC   <- file.path(FITS, "strata_scores"); dir.create(SC, showWarnings = FALSE)
ID   <- qs2::qs_read(file.path(STORE, "identity.qs2"))
JG <- "out/cache/jointgam"; EBM <- "out/cache/ebm"

# --- design, as the generator --------------------------------------------------
ecfg   <- yaml::read_yaml("config/attribution_eval.yml")
xcfg   <- yaml::read_yaml("config/external.yml")
bundle <- load_bundle(as.character(cfg_req(xcfg, "bundle")), verbose = FALSE)
cfg    <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
cfg_t  <- tar_read(cfg)
sigs   <- as.character(unlist(cfg$signals))
folds  <- tar_read(folds); tr <- tar_read(train_ids); y_tr <- as.integer(tar_read(y_train))
tabs_m <- tar_read(tabs)
test_ids <- folds$stay_id[folds$split == "test"]
tabs_e <- load_tables(cfg$paths, cfg, site = "eicu", verbose = FALSE)
eicu_ids <- tabs_e$cohort$stay_id
if (!identical(attr_key_hash(as.character(test_ids)), ID$test_ids) || !identical(attr_key_hash(as.character(eicu_ids)), ID$eicu_ids)) stop("the stays differ from the store's", call. = FALSE)
y_of <- function(tabs, ids) { v <- as.integer(tabs$cohort$mortality[match(as.character(ids), as.character(tabs$cohort$stay_id))]); if (anyNA(v)) stop("missing outcomes", call. = FALSE); v }
Y <- list(train = y_tr, test = y_of(tabs_m, test_ids), eicu = y_of(tabs_e, eicu_ids))
if (!identical(y_of(tabs_m, tr), y_tr)) stop("training outcomes differ from y_train", call. = FALSE)
SITES <- list(train = list(tabs = tabs_m, ids = tr), test = list(tabs = tabs_m, ids = test_ids), eicu = list(tabs = tabs_e, ids = eicu_ids))

BAG_SEED0 <- as.integer(cfg_req(ecfg, "bootstrap", "seed_base")); B_BOOT <- as.integer(cfg_req(ecfg, "bootstrap", "b"))
bag_of <- function(seed) {                    # verbatim, tests/stab_ladder_gen.R
  gcol <- resample_cols(tabs_m$cohort, cfg)$fold_group
  src  <- if (gcol %in% names(folds)) folds else tabs_m$cohort
  gid  <- as.character(src[[gcol]])[match(tr, src$stay_id)]
  ug   <- unique(gid)
  with_seed(seed, { keep <- sample(ug, length(ug), replace = TRUE); gid %in% unique(keep) })
}
BAG <- lapply(seq_len(B_BOOT), function(b) bag_of(BAG_SEED0 + b))
if (!identical(vapply(BAG, function(v) attr_key_hash(which(v)), ""), ID$bags)) stop("the bags differ from the store's", call. = FALSE)
BAGS <- setdiff(seq_len(B_BOOT), ID$skipped)
fpath <- function(arm, b) file.path(FITS, arm, sprintf(if (arm == "shap") "b%02d_s0.ubj" else "b%02d_s0.qs2", b))   # boosters are saved raw
spath <- function(arm, b) file.path(SC, sprintf("%s_b%02d.qs2", arm, b))

# =================================================================================
# WORKER: one score vector per site and bag
# =================================================================================
if (!SUMM) {
  t_all <- start_timer()
  if (ARM == "score") {
    PRI <- lapply(stats::setNames(sigs, sigs), function(sg) priors_for(bundle$priors, sg, "final", NA_integer_))
    MD  <- vapply(sigs, function(sg) if (is.na(spec_source(sg, "full", cfg))) "full" else "meas", "")
    ND <- RW <- list()
    for (st in names(SITES)) for (sg in sigs) {
      nd <- signal_frame(sg, MD[[sg]], SITES[[st]]$tabs, cfg, PRI[[sg]], stay_ids = SITES[[st]]$ids, stage = "predict")
      ND[[st]][[sg]] <- nd; RW[[st]][[sg]] <- match(as.character(nd$stay_id), as.character(SITES[[st]]$ids))
    }
    score_of <- function(models) lapply(stats::setNames(names(SITES), names(SITES)), function(st) {
      s <- numeric(length(SITES[[st]]$ids))
      for (sg in sigs) {
        eta <- as.numeric(stats::predict(models[[sg]], newdata = ND[[st]][[sg]], type = "link", discrete = FALSE))
        i <- RW[[st]][[sg]]; s[i] <- s[i] + eta - logit(PRI[[sg]]$p_bar)
      }
      s
    })
    # Check: the whole-set models give the stored llr_sum on the training stays? They are final fits, not
    # out of fold, so instead check bag 1 against the store's test contributions (row sums).
    loader <- function(b) qs2::qs_read(fpath("score", b))
  } else if (ARM == "shap") {
    FN <- bundle$xgb$xgb_feat$feature_names
    DM <- lapply(SITES, function(S) xgboost::xgb.DMatrix(align_design(xgb_design_feat(S$tabs, cfg, bundle$priors, S$ids, role = "final",
                                                                     fold = NA_integer_, feature_names = FN), FN), missing = NA))
    score_of <- function(bst) lapply(DM, function(d) as.numeric(stats::predict(bst, d, outputmargin = TRUE)))
    loader <- function(b) { f <- fpath("shap", b); xgboost::xgb.load.raw(readBin(f, "raw", file.size(f))) }
  } else {
    feat <- readLines(file.path(EBM, "feature_names.txt"))
    e <- new.env(parent = globalenv())
    e$cfg <- cfg_t; e$bs_ <- bam_settings(cfg_t); e$feat <- feat; e$y_all <- y_tr
    e$UNGATED <- ARM == "jgam_ungated"; e$VARIANT <- if (e$UNGATED) "ungated" else "gated"
    sys.source("tests/jointgam_lib.R", envir = e)
    rd <- function(n) { X <- as.data.frame(arrow::read_parquet(file.path(EBM, paste0(n, ".parquet")))); if (!identical(names(X), feat)) stop(n, ": columns", call. = FALSE); X }
    XS <- list(train = rd("design_final"), test = rd("design_test"), eicu = rd("design_eicu"))
    GS <- lapply(names(XS), function(st) e$gates_of(XS[[st]], st)); names(GS) <- names(XS)
    score_of <- function(fo) lapply(stats::setNames(names(XS), names(XS)), function(st) e$predict_joint(fo, XS[[st]], GS[[st]], seq_len(nrow(XS[[st]]))))
    loader <- function(b) qs2::qs_read(fpath(ARM, b))
  }
  # --reverse: a second worker on the same arm, from the last bag down; the two meet in the middle.
  # --start N: a third worker, from bag N upward. A bag another worker is still scoring is scored twice, identically.
  ORDER <- if ("--reverse" %in% args) rev(BAGS) else BAGS
  if (!is.null(.opt("--start"))) ORDER <- BAGS[BAGS >= as.integer(.opt("--start"))]
  for (b in ORDER) {
    if (file.exists(spath(ARM, b))) next
    if (!file.exists(fpath(ARM, b))) stop(ARM, " bag ", b, ": no saved fit", call. = FALSE)
    s <- score_of(loader(b))
    # The scores must order the test stays as the stored contributions do: row sums differ from the
    # score only by a constant (the intercept or SHAP's bias).
    C <- qs2::qs_read(file.path(STORE, "contrib", ARM, "test", sprintf("b%02d_s0.qs2", b)))
    d <- s$test - rowSums(C); spread <- diff(range(d))
    if (!is.finite(spread) || spread > 1e-3) stop(ARM, " bag ", b, ": scores differ from the stored contributions by more than a constant (", format(spread), ")", call. = FALSE)
    qs2::qs_save(s, spath(ARM, b))
    cat(sprintf("  %s bag %2d: done (%.1f min); offset spread vs stored contributions %.1e\n", ARM, b, t_all()$elapsed_sec / 60, spread))
  }
  cat(sprintf("%s done (%.1f min)\n", ARM, t_all()$elapsed_sec / 60))
  quit(save = "no")
}

# =================================================================================
# SUMMARY
# =================================================================================
strata_of <- function(s_ref, s, size) {
  br <- stats::quantile(s_ref, probs = seq(0, 1, length.out = 21L), names = FALSE)
  br <- unique(br); br[1] <- -Inf; br[length(br)] <- Inf
  if (length(br) != 21L) stop("tied cut points", call. = FALSE)
  (cut(s, br, include.lowest = TRUE, labels = FALSE) - 1L) %/% size + 1L
}
ROWS <- list(); ASSIGN <- list(); REF <- list(); AUC <- list()
for (a in ARMS) {
  miss <- BAGS[!file.exists(vapply(BAGS, function(b) spath(a, b), ""))]
  if (length(miss)) stop(a, ": no scores for bags ", paste(miss, collapse = ","), call. = FALSE)
  for (b in BAGS) {
    s <- qs2::qs_read(spath(a, b)); oob <- !BAG[[b]]
    # The worker checked the test scores; check eICU too, and that both sites share one offset.
    off <- vapply(c("test", "eicu"), function(st) {
      d <- s[[st]] - rowSums(qs2::qs_read(file.path(STORE, "contrib", a, st, sprintf("b%02d_s0.qs2", b))))
      if (diff(range(d)) > 1e-3) stop(a, " bag ", b, " ", st, ": scores differ from the stored contributions by more than a constant", call. = FALSE)
      mean(d)
    }, 0)
    if (abs(diff(off)) > 1e-3) stop(a, " bag ", b, ": test and eICU offsets differ", call. = FALSE)
    # Discrimination of the retrained fit, on the continuous score (R/09_metrics.R .auroc).
    AUC[[length(AUC) + 1L]] <- data.frame(arm = a, boot_id = b, site = c("oob", "test", "eicu"),
      auroc = c(.auroc(s$train[oob], Y$train[oob]), .auroc(s$test, Y$test), .auroc(s$eicu, Y$eicu)))
    for (part in c("quartile", "decile")) {
      size <- if (part == "quartile") 5L else 2L
      k_ref <- strata_of(s$train[oob], s$train[oob], size)
      ref <- tapply(Y$train[oob], k_ref, mean)
      if (part == "quartile") {
        cuts <- stats::quantile(s$train[oob], probs = c(0.25, 0.5, 0.75), names = FALSE)   # the bin edges 5, 10 and 15 of 20
        REF[[length(REF) + 1L]] <- data.frame(arm = a, boot_id = b, stratum = 1:4, reference = as.numeric(ref),
                                              upper_cut = c(cuts, NA), n_oob = sum(oob))
      }
      for (st in c("test", "eicu")) {
        k <- strata_of(s$train[oob], s[[st]], size)
        obs <- tapply(Y[[st]], factor(k, seq_along(ref)), mean)
        shr <- tabulate(k, length(ref)) / length(k)
        gap <- obs - ref
        ROWS[[length(ROWS) + 1L]] <- data.frame(arm = a, boot_id = b, site = st, partition = part, n_oob = sum(oob),
          max_abs_gap = max(abs(gap), na.rm = TRUE), top_share = shr[length(shr)], top_observed = obs[length(obs)],
          top_reference = ref[length(ref)], stringsAsFactors = FALSE)
        if (part == "quartile") ASSIGN[[paste(a, st, b)]] <- k
      }
    }
  }
}
res <- do.call(rbind, ROWS); refs <- do.call(rbind, REF); aucs <- do.call(rbind, AUC)
# Quartile assignment agreement between bags: same quartile, and within one.
AG <- list()
for (a in ARMS) for (st in c("test", "eicu")) {
  pr <- utils::combn(BAGS, 2L)
  for (q in seq_len(ncol(pr))) {
    k1 <- ASSIGN[[paste(a, st, pr[1, q])]]; k2 <- ASSIGN[[paste(a, st, pr[2, q])]]
    AG[[length(AG) + 1L]] <- data.frame(arm = a, site = st, bag_i = pr[1, q], bag_j = pr[2, q],
                                        same_quartile = mean(k1 == k2), within_one = mean(abs(k1 - k2) <= 1L))
  }
}
ag <- do.call(rbind, AG)
q3 <- function(v) c(median = stats::median(v), p025 = unname(stats::quantile(v, 0.025)), p975 = unname(stats::quantile(v, 0.975)))
summ <- do.call(rbind, lapply(split(res, paste(res$arm, res$site, res$partition)), function(d) do.call(rbind, lapply(
  c("max_abs_gap", "top_share", "top_observed", "top_reference"), function(m) data.frame(arm = d$arm[1], site = d$site[1], partition = d$partition[1],
    metric = m, n_bags = nrow(d), t(q3(100 * d[[m]])))))))
agsum <- do.call(rbind, lapply(split(ag, paste(ag$arm, ag$site)), function(d) do.call(rbind, lapply(c("same_quartile", "within_one"), function(m)
  data.frame(arm = d$arm[1], site = d$site[1], metric = m, n_pairs = nrow(d), t(q3(d[[m]])))))))
# The ladder's ordering, bag by bag: the share of bags in which each step widens the largest quartile gap.
STEPS <- list(c("jgam_gated", "score"), c("jgam_ungated", "jgam_gated"), c("shap", "jgam_ungated"))
order_s <- do.call(rbind, lapply(STEPS, function(stp) do.call(rbind, lapply(c("test", "eicu"), function(st) do.call(rbind, lapply(c("quartile", "decile"), function(part) {
  A <- res[res$arm == stp[1] & res$site == st & res$partition == part, ]; B <- res[res$arm == stp[2] & res$site == st & res$partition == part, ]
  A <- A[order(A$boot_id), ]; B <- B[order(B$boot_id), ]; stopifnot(identical(A$boot_id, B$boot_id))
  dv <- 100 * (A$max_abs_gap - B$max_abs_gap)
  data.frame(step = paste(stp[1], "minus", stp[2]), site = st, partition = part, n_bags = length(dv), t(q3(dv)), share_wider = mean(dv > 0))
}))))))

run <- new_run("stabstrata", cfg, note = sprintf("strata stability on the ladder's bags (out-of-bag reference strata), store %s; fits nothing", basename(STORE)))
save_table(run, res, "strata_by_bag")
save_table(run, summ, "strata_stability")
save_table(run, ag, "assignment_pairs")
save_table(run, refs, "reference_by_bag")
save_table(run, aucs, "auroc_by_bag")
save_table(run, agsum, "assignment_agreement")
save_table(run, order_s, "ladder_steps_by_bag")
finalize_run(run, extra = list(store = STORE, n_bags = length(BAGS), reference = "out-of-bag training stays, 20 equal-count bins"))
options(width = 200)
r2 <- function(d) { d[] <- lapply(d, function(v) if (is.numeric(v)) round(v, 2) else v); d }
cat("\n== largest stratum gap (points) and top stratum, over bags: median [2.5th, 97.5th]\n")
w <- summ; w$value <- sprintf("%.1f [%.1f, %.1f]", w$median, w$p025, w$p975); w$arm <- factor(w$arm, ARMS)
wide <- stats::reshape(w[, c("arm", "site", "partition", "metric", "value")], idvar = c("arm", "site", "partition"), timevar = "metric", direction = "wide")
names(wide) <- sub("^value\\.", "", names(wide))
print(wide[order(wide$partition, wide$site, wide$arm), ], row.names = FALSE)
cat("\n== each ladder step, bag by bag: change in the largest gap (points)\n")
print(r2(order_s), row.names = FALSE)
cat("\n== quartile assignment agreement between bags (pairs)\n")
agsum$arm <- factor(agsum$arm, ARMS)
print(r2(agsum[order(agsum$metric, agsum$site, agsum$arm), ]), row.names = FALSE)
cat("\n== AUROC of the retrained fits, over bags: median [2.5th, 97.5th]\n")
print(do.call(rbind, lapply(split(aucs, paste(aucs$arm, aucs$site)), function(d) data.frame(arm = d$arm[1], site = d$site[1], t(round(q3(d$auroc), 3))))), row.names = FALSE)
cat("\n== reference mortality (%) and upper cut point of each quartile, over bags: median / 2.5th / 97.5th\n")
print(do.call(rbind, lapply(split(refs, paste(refs$arm, refs$stratum)), function(d) data.frame(arm = d$arm[1], stratum = d$stratum[1],
  reference = paste(round(100 * q3(d$reference), 1), collapse = " / "),
  upper_cut = if (anyNA(d$upper_cut)) "" else paste(signif(q3(d$upper_cut), 3), collapse = " / ")))), row.names = FALSE)
cat("\nrun:", run$path, "\n")
