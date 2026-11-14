# tests/stab_ladder_gen.R ----------------------------------------------------------
# THE STABILITY LADDER: GENERATOR (author, 2026-10-30). A new analysis, outside
# the attribution store's contract: fitted models are saved, and the replicate
# coordinates are this script's own.
#
# QUESTION. If the package models' strata fail to transport because they read
# missingness as evidence, is the instability of their explanations from the
# same cause? Four arms are refitted on the same bootstrap bags and applied to
# the MIMIC-IV test set and to eICU:
#   score         the 19 channel models (`full` where paired, `meas` otherwise),
#                 fitted one by one through fit_one(), as the score is.
#   jgam_gated    the joint GAM under the score's missingness rule.
#   jgam_ungated  the joint GAM reading missingness as the EBM does.
#   shap          the booster on the constructed variables, read through SHAP.
# Score to gated isolates joint fitting; gated to ungated isolates reading
# missingness; ungated to SHAP isolates the model class.
#
# THE 31 COLUMNS (author, 2026-10-30). Every arm emits an n x 31 contribution
# matrix: the 19 channels, then the 12 paired interventions. Intervention
# evidence is kept whole, never redistributed:
#   score         channel = its model's intercept - logit(p_bar) + its
#                 measurement terms; intervention = the sum of that
#                 intervention's terms over every channel model that holds a
#                 copy (vasopressor: the MAP copy plus the heart-rate copy), as
#                 the score's sum adds them. Per channel model, the parts sum to
#                 its L (asserted).
#   joint GAMs    each block to its owner, once; the ungated arm's "unmeasured"
#                 constant to its channel. The single intercept belongs to no
#                 column (constant across stays, like SHAP's bias).
#   shap          each feature to the channel or intervention its name begins
#                 with; the bias column dropped. Additivity checked.
# Which entries are structurally zero by the model's own rule is saved as a mask
# per arm and site (masks/): a channel is available where measured; an
# intervention, for the score and the gated GAM, where any paired channel is
# measured, and for the ungated GAM and SHAP always. The consumer masks with it.
#
# THE BAGS. bag b = bag_of(bootstrap.seed_base + b), b = 1..40, the shared bags
# of the internal store (membership hashes checked against it); bags 20 and 26,
# excluded there, are skipped for every arm, leaving 38. Whole-bag fits, no
# folds, frozen FINAL bundle priors (the joint GAMs read the design built with
# them, out/cache/ebm/design_final). bag 0 is the unperturbed model: the bundle's
# channel models and booster, and the joint GAMs' cached final fits. SHAP has a
# second booster per bag (seed 1 = cfg$seed + levels.seed.seed_base), as in the
# eICU store, for the combined bag-and-seed reading.
#
# NO TEST OUTCOMES. The MIMIC-IV test stays' mortality is overwritten with 0 in
# the tables this script holds, before any test row is touched, so no test
# object can carry an outcome. Nothing here is a performance measure.
#
# CHECKS. Rule 1: console output is counts, timings and differences only.
#   init      the store's design equals the bundle's; the bags equal the
#             internal store's; the joint designs' rows are the test and eICU
#             stays in order, and their gates equal the measured masks.
#   bag 0     score: the bundle's 19-column L at the test set equals the test
#             look's stored L, and at eICU the eICU store's anchor, to 0.
#             SHAP: the bundle booster at eICU equals the eICU store's anchor.
#             joint: with --verify-anchor, a whole-training-set refit through
#             tests/jointgam_lib.R reproduces the cached final fit.
#   bags      score and SHAP at eICU equal the eICU store's bag replicates
#             (attrextgen_20261006T082355): the same fits, so a difference
#             means this generator is a different machine.
#
#   Rscript tests/stab_ladder_gen.R --init
#   Rscript tests/stab_ladder_gen.R --store out/runs/stabgen_... --arm score --bags 0,1,2
#   Rscript tests/stab_ladder_gen.R --store ... --arm jgam_gated --bags 0 --verify-anchor
#   Rscript tests/stab_ladder_gen.R --store ... --close
# The joint arms run on the OpenBLAS copy of R (JOINT_R below) with
# OPENBLAS_NUM_THREADS=2; the score and SHAP arms on standard R. Each refuses
# the other.
# Arms: score, shap, jgam_gated, jgam_ungated. Bags default to 0 and the 38.
# Workers on different arms or bags may run at once; every file is per
# (arm, bag, seed) and a finished one is skipped on rerun.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(targets); library(mgcv); library(xgboost); library(arrow); library(yaml); library(qs2)
})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

GEN_VERSION <- "stab_ladder_gen 2026-10-30"
args <- commandArgs(trailingOnly = TRUE)
.opt <- function(nm, default = NULL) { i <- which(args == nm); if (!length(i) || i[1] == length(args)) default else args[i[1] + 1L] }
INIT   <- "--init" %in% args
CLOSE  <- "--close" %in% args
VERIFY <- "--verify-anchor" %in% args
STORE  <- .opt("--store")
ARM    <- .opt("--arm")
ARMS   <- c("score", "shap", "jgam_gated", "jgam_ungated")
# Input runs (options added 2026-11-13): --attrgen (the internal bags), --attrextgen (the eICU store the score
# and SHAP arms reproduce) and --test (the test run whose l_test the anchor reproduces). The defaults are the
# runs of record. --init records the resolved inputs in the store (inputs.qs2); every later pass on that store
# reads them from there and stops if an option names a different run. A store made before 2026-11-13 has no
# inputs.qs2 and uses the options or the defaults, as before.
INPUT_DEFAULT <- list(attrgen = "out/runs/attrgen_20260909T115047", attrextgen = "out/runs/attrextgen_20261006T082355",
                      test = "out/runs/test_20261025T211320")
INPUTS <- local({
  saved <- if (!is.null(STORE) && file.exists(file.path(STORE, "inputs.qs2"))) qs2::qs_read(file.path(STORE, "inputs.qs2")) else NULL
  lapply(stats::setNames(names(INPUT_DEFAULT), names(INPUT_DEFAULT)), function(k) {
    given <- .opt(paste0("--", k))
    if (!is.null(saved) && !is.null(given) && !identical(given, saved[[k]]))
      stop("--", k, " ", given, " differs from the store's recorded input ", saved[[k]], call. = FALSE)
    given %||% saved[[k]] %||% INPUT_DEFAULT[[k]]
  })
})
INTERNAL  <- INPUTS$attrgen
EICU_REF  <- INPUTS$attrextgen
TEST_REF  <- file.path(INPUTS$test, "tables/l_test.rds")
JG        <- "out/cache/jointgam"
EBM       <- "out/cache/ebm"
TOL_EXACT <- 1e-8       # same code on the same data: measured 0 where checked before
TOL_JOINT <- 1e-6       # a whole-set joint refit against the cached one (multithreaded bam)
TOL_ADD   <- 1e-3       # SHAP additivity, float32 inside xgboost
# WHICH R EACH ARM RUNS ON (2026-10-30). The score and SHAP arms are checked
# exactly against the eICU store, fitted on R's reference BLAS, so they run on
# standard R. The joint arms run on the OpenBLAS copy of R at a fixed thread
# count: one and two OpenBLAS threads give joint fits that differ by about 1e-4
# in the fitted probability (measured 2026-10-30), so the count is part of the
# store's identity. The cached final joint fits were made at 2 threads.
JOINT_R       <- "out/cache/tools/R-4.5.2-openblas"
JOINT_THREADS <- "2"
if (!xor(INIT, !is.null(STORE))) stop("give --init, or --store <path>", call. = FALSE)
if (!INIT && !CLOSE && !isTRUE(ARM %in% ARMS)) stop("--arm must be one of ", paste(ARMS, collapse = ", "), call. = FALSE)

invisible(gc(reset = TRUE))
mem_gb <- function() { g <- gc(verbose = FALSE); sum(g[, ncol(g)]) / 1024 }
t_all <- start_timer()

# --- design: bundle, targets store, guards -----------------------------------------
ecfg   <- yaml::read_yaml("config/attribution_eval.yml")
xcfg   <- yaml::read_yaml("config/external.yml")
bpath  <- as.character(cfg_req(xcfg, "bundle"))
bundle <- load_bundle(bpath, verbose = FALSE)
cfg    <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
cfg_t  <- tar_read(cfg)
sigs   <- as.character(unlist(cfg$signals))
ivs    <- unique(unlist(lapply(sigs, interventions_of, cfg = cfg)))
COLS   <- c(sigs, ivs)
.dk <- intersect(BUNDLE_DESIGN_KEYS, names(cfg_t))
if (!identical(attr_key_hash(cfg_t[.dk]), attr_key_hash(bundle$cfg))) stop("the targets store's design differs from the bundle's", call. = FALSE)
priors_t <- tar_read(priors)
if (!identical(attr_key_hash(priors_final(priors_t)), attr_key_hash(bundle$priors))) stop("the targets store's final priors differ from the bundle's", call. = FALSE)

folds <- tar_read(folds)
tr    <- tar_read(train_ids)
y     <- as.integer(tar_read(y_train))
tabs_m <- tar_read(tabs)
test_ids <- folds$stay_id[folds$split == "test"]
# NO TEST OUTCOMES (header): overwritten before anything reads a test row.
.it <- match(as.character(test_ids), as.character(tabs_m$cohort$stay_id))
if (anyNA(.it)) stop("a test stay is not in the cohort table", call. = FALSE)
tabs_m$cohort$mortality[.it] <- 0L
if (!identical(as.integer(tabs_m$cohort$mortality[match(as.character(tr), as.character(tabs_m$cohort$stay_id))]), y)) {
  stop("the training outcomes changed when the test outcomes were blanked", call. = FALSE)
}
tabs_e <- load_tables(cfg$paths, cfg, site = "eicu", verbose = FALSE)
eicu_ids <- tabs_e$cohort$stay_id
SITES <- list(test = list(tabs = tabs_m, ids = test_ids), eicu = list(tabs = tabs_e, ids = eicu_ids))
for (s in names(SITES)) SITES[[s]]$meas <- measured_matrix(SITES[[s]]$tabs, cfg, SITES[[s]]$ids)

# --- the bags ---------------------------------------------------------------------
BAG_SEED0 <- as.integer(cfg_req(ecfg, "bootstrap", "seed_base"))
B_BOOT    <- as.integer(cfg_req(ecfg, "bootstrap", "b"))
SEED_B    <- as.integer(cfg_req(ecfg, "levels", "seed", "seed_base"))
bag_of <- function(seed) {                    # verbatim, tests/attr_external_bags.R L197
  gcol <- resample_cols(tabs_m$cohort, cfg)$fold_group
  src  <- if (gcol %in% names(folds)) folds else tabs_m$cohort
  gid  <- as.character(src[[gcol]])[match(tr, src$stay_id)]
  if (anyNA(gid)) stop("bag_of: `", gcol, "` is missing for some training stay", call. = FALSE)
  ug   <- unique(gid)
  with_seed(seed, { keep <- sample(ug, length(ug), replace = TRUE); gid %in% unique(keep) })
}
BAG <- lapply(seq_len(B_BOOT), function(b) bag_of(BAG_SEED0 + b))
BAG_HASH <- vapply(BAG, function(v) attr_key_hash(which(v)), "")
.ibt <- readRDS(file.path(INTERNAL, "diagnostics", "bootstrap_bags.rds"))
if (!identical(.ibt$membership_hash[match(seq_len(B_BOOT), .ibt$boot_id)], BAG_HASH)) stop("the bags differ from the internal store's", call. = FALSE)
.ix <- utils::read.csv(file.path(INTERNAL, "llr_bootstrap_excluded.csv"), stringsAsFactors = FALSE)
SKIP <- sort(unique(as.integer(.ix$boot_id[.ix$membership_hash == BAG_HASH[.ix$boot_id]])))
BAGS_ALL <- c(0L, setdiff(seq_len(B_BOOT), SKIP))
BAGS <- { a <- .opt("--bags"); if (is.null(a)) BAGS_ALL else as.integer(strsplit(a, ",", fixed = TRUE)[[1]]) }
if (!all(BAGS %in% BAGS_ALL)) stop("--bags outside 0 and the ", length(BAGS_ALL) - 1L, " bags", call. = FALSE)

# --- the joint GAMs' design ---------------------------------------------------------
feat <- readLines(file.path(EBM, "feature_names.txt"))
rows <- arrow::read_parquet(file.path(EBM, "rows.parquet"))
if (!identical(readRDS(file.path(EBM, "train_ids.rds")), as.character(tr)) || !identical(as.integer(rows$y), y)) stop("out/cache/ebm is not the store's training rows", call. = FALSE)
if (!identical(readRDS(file.path(EBM, "test_ids.rds")), as.character(test_ids))) stop("out/cache/ebm test rows are not the test stays in order", call. = FALSE)
if (!identical(readRDS(file.path(EBM, "eicu_ids.rds")), as.character(eicu_ids))) stop("out/cache/ebm eICU rows are not the eICU cohort in order", call. = FALSE)
read_design <- function(name) {
  X <- as.data.frame(arrow::read_parquet(file.path(EBM, paste0(name, ".parquet"))))
  if (!identical(names(X), feat)) stop(name, ": columns differ from feature_names.txt", call. = FALSE)
  X
}
#' The joint machinery for one variant, in its own environment.
jg_env <- function(ungated) {
  e <- new.env(parent = globalenv())
  e$cfg <- cfg_t; e$bs_ <- bam_settings(cfg_t); e$feat <- feat; e$y_all <- y
  e$UNGATED <- ungated; e$VARIANT <- if (ungated) "ungated" else "gated"
  sys.source("tests/jointgam_lib.R", envir = e)
  e
}

# --- identity: written by --init, verified by every worker --------------------------
md5 <- function(p) unname(tools::md5sum(p))
IDENT <- list(generator = GEN_VERSION, bundle = bpath, design = attr_key_hash(bundle$cfg),
              priors = attr_key_hash(bundle$priors), bags = BAG_HASH, skipped = SKIP,
              test_ids = attr_key_hash(as.character(test_ids)), eicu_ids = attr_key_hash(as.character(eicu_ids)),
              columns = COLS,
              ebm_cache = md5(file.path(EBM, c("design_final.parquet", "design_test.parquet", "design_eicu.parquet", "rows.parquet"))),
              jointgam_lib = md5("tests/jointgam_lib.R"),
              joint_r = JOINT_R, joint_threads = JOINT_THREADS,
              joint_blas = md5(file.path(JOINT_R, "bin", "x64", c("Rblas.dll", "Rlapack.dll"))))
if (INIT) {
  run <- new_run("stabgen", cfg, note = "stability ladder: whole-bag refits of four arms, applied to the MIMIC-IV test set and eICU; 31 columns")
  STORE <- run$path
  qs2::qs_save(IDENT, file.path(STORE, "identity.qs2"))
  qs2::qs_save(INPUTS, file.path(STORE, "inputs.qs2"))
  for (d in c("contrib", "masks", "diagnostics", "excluded")) dir.create(file.path(STORE, d), showWarnings = FALSE)
  dir.create(file.path("out/cache/stabladder", basename(STORE)), showWarnings = FALSE, recursive = TRUE)
} else {
  old <- qs2::qs_read(file.path(STORE, "identity.qs2"))
  bad <- names(IDENT)[!vapply(names(IDENT), function(k) identical(old[[k]], IDENT[[k]]), TRUE)]
  if (length(bad)) stop("the store was built under a different design: ", paste(bad, collapse = ", "), call. = FALSE)
  run <- structure(list(prefix = "stabgen", id = basename(STORE), path = STORE, started = Sys.time(), config = cfg,
                        log_file = file.path(STORE, "log.txt")), class = "llr_run")
}
FITS <- file.path("out/cache/stabladder", basename(STORE))
ON_OPENBLAS <- identical(normalizePath(R.home(), winslash = "/"), normalizePath(JOINT_R, winslash = "/", mustWork = TRUE))
if (!INIT && !CLOSE) {
  if (startsWith(ARM, "jgam_")) {
    if (!ON_OPENBLAS) stop("the joint arms run on ", JOINT_R, " (this is ", R.home(), ")", call. = FALSE)
    if (!identical(Sys.getenv("OPENBLAS_NUM_THREADS"), JOINT_THREADS)) stop("the joint arms run with OPENBLAS_NUM_THREADS=", JOINT_THREADS, call. = FALSE)
  } else if (ON_OPENBLAS) stop("the score and SHAP arms run on standard R, as the eICU store they reproduce", call. = FALSE)
}
cpath <- function(arm, site, b, s) file.path(STORE, "contrib", arm, site, sprintf("b%02d_s%d.qs2", b, s))
fpath <- function(arm, b, s) file.path(FITS, arm, sprintf("b%02d_s%d.qs2", b, s))
xpath <- function(arm, b) file.path(STORE, "excluded", sprintf("%s_b%02d.txt", arm, b))
diag_row <- function(arm, b, s, ...) {
  d <- data.frame(arm = arm, boot_id = b, seed_id = s, ..., r_max_gb = round(mem_gb(), 2),
                  blas = if (ON_OPENBLAS) paste0("openblas/", Sys.getenv("OPENBLAS_NUM_THREADS")) else "reference", stringsAsFactors = FALSE)
  saveRDS(d, file.path(STORE, "diagnostics", sprintf("%s_b%02d_s%d.rds", arm, b, s)))
  d
}
empty_C <- function(n) matrix(0, n, length(COLS), dimnames = list(NULL, COLS))
put <- function(C, arm, site, b, s) {
  if (!identical(dim(C), c(length(SITES[[site]]$ids), length(COLS))) || anyNA(C)) stop(arm, " ", site, " bag ", b, ": malformed contribution matrix", call. = FALSE)
  rownames(C) <- as.character(SITES[[site]]$ids)
  dir.create(dirname(cpath(arm, site, b, s)), showWarnings = FALSE, recursive = TRUE)
  qs2::qs_save(C, cpath(arm, site, b, s))
}
done <- function(arm, b, s) all(file.exists(vapply(names(SITES), function(st) cpath(arm, st, b, s), ""))) || file.exists(xpath(arm, b))
exclude <- function(arm, b, cause) {
  writeLines(substr(gsub("\\s+", " ", cause), 1L, 300L), xpath(arm, b))
  cat(sprintf("  %s bag %d EXCLUDED: %s\n", arm, b, substr(cause, 1L, 160L)))
}

# --- masks ------------------------------------------------------------------------
avail_iv <- function(meas) sapply(ivs, function(iv) {
  ps <- sigs[vapply(sigs, function(s) iv %in% interventions_of(s, cfg), TRUE)]
  rowSums(meas[, ps, drop = FALSE]) > 0
})
if (INIT) {
  for (site in names(SITES)) {
    me <- SITES[[site]]$meas
    gated <- cbind(me, avail_iv(me)); always <- cbind(me, matrix(TRUE, nrow(me), length(ivs), dimnames = list(NULL, ivs)))
    for (arm in ARMS) {
      m <- if (arm %in% c("score", "jgam_gated")) gated else always
      dimnames(m) <- list(as.character(SITES[[site]]$ids), COLS)
      qs2::qs_save(m, file.path(STORE, "masks", sprintf("%s_%s.qs2", arm, site)))
    }
    # The joint designs' gates must be the measured masks.
    e <- jg_env(FALSE); X <- read_design(paste0("design_", site)); g <- e$gates_of(X, site)
    G <- sapply(sigs, function(s) g[[s]] == 1)
    if (!identical(dim(G), dim(me)) || !all(G == unname(me[, sigs]))) stop(site, ": joint design gates differ from the measured mask", call. = FALSE)
    cat(sprintf("  %s: %d stays; joint design gates equal the measured mask\n", site, nrow(me)))
  }
  cat(sprintf("\nstore %s initialised; bags: 0 and %d (skipped %s)\n", STORE, length(BAGS_ALL) - 1L, paste(SKIP, collapse = ",")))
  quit(save = "no")
}

# --- CLOSE: coverage --------------------------------------------------------------
if (CLOSE) {
  cov <- do.call(rbind, lapply(ARMS, function(arm) {
    do.call(rbind, lapply(BAGS_ALL, function(b) { seeds <- if (arm == "shap" && b > 0L) 0:1 else 0L   # bag 0 is one unperturbed booster
      data.frame(arm = arm, boot_id = b,
      complete = all(vapply(seeds, function(s) all(file.exists(vapply(names(SITES), function(st) cpath(arm, st, b, s), ""))), TRUE)),
      excluded = file.exists(xpath(arm, b))) }))
  }))
  # The arms write different diagnostic columns: bind on their union.
  DG <- lapply(list.files(file.path(STORE, "diagnostics"), "[.]rds$", full.names = TRUE), readRDS)
  .dc <- unique(unlist(lapply(DG, names)))
  DG <- do.call(rbind, lapply(DG, function(d) { d[setdiff(.dc, names(d))] <- NA; d[.dc] }))
  save_table(run, cov, "coverage", subdir = "tables")
  if (!is.null(DG)) save_table(run, DG, "fit_diagnostics", subdir = "tables")
  common <- BAGS_ALL[vapply(BAGS_ALL, function(b) all(cov$complete[cov$boot_id == b]), TRUE)]
  print(aggregate(cbind(complete, excluded) ~ arm, cov, sum), row.names = FALSE)
  cat(sprintf("bags complete in all four arms: %d of %d (bag 0 %s)\n", sum(common > 0), length(BAGS_ALL) - 1L, if (0L %in% common) "present" else "MISSING"))
  ok <- all(cov$complete | cov$excluded)
  if (ok) finalize_run(run, extra = list(generator = GEN_VERSION, common_bags = sum(common > 0), fits = FITS)) else write_manifest(run, status = "running")
  quit(save = "no")
}

dir.create(file.path(FITS, ARM), showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== %s, bags %s ===\n", ARM, paste(BAGS, collapse = ",")))

# =================================================================================
# ARM score: 19 channel models, decomposed into the 31 columns
# =================================================================================
if (ARM == "score") {
  PRI <- lapply(stats::setNames(sigs, sigs), function(sg) priors_for(bundle$priors, sg, "final", NA_integer_))
  MD  <- vapply(sigs, function(sg) if (is.na(spec_source(sg, "full", cfg))) "full" else "meas", "")
  ND <- RW <- list()
  for (site in names(SITES)) for (sg in sigs) {
    nd <- signal_frame(sg, MD[[sg]], SITES[[site]]$tabs, cfg, PRI[[sg]], stay_ids = SITES[[site]]$ids, stage = "predict")
    ND[[site]][[sg]] <- nd
    RW[[site]][[sg]] <- match(as.character(nd$stay_id), as.character(SITES[[site]]$ids))
  }
  owner_of <- function(lab, sg) {
    v <- all.vars(str2lang(lab)); iv <- unique(sub("__.*$", "", v[grepl("__", v, fixed = TRUE)]))
    if (!length(iv)) return(sg)
    if (length(iv) != 1L || !iv %in% ivs || !iv %in% interventions_of(sg, cfg)) stop(sg, ": term ", lab, " has no single paired intervention", call. = FALSE)
    iv
  }
  #' One channel model at one site: its 31-column parts and its L.
  decomp <- function(b, sg, site) {
    nd <- ND[[site]][[sg]]
    tm <- stats::predict(b, newdata = nd, type = "terms", discrete = FALSE)
    eta <- as.numeric(stats::predict(b, newdata = nd, type = "link", discrete = FALSE))
    off <- attr(tm, "constant") - logit(PRI[[sg]]$p_bar)
    l <- eta - logit(PRI[[sg]]$p_bar)
    P <- matrix(0, nrow(nd), length(COLS), dimnames = list(NULL, COLS)); P[, sg] <- off
    for (j in colnames(tm)) { o <- owner_of(j, sg); P[, o] <- P[, o] + tm[, j] }
    d <- max(abs(rowSums(P) - l))
    if (!is.finite(d) || d > TOL_EXACT) stop(sg, " at ", site, ": the parts do not sum to L (", format(d), ")", call. = FALSE)
    list(P = P, l = l, d = d)
  }
  ref_eicu <- utils::read.csv(file.path(EICU_REF, "manifest_replicates.csv"), stringsAsFactors = FALSE)
  ref_rep <- function(method, route, b, s = 0L) {
    k <- ref_eicu$key[ref_eicu$method == method & ref_eicu$route == route & ref_eicu$boot_id == b & ref_eicu$seed_id == s]
    if (length(k) != 1L) return(NULL)
    qs2::qs_read(file.path(EICU_REF, "replicates", paste0(k, ".qs2")))
  }
  for (b in BAGS) {
    if (done("score", b, 0L)) next
    t0 <- start_timer(); models <- list(); fail <- NULL
    for (sg in sigs) {
      if (b == 0L) { models[[sg]] <- bundle$models[[paste0(sg, "/", MD[[sg]])]]; next }
      r <- try(fit_one(sg, MD[[sg]], tabs_m, cfg, PRI[[sg]], fit_ids = tr[BAG[[b]]], keep_model = TRUE, role = "final"), silent = TRUE)
      if (inherits(r, "try-error")) { fail <- paste0(sg, ": ", conditionMessage(attr(r, "condition"))); break }
      models[[sg]] <- strip_gam(r$model); rm(r)
    }
    if (!is.null(fail)) { exclude("score", b, fail); next }
    if (b > 0L) qs2::qs_save(models, fpath("score", b, 0L))
    chk <- c(sum_parts = 0, eicu_vs_store = NA_real_, test_vs_look = NA_real_)
    for (site in names(SITES)) {
      n <- length(SITES[[site]]$ids); C <- empty_C(n)
      L19 <- matrix(0, n, length(sigs), dimnames = list(as.character(SITES[[site]]$ids), sigs))
      for (sg in sigs) {
        z <- decomp(models[[sg]], sg, site); i <- RW[[site]][[sg]]
        C[i, ] <- C[i, ] + z$P; L19[i, sg] <- z$l; chk[["sum_parts"]] <- max(chk[["sum_parts"]], z$d)
      }
      if (site == "eicu") {
        ref <- ref_rep("llr_full", if (b == 0L) "ladder" else "bootstrap", b)
        if (is.null(ref)) stop("no eICU store replicate for bag ", b, call. = FALSE)
        chk[["eicu_vs_store"]] <- max(abs(ref[rownames(L19), sigs] - L19))
      }
      if (site == "test" && b == 0L) {
        lt <- readRDS(TEST_REF); dd <- 0
        for (sg in sigs) {
          z <- lt[lt$signal == sg & lt$model == MD[[sg]], ]
          if (!nrow(z) && MD[[sg]] == "full") z <- lt[lt$signal == sg & lt$model == "meas", ]   # l_long alias trap
          i <- match(as.character(z$stay_id), rownames(L19))
          if (anyNA(i) || nrow(z) != sum(SITES$test$meas[, sg])) stop("test look L does not cover ", sg, call. = FALSE)
          dd <- max(dd, abs(L19[i, sg] - z$l))
        }
        chk[["test_vs_look"]] <- dd
      }
      put(C, "score", site, b, 0L)
    }
    bad <- c(chk[["eicu_vs_store"]] > TOL_EXACT, isTRUE(chk[["test_vs_look"]] > TOL_EXACT))
    dg <- diag_row("score", b, 0L, elapsed_sec = round(t0()$elapsed_sec, 1), sum_parts_max = chk[["sum_parts"]],
                   eicu_vs_store = chk[["eicu_vs_store"]], test_vs_look = chk[["test_vs_look"]])
    cat(sprintf("  bag %2d: %.1f min; parts-L %.1e; eICU vs store %.1e%s\n", b, t0()$elapsed_sec / 60, chk[["sum_parts"]],
                chk[["eicu_vs_store"]], if (b == 0L) sprintf("; test vs test look %.1e", chk[["test_vs_look"]]) else ""))
    if (any(bad)) stop("bag ", b, ": the score does not reproduce its reference. Identify the branch before going on.", call. = FALSE)
    rm(models); invisible(gc(verbose = FALSE))
  }
}

# =================================================================================
# ARM shap: the booster, two seeds per bag, grouped to the 31 columns
# =================================================================================
if (ARM == "shap") {
  Xtr <- xgb_design_feat(tabs_m, cfg, bundle$priors, tr, role = "final")
  FN  <- bundle$xgb$xgb_feat$feature_names
  if (!identical(colnames(Xtr), FN)) stop("the training booster design's columns differ from the bundle booster's", call. = FALSE)
  GRP <- sub("__.*$", "", FN)
  if (!all(GRP %in% COLS)) stop("booster features outside the 31 columns: ", paste(unique(setdiff(GRP, COLS)), collapse = ","), call. = FALSE)
  DM <- lapply(SITES, function(S) xgboost::xgb.DMatrix(align_design(xgb_design_feat(S$tabs, cfg, bundle$priors, S$ids, role = "final",
                                                                      fold = NA_integer_, feature_names = FN), FN), missing = NA))
  grp <- patient_group_of(tabs_m$cohort, cfg, tr)
  shap31 <- function(booster, site) {
    ctr <- stats::predict(booster, DM[[site]], predcontrib = TRUE)
    if (ncol(ctr) != length(FN) + 1L) stop("SHAP: unexpected column count", call. = FALSE)
    add <- max(abs(rowSums(ctr) - stats::predict(booster, DM[[site]], outputmargin = TRUE)))
    if (!is.finite(add) || add > TOL_ADD) stop("SHAP additivity failed (", format(add), ")", call. = FALSE)
    S <- ctr[, seq_along(FN), drop = FALSE]; C <- empty_C(nrow(S))
    for (g in COLS) { j <- which(GRP == g); if (length(j)) C[, g] <- rowSums(S[, j, drop = FALSE]) }
    list(C = C, add = add)
  }
  ref_eicu <- utils::read.csv(file.path(EICU_REF, "manifest_replicates.csv"), stringsAsFactors = FALSE)
  for (b in BAGS) for (s in if (b == 0L) 0L else 0:1) {
    if (done("shap", b, s)) next
    t0 <- start_timer()
    if (b == 0L) { bst <- bundle$xgb$xgb_feat$booster; best <- bundle$xgb$xgb_feat$best_iter } else {
      r <- try(xgb_fit_full(Xtr[BAG[[b]], , drop = FALSE], y[BAG[[b]]], cfg, seed = cfg$seed + SEED_B * s, group = grp[BAG[[b]]]), silent = TRUE)
      if (inherits(r, "try-error")) { exclude("shap", b, conditionMessage(attr(r, "condition"))); break }
      bst <- r$booster; best <- r$best_iter
      writeBin(xgboost::xgb.save.raw(bst), sub("[.]qs2$", ".ubj", fpath("shap", b, s)))
    }
    add <- 0; dref <- NA_real_
    for (site in names(SITES)) {
      z <- shap31(bst, site); add <- max(add, z$add)
      if (site == "eicu") {
        rt <- if (b == 0L) "ladder" else if (s == 0L) "bootstrap" else "bootstrap_seeded"
        k <- ref_eicu$key[ref_eicu$method == "shap_xgb_feat" & ref_eicu$route == rt & ref_eicu$boot_id == b & ref_eicu$seed_id == s]
        if (length(k) != 1L) stop("no eICU store SHAP replicate for bag ", b, " seed ", s, call. = FALSE)
        ref <- qs2::qs_read(file.path(EICU_REF, "replicates", paste0(k, ".qs2")))
        dref <- max(abs(ref[, sigs] - z$C[, sigs]))
      }
      put(z$C, "shap", site, b, s)
    }
    diag_row("shap", b, s, elapsed_sec = round(t0()$elapsed_sec, 1), best_iter = as.integer(best), additivity_max = add, eicu_vs_store = dref)
    cat(sprintf("  bag %2d seed %d: %.1f min; additivity %.1e; eICU vs store %.1e\n", b, s, t0()$elapsed_sec / 60, add, dref))
    if (!is.finite(dref) || dref > TOL_EXACT) stop("bag ", b, " seed ", s, ": SHAP does not reproduce the eICU store. Identify the branch before going on.", call. = FALSE)
  }
}

# =================================================================================
# ARMS jgam_gated / jgam_ungated: the joint GAM on each bag
# =================================================================================
if (startsWith(ARM, "jgam_")) {
  ungated <- ARM == "jgam_ungated"; v <- if (ungated) "ungated" else "gated"
  e  <- jg_env(ungated)
  Xf <- read_design("design_final"); gf <- e$gates_of(Xf, "design_final")
  XS <- lapply(stats::setNames(names(SITES), names(SITES)), function(st) read_design(paste0("design_", st)))
  GS <- lapply(names(XS), function(st) e$gates_of(XS[[st]], st)); names(GS) <- names(XS)
  for (b in BAGS) {
    if (done(ARM, b, 0L) && !(b == 0L && VERIFY)) next
    t0 <- start_timer(); dref <- NA_real_
    if (b == 0L) {
      fc <- file.path(JG, v, "fits.rds")
      if (!file.exists(fc)) stop("bag 0 needs the cached final fit ", fc, "; the ", v, " run has not finished", call. = FALSE)
      fo <- readRDS(fc)$final
      if (VERIFY) {
        fr <- e$fit_joint(Xf, gf, seq_along(y), "final", NA_integer_)
        dref <- max(abs(e$predict_joint(fr, XS$test, GS$test, seq_len(nrow(XS$test))) - e$predict_joint(fo, XS$test, GS$test, seq_len(nrow(XS$test)))))
        cat(sprintf("  anchor: whole-set refit vs cached final fit, test linear predictor max|d| = %.2e (%.0f s)\n", dref, fr$seconds))
        if (!is.finite(dref) || dref > TOL_JOINT) stop("the whole-set refit does not reproduce the cached ", v, " fit", call. = FALSE)
        pt <- arrow::read_parquet(file.path(JG, v, "pred_test.parquet"))$p
        dp <- max(abs(stats::plogis(e$predict_joint(fo, XS$test, GS$test, seq_len(nrow(XS$test)))) - pt))
        if (dp > TOL_EXACT) stop("the cached fit does not reproduce its own stored test predictions (", format(dp), ")", call. = FALSE)
      }
    } else {
      fo <- try(e$fit_joint(Xf, gf, which(BAG[[b]]), "bag", b), silent = TRUE)
      if (inherits(fo, "try-error")) { exclude(ARM, b, conditionMessage(attr(fo, "condition"))); next }
      qs2::qs_save(fo, fpath(ARM, b, 0L))
    }
    for (site in names(SITES)) put(e$contrib_joint(fo, XS[[site]], GS[[site]], seq_len(nrow(XS[[site]])), COLS), ARM, site, b, 0L)
    diag_row(ARM, b, 0L, elapsed_sec = round(t0()$elapsed_sec, 1), fit_sec = round(fo$seconds, 1), edf = round(fo$edf, 2),
             converged = fo$converged, n_fit = fo$n, anchor_refit_max = dref)
    cat(sprintf("  bag %2d: %.1f min (fit %.0f s), edf %.1f, converged %s, R max %.1f GB\n", b, t0()$elapsed_sec / 60, fo$seconds, fo$edf, fo$converged, mem_gb()))
    rm(fo); invisible(gc(verbose = FALSE))
  }
}
cat(sprintf("\n%s done (%.1f min, R max %.2f GB)\n", ARM, t_all()$elapsed_sec / 60, mem_gb()))
