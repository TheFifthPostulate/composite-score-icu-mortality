# tests/jointgam_copies.R -----------------------------------------------------------
# CONTROL (ii) FOR THE FIRST LADDER STEP (author, 2026-10-30). A comparator
# variant only: it never enters a bundle or a score.
#
# QUESTION. The gated joint GAM differs from the score in three ways (Codex
# review, 2026-10-30): joint estimation; one function per shared intervention,
# gated by "any paired channel measured", where the score has one copy per paired
# channel, each gated by its own channel; and no channel level. Control (i)
# (tests/score_centered_strata.R) showed the channel levels do not move the
# score's transport. This control gives the joint GAM the score's intervention
# copies, so that against the gated joint GAM it isolates the sharing, and
# against the score it isolates joint estimation.
#
# THE COPIES. Every intervention paired with more than one channel (vasopressor,
# invasive ventilation, the three sedatives, diuretic, RRT) enters once per
# paired channel: the same design column, gated by that channel alone and
# centred over that channel's measured fitting stays, as each channel model holds
# it. Unpaired-sharing terms and every measurement term are exactly the gated
# joint GAM's (tests/jointgam_lib.R, unchanged and sourced as is).
#
# IDENTIFIABILITY. In one model, copies of one function of one variable are told
# apart only by the stays where their channels' measurement differs (heart rate
# and MBP: about 10 training stays). Unregularised, their split is not
# estimable. Hence:
#   - a smooth copy carries, besides its own wiggliness penalty, a penalty on its
#     null space (mgcv's `select` construction, for copies only), so the split
#     between near-collinear copies is resolved by shrinkage, its weight
#     estimated with the smoothing parameters;
#   - the admission indicators of the copies enter as one matrix term whose
#     penalty acts on the DIFFERENCES between copies of one intervention only
#     (I - J/k per intervention, block diagonal, one weight shared across the
#     interventions: bam's setup rejects several penalties on one matrix term):
#     the copies may differ where the data say so, and each intervention's
#     common effect is unpenalised, as in the gated arm.
# These are declared departures from the gated arm, limited to the copy blocks.
#
# CHECKS (aggregate only, hard rule 1):
#   - --verify: this script's fitting function, given the gated arm's blocks (no
#     copies), reproduces the cached gated final fit (out/cache/jointgam/gated)
#     coefficient for coefficient. Run on the OpenBLAS copy of R, 2 threads, as
#     the gated arm was.
#   - every block passes the library's gate check; the manual linear predictor
#     equals bam's undiscretised prediction on the fitting rows.
#   - the summary recomputes the gated arm's strata from its cached predictions
#     and must reproduce the reported bins (ebmbins_20261030T151631) exactly.
#
#   out/cache/tools/R-4.5.2-openblas/bin/Rscript.exe tests/jointgam_copies.R --verify
#   out/cache/tools/R-4.5.2-openblas/bin/Rscript.exe tests/jointgam_copies.R
#   Rscript tests/jointgam_copies.R --summarise
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(arrow); library(mgcv); library(yaml) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
args   <- commandArgs(trailingOnly = TRUE)
VERIFY <- "--verify" %in% args
SUMM   <- "--summarise" %in% args
EBM <- "out/cache/ebm"; GATED <- "out/cache/jointgam/gated"
OUT <- "out/cache/jointgam/copies"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
# --ebmbins <run>: the reported bins to reproduce (option added 2026-11-13; default is the run of record).
.opt <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
BINS_REF <- .opt("--ebmbins", "out/runs/ebmbins_20261030T151631")

cfg  <- tar_read(cfg); bs_ <- bam_settings(cfg)
rows <- arrow::read_parquet(file.path(EBM, "rows.parquet")); y <- rows$y; fold <- rows$fold
feat <- readLines(file.path(EBM, "feature_names.txt"))
read_design <- function(name) {
  X <- as.data.frame(arrow::read_parquet(file.path(EBM, paste0(name, ".parquet"))))
  if (!identical(names(X), feat)) stop(name, ": columns differ from feature_names.txt", call. = FALSE)
  X
}
UNGATED <- FALSE; VARIANT <- "copies"
source("tests/jointgam_lib.R")
y_all <- y
N_BINS <- as.integer(cfg_req(cfg, "metrics", "n_bins"))
p_bar  <- mean(y)
to_score <- function(p) { eps <- 1e-6; logit(pmin(pmax(p, eps), 1 - eps)) - logit(p_bar) }   # as tests/ebm_strata_bins.R

# =================================================================================
# SUMMARY: strata of the copies arm against the gated arm
# =================================================================================
if (SUMM) {
  tabs <- tar_read(tabs)
  t_ids <- readRDS(file.path(EBM, "test_ids.rds"))
  y_t <- as.integer(tabs$cohort$mortality[match(t_ids, as.character(tabs$cohort$stay_id))])
  xcfg <- yaml::read_yaml("config/external.yml")
  bundle <- load_bundle(cfg_req(xcfg, "bundle"), verbose = FALSE)
  cfg_e <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
  tabs_e <- load_tables(cfg_e$paths, cfg_e, site = "eicu", verbose = FALSE)
  if (!identical(readRDS(file.path(EBM, "eicu_ids.rds")), as.character(tabs_e$cohort$stay_id))) stop("eICU order", call. = FALSE)
  Y <- list(oof = y, test = y_t, eicu = as.integer(tabs_e$cohort$mortality))
  if (anyNA(unlist(Y))) stop("an outcome is missing", call. = FALSE)
  strata <- function(s, br, size) (cut(s, br, include.lowest = TRUE, labels = FALSE) - 1L) %/% size + 1L
  ROWS <- list()
  for (arm in c("gated", "copies")) {
    d <- if (arm == "gated") GATED else OUT
    S <- lapply(c(oof = "oof", test = "test", eicu = "eicu"), function(s) to_score(arrow::read_parquet(file.path(d, paste0("pred_", s, ".parquet")))$p))
    cp <- score_cutpoints(S$oof, n_bins = N_BINS)
    if (arm == "gated") for (s in names(S)) {      # check: the reported bins, exactly
      f <- file.path(BINS_REF, "tables", sprintf("risk_bins_jgamgated_%s.csv", if (s == "oof") "oof" else paste0(s, "_frozen")))
      ref <- utils::read.csv(f); rb <- risk_bins(S[[s]], Y[[s]], n_bins = N_BINS, breaks = cp)
      if (!identical(as.integer(rb$n), as.integer(ref$n)) || !identical(as.integer(rb$deaths), as.integer(ref$deaths))) stop("gated ", s, ": the reported bins are not reproduced", call. = FALSE)
    }
    br <- cp; br[1] <- -Inf; br[length(br)] <- Inf
    for (part in c("quartile", "decile")) {
      size <- if (part == "quartile") 5L else 2L
      ref <- tapply(Y$oof, strata(S$oof, br, size), mean)
      for (st in c("test", "eicu")) {
        k <- strata(S[[st]], br, size); obs <- tapply(Y[[st]], factor(k, seq_along(ref)), mean)
        ROWS[[length(ROWS) + 1L]] <- data.frame(arm = arm, site = st, partition = part,
          max_abs_difference = max(abs(obs - ref)), top_share = mean(k == length(ref)),
          top_observed = obs[[length(obs)]], top_reference = ref[[length(ref)]], auroc = .auroc(S[[st]], Y[[st]]))
      }
    }
    ROWS[[length(ROWS) + 1L]] <- data.frame(arm = arm, site = "oof", partition = "--", max_abs_difference = NA, top_share = NA,
      top_observed = NA, top_reference = NA, auroc = .auroc(S$oof, Y$oof))
  }
  res <- do.call(rbind, ROWS)
  cat("check passed: the gated arm's reported bins are reproduced at all three sites\n")
  run <- new_run("jgamcopies", cfg, note = "control (ii) for the first ladder step: joint GAM with the score's per-channel intervention copies, against the gated joint GAM; strata on frozen out-of-fold cut points")
  save_table(run, res, "strata_summary")
  fits <- readRDS(file.path(OUT, "fits.rds"))
  save_table(run, do.call(rbind, lapply(names(fits), function(n) data.frame(fit = n, seconds = fits[[n]]$seconds, edf = fits[[n]]$edf,
    converged = fits[[n]]$converged))), "fits")
  finalize_run(run, extra = list(gated = GATED, copies = OUT, bins_reference = BINS_REF))
  options(width = 200)
  res[] <- lapply(res, function(v) if (is.numeric(v)) round(v, 3) else v)
  print(res[order(res$partition, res$site, res$arm), ], row.names = FALSE)
  cat("\nrun:", run$path, "\n")
  quit(save = "no")
}

# =================================================================================
# THE COPIES: blocks, gates, the null-space-penalised smooth class
# =================================================================================
shared <- names(blocks)[vapply(blocks, function(b) b$kind == "intervention" && length(b$gate_of) > 1L, TRUE)]
make_cb <- function(copies) {
  cb <- list()
  for (nm in names(blocks)) {
    b <- blocks[[nm]]
    if (copies && nm %in% shared) {
      for (sg in b$gate_of) { bb <- b; bb$gate_of <- sg; bb$gkey <- paste0(nm, "@", sg); bb$copy_of <- nm; cb[[bb$gkey]] <- bb }
    } else {
      b$gkey <- if (b$kind == "measurement") b$owner else nm; cb[[nm]] <- b
    }
  }
  cb
}
gates_c <- function(X, name, cb) {
  g <- gates_of(X, name)                       # the library's gates, untouched
  for (k in names(cb)) if (!is.null(cb[[k]]$copy_of)) g[[k]] <- g[[cb[[k]]$gate_of]]
  g
}
smooth.construct.gzc.smooth.spec <- function(object, data, knots) {
  object <- smooth.construct.gz.smooth.spec(object, data, knots)
  es <- eigen(object$S[[1]], symmetric = TRUE)
  nul <- es$values < max(es$values) * 1e-9
  if (any(nul)) {
    U <- es$vectors[, nul, drop = FALSE]
    object$S <- list(object$S[[1]], tcrossprod(U))
    object$rank <- c(sum(!nul), sum(nul)); object$null.space.dim <- 0L
  }
  object
}
registerS3method("smooth.construct", "gzc.smooth.spec", smooth.construct.gzc.smooth.spec, envir = asNamespace("mgcv"))

#' The library's fit_joint(), with copy blocks: same blocks, gate checks, order and
#' acceptance; copies as declared in the header. Without copies the formula is
#' fit_joint()'s exactly (checked by --verify).
fit_copies <- function(X, g, fit, role, fold_id, cb) {
  t0 <- proc.time()[["elapsed"]]
  dat <- list(y = y_all[fit]); par <- list(); grp <- list(); smo <- list()
  for (k in names(cb)) {
    b <- cb[[k]]; gt <- g[[b$gkey]][fit]
    if (sum(gt) < 50) stop(k, ": fewer than 50 gated fitting stays", call. = FALSE)
    x <- X[[b$col]][fit]; bb <- make_block(b, x, gt); M <- block_matrix(bb, x, gt)
    if (any(M[gt == 0, ] != 0)) stop(k, ": block is not zero on ungated rows", call. = FALSE)
    if (max(abs(colMeans(M[gt == 1, , drop = FALSE]))) > 1e-8) stop(k, ": block is not centred over its gated rows", call. = FALSE)
    if (bb$type == "binary") {
      if (is.null(b$copy_of)) {
        par[[k]] <- c(bb, list(key = k)); dat[[paste0("P", length(par))]] <- as.numeric(M)
      } else {
        grp[[b$copy_of]][[k]] <- c(bb, list(key = k))
      }
    } else {
      smo[[k]] <- c(bb, list(key = k, ncol = ncol(M), copy = !is.null(b$copy_of)))
      i <- length(smo); x[gt == 0] <- bb$fill
      dat[[paste0("V", i)]] <- x; dat[[paste0("G", i)]] <- gt
    }
  }
  # All admission copies form ONE matrix term (bam's discrete path takes one penalised matrix term, as the
  # ungated arm's MISSM), with one penalty per intervention, each acting on that intervention's differences.
  pp <- list()
  if (length(grp)) {
    keys <- unlist(lapply(grp, names), use.names = FALSE); nk <- length(keys)
    dat$BC <- do.call(cbind, lapply(keys, function(k) { b <- cb[[k]]; gt <- g[[b$gkey]][fit]
      as.numeric(block_matrix(fo_built_of(grp, k), X[[b$col]][fit], gt)) }))
    # One penalty, block diagonal over the interventions, with ONE weight: bam's setup rejects several
    # penalties on one matrix term (2026-10-30, Sl.setup), as the ungated arm's MISSM carries one.
    at <- 0L; S <- matrix(0, nk, nk)
    for (j in seq_along(grp)) {
      kk <- length(grp[[j]])
      S[at + seq_len(kk), at + seq_len(kk)] <- diag(kk) - matrix(1 / kk, kk, kk); at <- at + kk
    }
    pp$BC <- list(S)
  }
  rhs <- c(if (length(par)) paste0("P", seq_along(par)), if (length(grp)) "BC",
           if (length(smo)) sprintf("s(V%d, by = G%d, bs = \"%s\", xt = list(sm = smo[[%d]]$sm, Z = smo[[%d]]$Z))",
                                    seq_along(smo), seq_along(smo), ifelse(vapply(smo, `[[`, TRUE, "copy"), "gzc", "gz"),
                                    seq_along(smo), seq_along(smo)))
  f <- stats::as.formula(paste("y ~", paste(rhs, collapse = " + ")), env = environment())
  b <- mgcv::bam(f, data = dat, family = bs_$family, method = bs_$method, discrete = bs_$discrete,
                 nthreads = bs_$nthreads, select = bs_$select, gamma = bs_$gamma, gc.level = bs_$gc_level,
                 paraPen = if (length(pp)) pp else NULL, na.action = stats::na.fail)
  diag <- list(signal = "joint", model = VARIANT, role = role, fold = fold_id, converged = .bam_converged(b))
  .accept_outcome_fit(b, diag, cfg)
  # the coefficient layout, in bam's order: intercept, P terms, BC groups, smooths
  lay <- c(lapply(par, function(z) list(type = "binary", keys = z$key, ncol = 1L)),
           if (length(grp)) list(list(type = "group", keys = unlist(lapply(grp, names), use.names = FALSE),
                                      ncol = sum(lengths(grp)))),
           lapply(smo, function(z) list(type = "smooth", keys = z$key, ncol = z$ncol)))
  cf <- stats::coef(b)
  if (length(cf) != 1L + sum(vapply(lay, `[[`, 1L, "ncol"))) stop("coefficient count does not match the blocks", call. = FALSE)
  built <- c(par, unlist(unname(grp), recursive = FALSE), smo)
  fo <- list(lay = lay, built = built, coef = cf, edf = sum(b$edf), sp = b$sp, converged = diag$converged,
             n = length(fit), seconds = proc.time()[["elapsed"]] - t0)
  eta_b <- as.numeric(stats::predict(b, newdata = dat, type = "link", discrete = FALSE))
  eta_m <- predict_copies(fo, X, g, fit, cb)
  if (max(abs(eta_m - eta_b)) > 1e-8) stop("manual linear predictor differs from bam's (", format(max(abs(eta_m - eta_b))), ")", call. = FALSE)
  fo
}
fo_built_of <- function(grp, k) { for (z in grp) if (k %in% names(z)) return(z[[k]]); stop("no block ", k, call. = FALSE) }
predict_copies <- function(fo, X, g, idx, cb) {
  eta <- rep(unname(fo$coef[1]), length(idx)); at <- 1L
  for (L in fo$lay) {
    M <- do.call(cbind, lapply(L$keys, function(k) { b <- cb[[k]]; block_matrix(fo$built[[k]], X[[b$col]][idx], g[[b$gkey]][idx]) }))
    eta <- eta + as.numeric(M %*% fo$coef[at + seq_len(L$ncol)]); at <- at + L$ncol
  }
  if (anyNA(eta)) stop("a prediction is NA", call. = FALSE)
  eta
}

# =================================================================================
# --verify: without copies, this fitting function is the gated arm's
# =================================================================================
Xf <- read_design("design_final")
if (VERIFY) {
  cb0 <- make_cb(FALSE); gf <- gates_c(Xf, "design_final", cb0)
  fo <- fit_copies(Xf, gf, seq_along(y), "verify", NA_integer_, cb0)
  ref <- readRDS(file.path(GATED, "fits.rds"))$final$coef
  d <- if (length(ref) == length(fo$coef)) max(abs(unname(ref) - unname(fo$coef))) else Inf
  cat(sprintf("verify: %d coefficients, max |difference| against the cached gated final fit %.2e (%.0f s)\n", length(fo$coef), d, fo$seconds))
  if (d > 1e-8) stop("this fitting function does not reproduce the gated arm", call. = FALSE)
  quit(save = "no")
}

# =================================================================================
# the arm: five fold fits, the final fit, predictions
# =================================================================================
cb <- make_cb(TRUE)
cat(sprintf("blocks: %d (gated arm %d); shared interventions split into copies: %s\n", length(cb), length(blocks), paste(shared, collapse = ", ")))
fits <- list(); p_oof <- rep(NA_real_, length(y))
for (f in sort(unique(fold))) {
  Xk <- read_design(paste0("design_fold", f)); gk <- gates_c(Xk, paste0("design_fold", f), cb)
  fo <- fit_copies(Xk, gk, which(fold != f), "oof", as.integer(f), cb)
  ho <- which(fold == f); p_oof[ho] <- stats::plogis(predict_copies(fo, Xk, gk, ho, cb))
  fits[[paste0("fold", f)]] <- fo
  cat(sprintf("fold %d: %.0f s, edf %.1f, converged %s\n", f, fo$seconds, fo$edf, fo$converged))
}
if (anyNA(p_oof)) stop("a training stay was never scored out of fold", call. = FALSE)
arrow::write_parquet(data.frame(row = rows$row, p = p_oof), file.path(OUT, "pred_oof.parquet"))
gf <- gates_c(Xf, "design_final", cb)
fo <- fit_copies(Xf, gf, seq_along(y), "final", NA_integer_, cb); fits$final <- fo
cat(sprintf("final: %.0f s, edf %.1f, converged %s\n", fo$seconds, fo$edf, fo$converged))
for (site in c("test", "eicu")) {
  Xs <- read_design(paste0("design_", site)); gs <- gates_c(Xs, paste0("design_", site), cb)
  arrow::write_parquet(data.frame(p = stats::plogis(predict_copies(fo, Xs, gs, seq_len(nrow(Xs)), cb))), file.path(OUT, paste0("pred_", site, ".parquet")))
  cat(sprintf("%s: %d stays predicted\n", site, nrow(Xs)))
}
saveRDS(lapply(fits, function(z) z[c("coef", "edf", "sp", "converged", "n", "seconds", "lay")]), file.path(OUT, "fits.rds"))
cat("done\n")
