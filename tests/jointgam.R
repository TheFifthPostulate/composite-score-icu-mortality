# tests/jointgam.R -----------------------------------------------------------------
# THE JOINT GAM COMPARATOR: one logistic GAM over every training stay, holding
# the terms of all 19 channels' `full` formulas, fitted jointly (author,
# 2026-10-30). A comparator arm only: it never enters a bundle and never scores
# for the paper's main model.
#
# WHY IT EXISTS. The EBM comparator reads missingness as evidence: every
# measurement covariate has a missing bin with a learned score, and every
# intervention covariate is read for every stay. The score does neither. An
# unmeasured channel contributes exactly zero, and each channel's curves are
# centred over the stays where it was measured. So the EBM-versus-score gap
# mixes two differences. The joint GAM separates them:
#
#   gated   (default)  the score's missingness rule, fitted jointly. Every term
#           of channel j is multiplied by "channel j measured" and centred over
#           the fitting stays where it was measured. An intervention enters ONCE,
#           gated on "any of its paired channels measured" and centred there. The
#           only free constant is one intercept, so "not measured" can carry no
#           evidence anywhere.
#   ungated (--ungated) the EBM's information, in the same model class:
#           - the measurement curves exactly as in the gated arm (zero when
#             unmeasured, centred over measured fitting stays);
#           - plus, per channel, one constant added whenever the channel is
#             unmeasured: the column is 1[unmeasured]. Its coefficient is the
#             log-odds shift of "not measured" against a typical measured stay.
#             Agreed design: one constant per group of EXACTLY co-missing
#             channels. No two channels are exactly co-missing in training
#             (closest: heart rate and SpO2 differ in 9 of 41,250 stays,
#             2026-10-30), so every channel has its own. The 19 constants
#             share ONE ridge penalty whose weight is estimated with the
#             smoothing parameters (author, 2026-10-30): unpenalised, the 9
#             stays missing SpO2 alone (all survivors) were separated, and the
#             SpO2 and heart-rate constants ran to -26 and +26, which gave
#             about 1,150 eICU stays a risk below 1e-6 or above 0.999;
#           - intervention terms read for EVERY stay (gate 1), centred over all
#             fitting stays, once per intervention, as the EBM reads them.
#           The EBM's several missing bins per channel always fire together,
#           so they act as one constant per channel; this arm states it.
#
# THE CONSTRUCTION, AND THE TRAP IT AVOIDS. mgcv's `s(x, by = m)` does not give
# this term. With a numeric or ordered-factor `by` the smooth is not centred, so
# its level is free and absorbs the missingness effect; with an unordered factor
# the second level's smooth learns it directly (both shown on simulated data,
# 2026-10-30). So each block is built by hand: `smoothCon()` with the centring
# constraint absorbed explicitly over the gated fitting rows; `PredictMat()` for
# every row, times the gate; penalised through `paraPen` with the smooth's own
# penalty. A binary term is (a - its mean over the gated fitting rows) * gate,
# unpenalised, as the channel formulas leave it.
#
# INPUT. The EBM's design (tests/ebm_export.R, out/cache/ebm/): exactly the
# covariates `signal_frame()` hands each channel's bam(), fold-aware, NA where a
# channel is unmeasured, interventions once and cohort-wide. Reading it rather
# than re-deriving it guarantees the joint GAM and the EBM see the same values.
# The smooth bases, k and bam settings are each channel's own (build_formula(),
# bam_settings()).
#
# CHECKS, all aggregate (hard rule 1):
#   - every channel's measurement columns are missing together, intervention
#     columns never; every term of every `full` formula maps to one column.
#   - gate check on every fit: each block is exactly 0 on ungated rows and its
#     columns average 0 over the gated fitting rows.
#   - the manual linear predictor equals bam's own undiscretised prediction on
#     the fitting rows (coefficient order and the gz basis).
#   - --check: one paired and one unpaired channel, fold 1, refitted on its
#     measured rows two ways: through its own formula (must equal the stored
#     l_oof, which proves the design columns ARE the channel frames) and through
#     the hand-built blocks (must match within a tolerance; the channel fits
#     discretise covariates, the hand-built blocks do not).
#   - `.accept_outcome_fit()` on every fit before it predicts.
#   The --check comparisons are printed, not asserted. Result of 2026-10-30
#   (mbp, fold 1): own formula vs stored l_oof 1.2e-12; hand blocks vs own
#   formula undiscretised 5.6e-3 under bam. The gap is bam's fitting path, not
#   the basis: both refitted with gam(method = "REML") agree to 4.4e-10.
#
# DISCRETISED, LIKE THE CHANNEL FITS (2026-10-30). The first version entered each
# block as a penalised parametric matrix (paraPen); with no s() terms bam ignored
# `discrete = TRUE`, and a fold took 10 minutes to over 1.5 hours. Each smooth
# block is now an s() term of the class "gz" (tests/jointgam_lib.R), carrying the
# same hand-built basis, centring and penalty, so bam fits it by its discrete
# method under bam.discrete. A final fit went from 718 s to 92 s. It does NOT
# reproduce the channel models' discretisation error: --check (2026-10-30, fold
# 1) puts the gz fit within 6.9e-3 (mbp) and 2.7e-3 (lactate) of the channel
# formula fitted undiscretised, while the channel formula's own discretised fit
# sits 0.32 (mbp, sd 0.015) and 0.11 (lactate, sd 0.002) from its undiscretised
# fit. Prediction stays undiscretised. Recorded as `discrete_effective` in
# fits.rds. Runs of record use the OpenBLAS copy of R (see
# docs/jointgam_comparator_20261030.md).
#
# OUTPUT. out/cache/jointgam/<variant>/pred_{oof,test,eicu}.parquet in the EBM's
# format, so tests/ebm_score.R can score them; fits.rds holds stripped fits
# (smooth objects without their model matrices, coefficients, edf, convergence).
#
#   Rscript tests/jointgam.R --check        reproduction + gate checks only
#   Rscript tests/jointgam.R --benchmark    time one final fit, predict nothing
#   Rscript tests/jointgam.R                five fold fits + final fit + predictions
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(arrow); library(mgcv) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args    <- commandArgs(trailingOnly = TRUE)
CHECK   <- "--check" %in% args
BENCH   <- "--benchmark" %in% args
UNGATED <- "--ungated" %in% args
if (UNGATED && CHECK) stop("--check reproduces the channel fits, which are gated; run it without --ungated", call. = FALSE)
# --single: THE MACHINERY CONTROL for the shape comparison (author, 2026-10-30).
# Each channel's own `full` model, on its own measured stays, fitted alone but
# through this script's hand-built blocks and undiscretised bam. Against the
# stored channel model it isolates the fitting machinery; against the gated
# joint GAM it isolates the joint fit. Final full-train fits only.
SINGLE  <- "--single" %in% args
if (SINGLE && (UNGATED || CHECK || BENCH)) stop("--single runs alone", call. = FALSE)
VARIANT <- if (UNGATED) "ungated" else if (SINGLE) "single" else "gated"
EBM <- "out/cache/ebm"
OUT <- file.path("out/cache/jointgam", VARIANT)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cfg  <- tar_read(cfg)
bs_  <- bam_settings(cfg)
rows <- arrow::read_parquet(file.path(EBM, "rows.parquet"))
y    <- rows$y
fold <- rows$fold
feat <- readLines(file.path(EBM, "feature_names.txt"))
read_design <- function(name) {
  X <- as.data.frame(arrow::read_parquet(file.path(EBM, paste0(name, ".parquet"))))
  if (!identical(names(X), feat)) stop(name, ": columns differ from feature_names.txt", call. = FALSE)
  X
}

# --- the machinery: block map, gates, blocks, fit, predict (tests/jointgam_lib.R) ---
source("tests/jointgam_lib.R")

y_all <- y
# The cache must be the store's: same training stays, same order.
if (!identical(readRDS(file.path(EBM, "train_ids.rds")), as.character(tar_read(train_ids)))) {
  stop("out/cache/ebm is not in the targets store's training order; re-run tests/ebm_export.R", call. = FALSE)
}

# --- --check: reproduction against the stored channel fits -------------------------
if (CHECK) {
  priors <- tar_read(priors); l_oof <- tar_read(l_oof); tr <- as.character(tar_read(train_ids))
  ids <- readRDS(file.path(EBM, "train_ids.rds"))
  if (!identical(ids, tr)) stop("the EBM design is not in the store's training order", call. = FALSE)
  X1 <- read_design("design_fold1"); g1 <- gates_of(X1, "design_fold1")
  for (sg in c("mbp", "lactate")) {
    model <- if (length(interventions_of(sg, cfg))) "full" else "meas"
    pri <- priors_for(priors, sg, "oof", 1L)
    fit <- which(fold != 1 & g1[[sg]] == 1); ho <- which(fold == 1 & g1[[sg]] == 1)
    ref <- l_oof[l_oof$signal == sg & l_oof$model == model & l_oof$fold == 1, ]
    ref <- ref$l[match(ids[ho], as.character(ref$stay_id))]
    if (anyNA(ref)) stop(sg, ": stored l_oof does not cover the held-out measured stays", call. = FALSE)
    # (1) through the channel's own formula, on frames rebuilt from the design
    use <- names(blocks)[vapply(blocks, function(b) sg %in% b$gate_of, TRUE)]
    fr <- function(idx) { d <- data.frame(mortality = y[idx])
      for (nm in use) d[[blocks[[nm]]$var]] <- X1[[nm]][idx]; d }
    b1 <- .bam_fit(build_formula(sg, "full", cfg), fr(fit), cfg)
    l1 <- as.numeric(predict(b1, newdata = fr(ho), type = "link", discrete = FALSE)) - logit(pri$p_bar)
    # (2) through the hand-built blocks, every gate 1 on the channel's measured rows
    g_one <- lapply(g1, function(v) v * 0 + (seq_along(v) %in% c(fit, ho)))
    fo <- fit_joint(X1, g_one, fit, "check", 1L, use = use)
    l2 <- predict_joint(fo, X1, g_one, ho) - logit(pri$p_bar)
    # (3) the channel's own formula without discretisation: the hand-built
    # blocks are never discretised, so THIS is what they must reproduce tightly.
    cfg_nd <- cfg; cfg_nd$bam$discrete <- FALSE
    b3 <- .bam_fit(build_formula(sg, "full", cfg), fr(fit), cfg_nd)
    l3 <- as.numeric(predict(b3, newdata = fr(ho), type = "link", discrete = FALSE)) - logit(pri$p_bar)
    cat(sprintf(paste0("%-8s %d fit / %d held-out stays (hand blocks fitted with discrete = %s)\n",
                       "  own formula vs stored l_oof:                       max |diff| %.2e\n",
                       "  hand blocks vs own formula (stored, discretised):  max |diff| %.2e, sd %.2e, cor %.6f\n",
                       "  hand blocks vs own formula undiscretised:          max |diff| %.2e\n",
                       "  own formula, discretised vs undiscretised:         max |diff| %.2e, sd %.2e\n"),
                sg, length(fit), length(ho), bs_$discrete, max(abs(l1 - ref)),
                max(abs(l2 - ref)), stats::sd(l2 - ref), stats::cor(l2, ref), max(abs(l2 - l3)),
                max(abs(l1 - l3)), stats::sd(l1 - l3)))
  }
  quit(save = "no")
}

# --- --single: the machinery control, one channel at a time ------------------------
if (SINGLE) {
  Xf <- read_design("design_final"); gf <- gates_of(Xf, "design_final")
  fits <- list()
  for (sg in cfg$signals) {
    use <- names(blocks)[vapply(blocks, function(b) sg %in% b$gate_of, TRUE)]
    # every block, interventions included, gated on THIS channel: the channel
    # model's own copy of its intervention terms, fitted on its measured stays
    g_one <- lapply(gf, function(v) gf[[sg]])
    fo <- fit_joint(Xf, g_one, which(gf[[sg]] == 1), "single", NA_integer_, use = use)
    fits[[sg]] <- fo
    cat(sprintf("%-18s %d stays, %d blocks: %.0f s, edf %.1f, converged %s\n",
                sg, fo$n, length(use), fo$seconds, fo$edf, fo$converged))
  }
  saveRDS(fits, file.path(OUT, "fits.rds"))
  cat("done\n")
  quit(save = "no")
}

# --- --benchmark: one final fit --------------------------------------------------
Xf <- read_design("design_final"); gf <- gates_of(Xf, "design_final")
cat(sprintf("blocks %d (%d smooth); training stays %d\n", length(blocks),
            sum(startsWith(vapply(blocks, `[[`, "", "term"), "s(")), length(y)))
if (BENCH) {
  fo <- fit_joint(Xf, gf, seq_along(y), "final", NA_integer_)
  cat(sprintf("final fit: %.0f s, %d coefficients, edf %.1f, converged %s\n",
              fo$seconds, length(fo$coef), fo$edf, fo$converged))
  if (!is.null(fo$miss)) print(within(fo$miss, { coef <- round(coef, 3); se <- round(se, 3) }), row.names = FALSE)
  quit(save = "no")
}

# --- the arm: five fold fits, one final fit, predictions ---------------------------
# CHECKPOINTS (author, 2026-10-30). Each fold's fit and held-out predictions are
# saved the moment the fold finishes, and the final fit likewise, so a fit that
# fails acceptance late loses only itself. A checkpoint is reused only if it was
# written under the same library, design files, settings and variant.
CKPT <- file.path(OUT, "checkpoints"); dir.create(CKPT, showWarnings = FALSE)
CK_ID <- list(variant = VARIANT, lib = unname(tools::md5sum("tests/jointgam_lib.R")),
              design = unname(tools::md5sum(file.path(EBM, c(paste0("design_fold", sort(unique(fold)), ".parquet"),
                                                         "design_final.parquet", "rows.parquet")))),
              bam = bs_[setdiff(names(bs_), "family")])
ckpt <- function(name, make) {
  p <- file.path(CKPT, paste0(name, ".rds"))
  if (file.exists(p)) {
    z <- readRDS(p)
    if (identical(z$id, CK_ID)) { cat(sprintf("%s: reused checkpoint\n", name)); return(z$value) }
    cat(sprintf("%s: checkpoint written under a different library, design or setting; refitting\n", name))
  }
  v <- make(); saveRDS(list(id = CK_ID, value = v), p); v
}
fits <- list(); p_oof <- rep(NA_real_, length(y))
for (f in sort(unique(fold))) {
  z <- ckpt(paste0("fold", f), function() {
    Xk <- read_design(paste0("design_fold", f)); gk <- gates_of(Xk, paste0("design_fold", f))
    fo <- fit_joint(Xk, gk, which(fold != f), "oof", as.integer(f))
    ho <- which(fold == f)
    list(fo = fo, ho = ho, p = stats::plogis(predict_joint(fo, Xk, gk, ho)))
  })
  p_oof[z$ho] <- z$p
  fits[[paste0("fold", f)]] <- z$fo
  cat(sprintf("fold %d: %.0f s, edf %.1f, converged %s, PQL iterations %s\n", f, z$fo$seconds, z$fo$edf, z$fo$converged,
              paste(z$fo$iter, collapse = ",")))
}
if (anyNA(p_oof)) stop("a training stay was never scored out of fold", call. = FALSE)
arrow::write_parquet(data.frame(row = rows$row, p = p_oof), file.path(OUT, "pred_oof.parquet"))

fo <- ckpt("final", function() fit_joint(Xf, gf, seq_along(y), "final", NA_integer_))
fits$final <- fo
cat(sprintf("final: %.0f s, edf %.1f, converged %s, PQL iterations %s\n", fo$seconds, fo$edf, fo$converged, paste(fo$iter, collapse = ",")))
if (!is.null(fo$miss)) print(within(fo$miss, { coef <- round(coef, 3); se <- round(se, 3) }), row.names = FALSE)
for (site in c("test", "eicu")) {
  Xs <- read_design(paste0("design_", site)); gs <- gates_of(Xs, paste0("design_", site))
  arrow::write_parquet(data.frame(p = stats::plogis(predict_joint(fo, Xs, gs, seq_len(nrow(Xs))))),
                       file.path(OUT, paste0("pred_", site, ".parquet")))
  cat(sprintf("%s: %d stays predicted\n", site, nrow(Xs)))
}
saveRDS(fits, file.path(OUT, "fits.rds"))
cat("done\n")
