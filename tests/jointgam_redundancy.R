# tests/jointgam_redundancy.R ------------------------------------------------------
# WHY DO THE CURVES FLATTEN IN THE JOINT FIT? (author, 2026-10-30)
#
# In the gated joint GAM, 61 of 66 measurement curves are flatter than in the
# stand-alone channel models, including temperature, whose largest correlation
# with any single channel is 0.24. The candidate mechanism: a channel's
# stand-alone evidence is partly predictable from the OTHER channels together,
# and a joint fit keeps only the part the others do not carry. A pairwise
# correlation cannot see redundancy with a combination of channels; an R^2 on
# all of them can.
#
# Per channel j, on the training stays where j is measured, with the stored
# out-of-fold stand-alone weights of evidence (the paired model where fitted,
# the measurement model otherwise; unmeasured channels are 0, as in the score):
#   r2_linear    R^2 of L_j on the other 18 L's (least squares)
#   r2_flexible  the same with a natural cubic spline (3 df) per other channel
#   max_abs_r    the largest absolute pairwise correlation with one other L
# against the flattening of j's measurement terms, the median over its terms
# of the joint curve's swing over the stand-alone curve's (supported points,
# tests/jointgam_shape_compare.R, column AC_swing_ratio).
#
# EXPECTATION, written before running: channels with a higher R^2 flatten more
# (a negative rank correlation between r2 and the swing ratio), and temperature's
# R^2 is well above what its largest pairwise correlation suggests. If neither
# holds, redundancy does not explain the flattening and the joint fit's own
# smoothing (87 smoothing parameters estimated together) is the next suspect.
#
# Reads the targets store's l_oof and an aggregate table; prints per-channel
# aggregates only (hard rule 1).
#
#   Rscript tests/jointgam_redundancy.R [shape run directory]
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(splines) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
SHAPE <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[1] else "out/runs/jgamshape_20261030T142353" }
cfg <- tar_read(cfg); sigs <- unlist(cfg$signals)
lo  <- tar_read(l_oof); tr <- as.character(tar_read(train_ids))
L <- matrix(0, length(tr), length(sigs), dimnames = list(tr, sigs))
for (sg in sigs) {
  md <- if (is.na(spec_source(sg, "full", cfg))) "full" else "meas"
  z <- lo[lo$signal == sg & lo$model == md, ]
  if (!nrow(z) && md == "full") z <- lo[lo$signal == sg & lo$model == "meas", ]       # l_long alias trap
  i <- match(as.character(z$stay_id), tr)
  if (anyNA(i) || anyDuplicated(i)) stop(sg, ": l_oof rows do not map one-to-one to training stays", call. = FALSE)
  L[i, sg] <- z$l
}
meas <- measured_matrix(tar_read(tabs), cfg, tr)
# The stored sum must be reproduced, or this is not the score's L.
llr <- tar_read(oof_scores)$llr_sum
if (!is.null(llr) && max(abs(rowSums(L) - llr[tr])) > 1e-8) stop("rowSums(L) does not reproduce the stored llr_sum", call. = FALSE)

sh <- readRDS(file.path(SHAPE, "tables", "shape_comparison.rds"))
sh <- sh[sh$kind == "measurement", ]
res <- do.call(rbind, lapply(sigs, function(sg) {
  m <- meas[, sg]; y <- L[m, sg]; X <- L[m, setdiff(sigs, sg), drop = FALSE]
  r2l <- summary(stats::lm(y ~ X))$r.squared
  Xs <- do.call(cbind, lapply(seq_len(ncol(X)), function(j) if (stats::sd(X[, j]) > 0) splines::ns(X[, j], df = 3) else NULL))
  r2f <- summary(stats::lm(y ~ Xs))$r.squared
  r <- suppressWarnings(stats::cor(y, X)); r[is.na(r)] <- 0
  data.frame(channel = sg, n_measured = sum(m), r2_linear = r2l, r2_flexible = r2f, max_abs_r = max(abs(r)),
             partner = colnames(X)[which.max(abs(r))],
             swing_ratio_median = stats::median(sh$AC_swing_ratio[sh$channel == sg], na.rm = TRUE),
             n_terms = sum(sh$channel == sg), stringsAsFactors = FALSE)
}))
# SECOND CANDIDATE, added after the first expectation failed (2026-10-30, the R^2
# rank correlation with the swing ratio was +0.10 to +0.15): the joint fit
# rescales every channel because the stand-alone weights do not add up on the
# log-odds scale. The linear analogue of the joint fit is the logistic stack of
# the 19 L's; its coefficient on L_j is how much of j's stand-alone evidence a
# joint linear model keeps. Expectation: the stack coefficients track the
# swing ratios (positive rank correlation), and are low for temperature too.
y_tr <- as.integer(tar_read(y_train))
st <- stats::glm(y_tr ~ L, family = stats::binomial())
cf <- stats::coef(st)[-1]; names(cf) <- sigs
res$stack_coef <- unname(cf[res$channel])
res <- res[order(-res$r2_flexible), ]
options(width = 200)
r3 <- res; r3[] <- lapply(r3, function(v) if (is.numeric(v) && !is.integer(v)) round(v, 3) else v)
print(r3, row.names = FALSE)
ok <- is.finite(res$swing_ratio_median)
cat(sprintf("\nSpearman correlation across channels, R^2 vs median swing ratio: linear %.2f, flexible %.2f; max pairwise |r| vs swing ratio %.2f (n = %d)\n",
            stats::cor(res$r2_linear[ok], res$swing_ratio_median[ok], method = "spearman"),
            stats::cor(res$r2_flexible[ok], res$swing_ratio_median[ok], method = "spearman"),
            stats::cor(res$max_abs_r[ok], res$swing_ratio_median[ok], method = "spearman"), sum(ok)))
cat(sprintf("logistic stack coefficient vs median swing ratio: Spearman %.2f, Pearson %.2f; stack coefficients median %.2f [%.2f, %.2f]\n",
            stats::cor(res$stack_coef[ok], res$swing_ratio_median[ok], method = "spearman"),
            stats::cor(res$stack_coef[ok], res$swing_ratio_median[ok]),
            stats::median(res$stack_coef), min(res$stack_coef), max(res$stack_coef)))
run <- new_run("jgamredund", cfg, note = "redundancy of each channel's stand-alone evidence with the other channels, against joint-fit flattening; fits no GAM")
save_table(run, res, "redundancy_vs_flattening")
finalize_run(run, extra = list(shape_run = SHAPE))
cat("run:", run$path, "\n")
