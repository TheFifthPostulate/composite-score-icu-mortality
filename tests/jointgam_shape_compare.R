# tests/jointgam_shape_compare.R ---------------------------------------------------
# SHAPE FUNCTIONS: our stand-alone channel models against the gated joint GAM
# (tests/jointgam.R), with a machinery control (author, 2026-10-30).
#
# THREE CURVES PER TERM, all anchored at the card's quiet-day reference value
# (R/15_anchors.R, term_anchor()), so each is a contribution to the log-odds
# relative to the same point and the models' different centring drops out:
#   A  the stored channel model: the bundle's final model of the channel, the
#      curve the paper shows (anchored_curve(), asserted equal to the project's
#      stored anchored grid).
#   B  the machinery control: the same `full` formula on the same measured
#      stays, fitted alone through tests/jointgam.R's hand-built blocks and
#      undiscretised bam (`--single`).
#   C  the gated joint GAM's final fit: all 19 channels together, under the
#      score's missingness rule (`out/cache/jointgam/gated`).
#   A vs B is the fitting machinery only; B vs C is the joint fit only; A vs C
#   is what a reader compares.
#
# EXPECTATIONS, written 2026-10-30 before any joint or control curve existed,
# carried over unchanged from the EBM comparison (docs/ebm_comparator_20261028.md):
#   E1 the machinery moves the curves little: A and B agree far more closely
#      than B and C.
#   E2 in the joint fit, correlated channels share their signal, so the curves
#      of the two most correlated pairs outside the Glasgow components
#      (creatinine with urea nitrogen, r = 0.73; mean blood pressure with heart
#      rate, 0.71) are flatter or reshaped in C against B.
#   E3 the control, temperature (at most 0.24 with any channel), keeps its
#   curves: C agrees with B. (This failed for the EBM, which flattened 3 of 5
#   temperature terms; the gated GAM cannot flatten through missing bins.)
#   Those five channels are the headline figure. Every channel's terms are in
#   the table, and the measurement terms of all 19 in a second figure.
#
# Per term, on supported grid points (the project's rule: at least 50 stays and
# 5 deaths in the 5% window): correlation, swing (max - min) of the second
# curve over the first, and the largest absolute difference in nats, for each
# of A-B, B-C and A-C. Admission indicators are single coefficients, not
# curves, and are not compared. A shared intervention term has one curve in C
# and one per paired channel in A and B; it is compared once per channel and
# flagged.
#
# Fits nothing; reads model parameters and aggregate grids only (hard rule 1).
#
#   Rscript tests/jointgam_shape_compare.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(yaml); library(mgcv); library(ggplot2) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

HEADLINE <- c("creatinine", "bun", "mbp", "heart_rate", "temperature")
JG <- "out/cache/jointgam"
cfg    <- tar_read(cfg)
rc     <- yaml::read_yaml("config/internal.yml")
bundle <- load_bundle(rc$test_look$bundle, cfg = cfg, strict = TRUE, verbose = FALSE)
ANCHOR_N <- read.csv("paper/figs/coefs/anchor_n.csv", stringsAsFactors = FALSE)
source("paper/domain_order.R")                                 # CH_LAB
ACV <- readRDS("paper/figs/coefs/anchored_curves.rds")
fit_B <- readRDS(file.path(JG, "single", "fits.rds"))
fit_C <- readRDS(file.path(JG, "gated", "fits.rds"))$final
if (!identical(sort(names(fit_B)), sort(cfg$signals))) stop("the control fits do not cover the 19 channels", call. = FALSE)
if (is.null(fit_C)) stop("the gated run has no final fit", call. = FALSE)

TERM_LAB <- c(pi_minus = "propensity below range", pi_plus = "propensity above range", trend = "trend",
              q05_delta = "low tail (p5) extremity", q95_delta = "high tail (p95) extremity",
              value_min_delta = "low tail (min) extremity", value_max_delta = "high tail (max) extremity")

#' One block's curve from a stripped tests/jointgam.R fit: the block's basis at
#' x (gate 1) times its coefficients. Blocks are found by design-column name and
#' their coefficients by position, as predict_joint() reads them.
curve_val <- function(fo, col, x) {
  at <- 1L
  for (bb in fo$blocks) {
    if (identical(bb$name, col)) {
      if (bb$type != "smooth") stop(col, ": not a smooth block", call. = FALSE)
      M <- mgcv::PredictMat(bb$sm, data.frame(v = x)) %*% bb$Z
      return(as.numeric(M %*% fo$coef[at + seq_len(bb$ncol)]))
    }
    at <- at + bb$ncol
  }
  stop("no block ", col, " in the fit", call. = FALSE)
}
pair_stats <- function(u, v, sp, tag) {
  su <- diff(range(u[sp])); sv <- diff(range(v[sp]))
  r  <- if (sum(sp) > 2 && sd(u[sp]) > 0 && sd(v[sp]) > 0) stats::cor(u[sp], v[sp]) else NA_real_
  stats::setNames(data.frame(r, if (su > 0.01) sv / su else NA_real_, max(abs(u[sp] - v[sp]))),
                  paste0(tag, c("_correlation", "_swing_ratio", "_max_abs_diff")))
}

curves <- list(); rows <- list()
for (sg in cfg$signals) {
  md  <- if (is.na(spec_source(sg, "full", cfg))) "full" else "meas"
  b   <- bundle$models[[paste(sg, md, sep = "/")]]
  pri <- priors_for(bundle$priors, sg, "final", NA_integer_)
  n_anc <- ANCHOR_N$n_median[match(sg, ANCHOR_N$signal)]
  for (s in b$smooth) {
    var  <- s$term[1]
    kind <- if (grepl("__", var, fixed = TRUE)) "intervention" else "measurement"
    col  <- if (kind == "intervention") var else paste0(sg, "__", var)
    g <- ACV[ACV$key == paste(sg, md, sep = "/") & ACV$term == s$label, ]
    g <- g[order(g$grid_x), ]
    if (!nrow(g)) stop("no project grid for ", sg, " ", s$label, call. = FALSE)
    grid <- g$grid_x; sp <- g$supported %in% TRUE
    x0 <- term_anchor(var, pri, n_anc)
    A <- anchored_curve(b, s, grid, x0)
    if (max(abs(A$fit - g$fit)) > 1e-8) stop("curve A differs from the stored anchored curve for ", sg, " ", s$label, call. = FALSE)
    B <- curve_val(fit_B[[sg]], col, grid) - curve_val(fit_B[[sg]], col, x0)
    C <- curve_val(fit_C, col, grid) - curve_val(fit_C, col, x0)
    shared <- kind == "intervention" && sum(vapply(cfg$signals, function(z) sub("__.*$", "", var) %in% interventions_of(z, cfg), TRUE)) > 1
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(channel = sg, term = var, kind = kind, shared = shared, supported_share_of_grid = mean(sp),
                 swing_A = diff(range(A$fit[sp])), stringsAsFactors = FALSE),
      pair_stats(A$fit, B, sp, "AB"), pair_stats(B, C, sp, "BC"), pair_stats(A$fit, C, sp, "AC"))
    if (kind == "measurement") {
      lab <- sprintf("%s: %s", CH_LAB[[sg]], TERM_LAB[[var]] %||% var)
      run_id <- cumsum(c(1, diff(sp) != 0))
      curves[[length(curves) + 1L]] <- rbind(
        data.frame(channel = sg, panel = lab, x = grid, y = A$fit, lo = A$fit - 1.96 * A$se, hi = A$fit + 1.96 * A$se,
                   model = "A. Stand-alone channel model (ours)", supported = sp, seg = paste("a", run_id)),
        data.frame(channel = sg, panel = lab, x = grid, y = B, lo = NA_real_, hi = NA_real_,
                   model = "B. Same model, joint machinery (control)", supported = sp, seg = paste("b", run_id)),
        data.frame(channel = sg, panel = lab, x = grid, y = C, lo = NA_real_, hi = NA_real_,
                   model = "C. Gated joint GAM", supported = sp, seg = paste("c", run_id)))
    }
  }
}
st <- do.call(rbind, rows); cv <- do.call(rbind, curves)

COLS <- c("A. Stand-alone channel model (ours)" = "#2a78d6", "B. Same model, joint machinery (control)" = "#8f8c85",
          "C. Gated joint GAM" = "#c9722f")
draw <- function(d) {
  d$panel <- factor(d$panel, unique(d$panel))
  p <- ggplot(d, aes(x, y, colour = model, linetype = model)) +
    geom_hline(yintercept = 0, colour = "#52514e", linewidth = 0.3) +
    geom_ribbon(data = d[!is.na(d$lo) & d$supported, ], aes(ymin = lo, ymax = hi, group = seg), fill = "#2a78d6", alpha = 0.12, colour = NA, inherit.aes = TRUE) +
    geom_line(data = d[!d$supported, ], aes(group = seg), linewidth = 0.4, alpha = 0.25) +
    geom_line(data = d[d$supported, ], aes(group = seg), linewidth = 0.6) +
    facet_wrap(~ panel, scales = "free", ncol = 4) +
    scale_colour_manual(values = COLS) +
    scale_linetype_manual(values = c("A. Stand-alone channel model (ours)" = "solid", "B. Same model, joint machinery (control)" = "22", "C. Gated joint GAM" = "solid")) +
    labs(x = "Value of the random variable (faint: outside the supported range, fewer than 50 stays or 5 deaths nearby)",
         y = "Contribution to log-odds, relative to a quiet day (nats)", colour = NULL, linetype = NULL) +
    theme_minimal(base_size = 9) + theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 7.5))
  list(p = p, h = 2.3 * ceiling(nlevels(d$panel) / 4) + 0.6)
}

run <- new_run("jgamshape", cfg, note = "shape functions: stand-alone channel models (A), machinery control (B), gated joint GAM (C); fits nothing")
save_table(run, st, "shape_comparison")
dir.create(file.path(run$path, "figs"), showWarnings = FALSE)
for (w in list(list(name = "jgam_shapes_headline", d = cv[cv$channel %in% HEADLINE, ]), list(name = "jgam_shapes_all", d = cv))) {
  fg <- draw(w$d)
  for (ext in c("png", "pdf")) ggsave(file.path(run$path, "figs", paste0(w$name, ".", ext)), fg$p, width = 10, height = fg$h, dpi = 200, limitsize = FALSE)
}
finalize_run(run, extra = list(bundle = rc$test_look$bundle, control = file.path(JG, "single"), joint = file.path(JG, "gated"),
                                headline = paste(HEADLINE, collapse = ",")))

# Console: aggregates per comparison, then the headline channels' terms.
options(width = 220)
summ <- function(d, lab) {
  q <- function(v) sprintf("%.2f [%.2f, %.2f]", median(v, na.rm = TRUE), min(v, na.rm = TRUE), max(v, na.rm = TRUE))
  data.frame(terms = lab, n = nrow(d),
             AB_max_abs_diff = q(d$AB_max_abs_diff), BC_max_abs_diff = q(d$BC_max_abs_diff),
             BC_swing_ratio = q(d$BC_swing_ratio), BC_correlation = q(d$BC_correlation),
             BC_flatter = sprintf("%d of %d", sum(d$BC_swing_ratio < 1, na.rm = TRUE), sum(!is.na(d$BC_swing_ratio))))
}
cat("median [min, max] over terms, supported points only\n")
print(rbind(summ(st[st$kind == "measurement", ], "measurement"), summ(st[st$kind == "intervention", ], "intervention")), row.names = FALSE)
hd <- st[st$channel %in% HEADLINE & st$kind == "measurement", c("channel", "term", "swing_A", "AB_max_abs_diff", "AB_correlation",
                                                                 "BC_swing_ratio", "BC_correlation", "BC_max_abs_diff")]
hd[] <- lapply(hd, function(v) if (is.numeric(v)) round(v, 3) else v)
print(hd, row.names = FALSE)
cat("\nrun:", run$path, "\n")
