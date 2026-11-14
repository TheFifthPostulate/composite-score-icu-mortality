# tests/jointgam_atlas.R ----------------------------------------------------------
# THE ATLAS OF THE GATED JOINT GAM (author, 2026-10-30): the paper's atlas
# (paper/make_atlas.R), one page per channel and every term, with each panel
# drawing the score's stand-alone paired model (the paper's curve) against the
# gated joint GAM's (tests/jointgam.R, final fit), both with 95% bands.
#
# Same grid, support shading, anchors and vocabulary as the paper's atlas: the
# vocabulary (CH, UNIT, IV, describe(), rank_of(), the theme) is read from
# paper/make_atlas.R itself, not copied. Each joint curve is anchored at the
# paper's reference value of its variable; its band is the 95% interval of the
# difference from that reference, from the fit's Bayesian covariance (Vp).
# A shared intervention (vasopressor, invasive ventilation, the sedatives,
# diuretic, RRT) has ONE curve in the joint GAM, entered once and gated on any
# of its paired channels being measured; it is drawn on every page that holds
# it, and its title says so. The stand-alone model has its own copy per channel.
# Admission indicators are single coefficients with their intervals.
#
# Asserted: curve A equals the paper's stored anchored curve (it is read from
# it); the joint fit carries its covariance and its blocks match the design.
# Fits nothing; reads model parameters and aggregate grids only (hard rule 1).
#
#   Rscript tests/jointgam_atlas.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(mgcv); library(ggplot2); library(patchwork) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

# The paper atlas's vocabulary and theme, evaluated from its own source.
.want <- c("INK2", "GRID", "thm", "CH", "UNIT", "IV", "var_of", "describe", "rank_of", "wrap_title")
for (e in parse("paper/make_atlas.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% .want) eval(e, globalenv())
}
if (!all(vapply(.want, exists, TRUE))) stop("paper/make_atlas.R no longer defines ", paste(.want[!vapply(.want, exists, TRUE)], collapse = ", "), call. = FALSE)

cfg  <- tar_read(cfg)
JG   <- "out/cache/jointgam/gated"
fo   <- readRDS(file.path(JG, "fits.rds"))$final
if (is.null(fo$Vp)) stop("the gated final fit carries no covariance; refit with the current tests/jointgam_lib.R", call. = FALSE)
if (!identical(dim(fo$Vp), rep(length(fo$coef), 2L))) stop("covariance and coefficients disagree in size", call. = FALSE)
cv   <- readRDS("paper/figs/coefs/anchored_curves.rds")
ANCH <- read.csv("paper/figs/coefs/anchors.csv", stringsAsFactors = FALSE)
pc   <- read.csv("paper/figs/coefs/parametric_coefs.csv", stringsAsFactors = FALSE)
pc   <- pc[!duplicated(pc[, c("key", "term")]), ]
anchor_at <- function(key, term) ANCH$anchor[ANCH$key == key & ANCH$term == term][1]

# Positions of every block's coefficients in the joint fit (intercept first).
POS <- list(); at <- 1L
for (bb in fo$blocks) { POS[[bb$name]] <- at + seq_len(bb$ncol); at <- at + bb$ncol }
BLK <- stats::setNames(fo$blocks, vapply(fo$blocks, `[[`, "", "name"))
paired_with <- function(iv) cfg$signals[vapply(cfg$signals, function(s) iv %in% interventions_of(s, cfg), TRUE)]

#' The joint curve of design column `col` at x, anchored at x0, with its SE.
joint_curve <- function(col, x, x0) {
  bb <- BLK[[col]]; if (is.null(bb) || bb$type != "smooth") stop("no smooth block ", col, " in the joint fit", call. = FALSE)
  B  <- mgcv::PredictMat(bb$sm, data.frame(v = x))
  B0 <- mgcv::PredictMat(bb$sm, data.frame(v = x0))            # 1 row: the basis at the reference value
  D  <- sweep(B, 2L, as.numeric(B0)) %*% bb$Z                    # basis differences, centred columns
  i <- POS[[col]]
  list(fit = as.numeric(D %*% fo$coef[i]), se = sqrt(pmax(rowSums((D %*% fo$Vp[i, i, drop = FALSE]) * D), 0)))
}

LAB_A <- "Stand-alone paired model (score)"; LAB_C <- "Gated joint GAM"
COL <- stats::setNames(c("#2a78d6", "#c9722f"), c(LAB_A, LAB_C))

smooth_panel <- function(sg, key, term) {
  z <- cv[cv$key == key & cv$term == term, ]; z <- z[order(z$grid_x), ]
  if (!nrow(z)) stop("no paper curve for ", key, " ", term, call. = FALSE)
  var <- var_of(term); col <- if (grepl("__", var)) var else paste0(sg, "__", var)
  x0 <- anchor_at(key, term); jc <- joint_curve(col, z$grid_x, x0)
  run <- cumsum(c(1L, diff(as.integer(z$supported)) != 0L))
  d <- rbind(data.frame(grid_x = z$grid_x, fit = z$fit, lo = z$fit - 1.96 * z$se, hi = z$fit + 1.96 * z$se, supported = z$supported, run = run, m = LAB_A),
             data.frame(grid_x = z$grid_x, fit = jc$fit, lo = jc$fit - 1.96 * jc$se, hi = jc$fit + 1.96 * jc$se, supported = z$supported, run = run, m = LAB_C))
  d$m <- factor(d$m, c(LAB_A, LAB_C))
  half <- if (nrow(z) > 1L) stats::median(diff(z$grid_x)) / 2 else 0
  r <- rle(!z$supported); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  rects <- data.frame(xmin = z$grid_x[starts[r$values]] - half, xmax = z$grid_x[ends[r$values]] + half)
  lab <- describe(sg, term)
  shared <- if (grepl("__", var)) setdiff(paired_with(sub("__.*", "", var)), sg) else character(0)
  ttl <- paste0(wrap_title(lab$title), if (length(shared)) paste0("
joint curve shared with ", tolower(paste(CH[shared], collapse = ", "))) else "")
  g <- ggplot(d, aes(grid_x, fit, colour = m, fill = m))
  if (nrow(rects)) g <- g + geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "#f0efec")
  g <- g + geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
    geom_vline(xintercept = x0, colour = INK2, linewidth = 0.3, linetype = "22") +
    geom_ribbon(aes(ymin = lo, ymax = hi, group = m), colour = NA, alpha = 0.14)
  g <- if (z$kind[1] == "atoms") g + geom_point(aes(shape = supported), size = 1.2) +
         scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none")
       else g + geom_line(aes(linetype = supported, group = interaction(m, run)), linewidth = 0.6) +
         scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dotted"), guide = "none")
  list(plot = g + scale_colour_manual(values = COL, guide = "none") + scale_fill_manual(values = COL, guide = "none") +
         expand_limits(y = c(-0.25, 0.25)) + labs(title = ttl, x = lab$x, y = "Relative to reference (nats)") + thm,
       row = data.frame(channel = sg, term = term, shared = length(shared) > 0, swing_A = diff(range(z$fit[z$supported])),
                        swing_C = diff(range(jc$fit[z$supported])), max_abs_diff = max(abs(z$fit - jc$fit)[z$supported]),
                        share_C_inside_A_band = mean((jc$fit >= z$fit - 1.96 * z$se & jc$fit <= z$fit + 1.96 * z$se)[z$supported])))
}
coef_panel <- function(sg, key, term) {
  a <- pc[pc$key == key & pc$term == term, ]
  if (nrow(a) != 1L) stop("no paper coefficient for ", key, " ", term, call. = FALSE)
  i <- POS[[term]]; if (length(i) != 1L || BLK[[term]]$type != "binary") stop("no admission block ", term, " in the joint fit", call. = FALSE)
  d <- data.frame(m = factor(c(LAB_A, LAB_C), c(LAB_A, LAB_C)), fit = c(a$estimate, fo$coef[i]), se = c(a$se, sqrt(fo$Vp[i, i])))
  shared <- setdiff(paired_with(sub("__.*", "", term)), sg)
  ttl <- paste0(wrap_title(describe(sg, term)$title), if (length(shared)) paste0("
joint coefficient shared with ", tolower(paste(CH[shared], collapse = ", "))) else "")
  list(plot = ggplot(d, aes(m, fit, colour = m)) + geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
         geom_errorbar(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), width = 0.15) + geom_point(size = 1.8) +
         scale_colour_manual(values = COL, guide = "none") + expand_limits(y = c(-0.25, 0.25)) +
         labs(title = ttl, x = "Blue: stand-alone.  Orange: joint.", y = "Relative to reference (nats)") +
         thm + theme(axis.text.x = element_blank()),
       row = data.frame(channel = sg, term = term, shared = length(shared) > 0, swing_A = abs(a$estimate), swing_C = abs(fo$coef[i]),
                        max_abs_diff = abs(a$estimate - fo$coef[i]), share_C_inside_A_band = as.numeric(abs(fo$coef[i] - a$estimate) <= 1.96 * a$se)))
}

run <- new_run("jgamatlas", cfg, note = "atlas: stand-alone paired models against the gated joint GAM, every term; fits nothing")
FIGS <- file.path(run$path, "figs"); dir.create(FIGS, showWarnings = FALSE)
rows <- list(); pages <- list()
for (sg in names(CH)) {
  key <- if (is.na(spec_source(sg, "full", cfg))) paste0(sg, "/full") else paste0(sg, "/meas")
  terms <- unique(cv$term[cv$key == key]); paras <- unique(pc$term[pc$key == key])
  all_t <- c(terms, paras); all_t <- all_t[order(vapply(all_t, rank_of, numeric(1)))]
  out <- lapply(all_t, function(t) if (t %in% paras) coef_panel(sg, key, t) else smooth_panel(sg, key, t))
  rows <- c(rows, lapply(out, `[[`, "row"))
  p <- wrap_plots(lapply(out, `[[`, "plot"), ncol = 3) +
    plot_annotation(title = CH[[sg]],
                    subtitle = "Blue: the score's stand-alone paired model (the paper's atlas).  Orange: the gated joint GAM.
Bands: 95% intervals of the difference from the reference value (dotted line).  Admission indicators: present at admission against not.",
                    theme = theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(colour = INK2, size = 8)))
  h <- min(9.2, 2.35 * ceiling(length(out) / 3) + 0.9)
  f <- file.path(FIGS, paste0("jgam_atlas_", sg))
  ragg::agg_png(paste0(f, ".png"), width = 8, height = h, units = "in", res = 200); print(p); dev.off()
  pages[[sg]] <- list(p = p, h = h)
  cat(sprintf("%-18s %2d panels\n", sg, length(out)))
}
grDevices::cairo_pdf(file.path(FIGS, "jgam_atlas.pdf"), width = 8, height = 9.2, onefile = TRUE)
for (pg in pages) print(pg$p)
dev.off()
st <- do.call(rbind, rows); st$swing_ratio_C_over_A <- ifelse(st$swing_A > 0.01, st$swing_C / st$swing_A, NA_real_)
save_table(run, st, "atlas_terms")
finalize_run(run, extra = list(joint = JG, paper_curves = "paper/figs/coefs/anchored_curves.rds"))
options(width = 200)
cat("\nper term (supported points): swing of the joint curve over the stand-alone one, and whether it stays inside the stand-alone band\n")
agg <- function(d) sprintf("n %3d | swing ratio median %.2f | flatter %d | inside A's band (median share) %.2f",
                           nrow(d), median(d$swing_ratio_C_over_A, na.rm = TRUE), sum(d$swing_ratio_C_over_A < 1, na.rm = TRUE),
                           median(d$share_C_inside_A_band, na.rm = TRUE))
iv <- grepl("__", st$term)
cat("  measurement terms:  ", agg(st[!iv, ]), "\n  intervention terms: ", agg(st[iv, ]), "\n")
cat("\nrun:", run$path, "\n")
