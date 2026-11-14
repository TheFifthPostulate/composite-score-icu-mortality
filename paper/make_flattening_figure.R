# paper/make_flattening_figure.R ---------------------------------------------------
# Results 4.3, the joint curves (author, 2026-10-31). Flattening in the gated joint GAM does not follow redundancy.
#   A  per channel: R^2 of its stand-alone weight of evidence on the other 18 (flexible, tests/jointgam_redundancy.R)
#      against the median share of the stand-alone range that its joint measurement curves keep
#   B  Glasgow motor, lowest value: a redundant channel, flattened as expected
#   C  temperature, fever propensity: the least redundant channel, its curve gone in the joint fit
# Curves as in tests/jointgam_atlas.R: the paper's stored anchored stand-alone curve against the gated joint GAM's
# final fit, both anchored at the paper's reference value, 95% bands. Asserted: the redrawn swing ratios equal the
# stored ones (jgamshape run, AC_swing_ratio) to 1e-3 relative: the stored table evaluates the bundle model, the paper
# curve grid differs from it by about 3e-5 relative. Reads model parameters and aggregate grids only (hard rule 1).
#
#   Rscript paper/make_flattening_figure.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(mgcv); library(ggplot2); library(patchwork) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
.want <- c("CH", "UNIT", "IV", "describe", "var_of", "wrap_title")
for (e in parse("paper/make_atlas.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% .want) eval(e, globalenv())
}
cfg  <- tar_read(cfg)
fo   <- readRDS("out/cache/jointgam/gated/fits.rds")$final
cv   <- readRDS("paper/figs/coefs/anchored_curves.rds")
ANCH <- read.csv("paper/figs/coefs/anchors.csv", stringsAsFactors = FALSE)
# --jgamredund --jgamshape <run>: options added 2026-11-13; defaults are the runs of record.
args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
RED  <- read.csv(file.path(.opt("--jgamredund", "out/runs/jgamredund_20261030T145610"), "tables/redundancy_vs_flattening.csv"), stringsAsFactors = FALSE)
SHP  <- read.csv(file.path(.opt("--jgamshape",  "out/runs/jgamshape_20261030T142353"),  "tables/shape_comparison.csv"), stringsAsFactors = FALSE)
OUT  <- "paper/figs"

POS <- list(); at <- 1L
for (bb in fo$blocks) { POS[[bb$name]] <- at + seq_len(bb$ncol); at <- at + bb$ncol }
BLK <- stats::setNames(fo$blocks, vapply(fo$blocks, `[[`, "", "name"))
joint_curve <- function(col, x, x0) {
  bb <- BLK[[col]]; if (is.null(bb) || bb$type != "smooth") stop("no smooth block ", col, call. = FALSE)
  D <- sweep(mgcv::PredictMat(bb$sm, data.frame(v = x)), 2L, as.numeric(mgcv::PredictMat(bb$sm, data.frame(v = x0)))) %*% bb$Z
  i <- POS[[col]]
  list(fit = as.numeric(D %*% fo$coef[i]), se = sqrt(pmax(rowSums((D %*% fo$Vp[i, i, drop = FALSE]) * D), 0)))
}

INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e1"
LAB_A <- "Channel model (stand-alone)"; LAB_C <- "Gated joint GAM"
COL <- stats::setNames(c("#2a78d6", "#c9722f"), c(LAB_A, LAB_C))
thm <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2),
        plot.title = element_text(colour = INK, face = "bold", size = 10),
        legend.position = "bottom", legend.title = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA))

curve_panel <- function(sg, term, title) {
  key <- if (is.na(spec_source(sg, "full", cfg))) paste0(sg, "/full") else paste0(sg, "/meas")
  z <- cv[cv$key == key & cv$term == term, ]; z <- z[order(z$grid_x), ]
  if (!nrow(z)) stop("no stored curve for ", key, " ", term, call. = FALSE)
  x0 <- ANCH$anchor[ANCH$key == key & ANCH$term == term][1]
  jc <- joint_curve(paste0(sg, "__", var_of(term)), z$grid_x, x0)
  ratio <- diff(range(jc$fit[z$supported])) / diff(range(z$fit[z$supported]))
  ref <- SHP$AC_swing_ratio[SHP$channel == sg & SHP$term == var_of(term)]
  if (length(ref) != 1L || abs(ratio - ref) > 1e-3 * max(ref, 1e-3)) stop("swing ratio of ", sg, " ", term, " = ", ratio, ", stored ", ref, call. = FALSE)
  cat(sprintf("%s %s: joint curve keeps %.3g of the stand-alone range\n", sg, term, ratio))
  zz <- z[z$supported, ]; jj <- lapply(jc, `[`, z$supported)
  d <- rbind(data.frame(x = zz$grid_x, fit = zz$fit, lo = zz$fit - 1.96 * zz$se, hi = zz$fit + 1.96 * zz$se, m = LAB_A),
             data.frame(x = zz$grid_x, fit = jj$fit, lo = jj$fit - 1.96 * jj$se, hi = jj$fit + 1.96 * jj$se, m = LAB_C))
  d$m <- factor(d$m, c(LAB_A, LAB_C))
  ggplot(d, aes(x, fit, colour = m, fill = m)) +
    geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.15) + geom_line(linewidth = 0.8) +
    scale_colour_manual(values = COL) + scale_fill_manual(values = COL) +
    labs(title = title, x = describe(sg, term)$x, y = "Weight of evidence (nats)") + thm
}

ok <- is.finite(RED$swing_ratio_median)
rho <- stats::cor(RED$r2_flexible[ok], RED$swing_ratio_median[ok], method = "spearman")
cat(sprintf("Spearman, R^2 (flexible) vs median swing ratio: %.2f over %d channels\n", rho, sum(ok)))
RED$lab <- CH[RED$channel]
RED$hl <- ifelse(grepl("^gcs_", RED$channel), "Glasgow", ifelse(RED$channel == "temperature", "Temperature", "other"))
pA <- ggplot(RED[ok, ], aes(r2_flexible, swing_ratio_median)) +
  geom_point(aes(colour = hl), size = 2.2) +
  ggrepel::geom_text_repel(aes(label = lab, colour = hl), size = 2.6, seed = 1, min.segment.length = 0.2, segment.colour = GRID, box.padding = 0.3) +
  scale_colour_manual(values = c(Glasgow = "#2a78d6", Temperature = "#c9722f", other = INK2), guide = "none") +
  scale_x_continuous(limits = c(0, 0.8)) + scale_y_continuous(limits = c(0, 0.75)) +
  labs(title = "A. Flattening against overlap",
       x = "Overlap with the other 18 channels (R²)",
       y = "Fraction of range kept by the joint curves\n(median over terms)") + thm
pB <- curve_panel("gcs_motor", "s(value_min_delta)", "B. Glasgow motor (overlap 0.73)")
pC <- curve_panel("temperature", "s(pi_plus)", "C. Temperature (overlap 0.11)")
p <- pA + (pB / pC + plot_layout(guides = "collect") & theme(legend.position = "bottom")) + plot_layout(widths = c(1.15, 1))
ragg::agg_png(file.path(OUT, "F28_joint_flattening.png"), width = 7.5, height = 4.6, units = "in", res = 300); print(p); dev.off()
grDevices::cairo_pdf(file.path(OUT, "F28_joint_flattening.pdf"), width = 7.5, height = 4.6); print(p); dev.off()
cat("wrote F28_joint_flattening\n")
