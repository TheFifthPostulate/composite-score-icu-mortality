# paper/make_strata_figure.R -------------------------------------------------------
# Results 4.1: observed mortality of the score in its 20 frozen MIMIC-IV risk bins, at the three sites
# (MIMIC-IV out of fold, MIMIC-IV test, eICU), with Wilson 95% intervals. Light vertical lines mark the
# quartile boundaries (bins 1-5, 6-10, 11-15, 16-20). Same per-bin tables as F1 in make_figures.R.
# Reads aggregate run tables only.
#
#   Rscript paper/make_strata_figure.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(scales) })
# --internal --test --external <run>: options added 2026-11-13; defaults are the runs of record.
args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
I <- file.path(.opt("--internal", "out/runs/internal_20260909T112643"), "tables")
T <- file.path(.opt("--test",     "out/runs/test_20260912T215546"), "tables")
E <- file.path(.opt("--external", "out/runs/external_20260909T112824"), "tables")
OUT <- "paper/figs"
rd <- function(dir, name) read.csv(file.path(dir, paste0(name, ".csv")), stringsAsFactors = FALSE)
SITE_COL <- c("MIMIC-IV (out of fold)" = "#2a78d6", "MIMIC-IV test" = "#1baf7a", "eICU" = "#eb6834")
SITE_LEV <- names(SITE_COL)
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e1"
thm <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2),
        legend.position = "bottom", legend.title = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA))

cols <- c("bin", "n", "deaths", "obs_rate", "lo", "hi")
bins <- rbind(cbind(site = SITE_LEV[1], rd(I, "risk_bins_llr_sum_oof")[, cols]),
              cbind(site = SITE_LEV[2], rd(T, "risk_bins_llr_sum_frozen")[, cols]),
              cbind(site = SITE_LEV[3], rd(E, "risk_bins_llr_sum_frozen")[, cols]))
if (!all(table(bins$site) == 20L)) stop("expected 20 bins per site", call. = FALSE)
bins$site <- factor(bins$site, SITE_LEV)
lab <- data.frame(x = c(3, 8, 13, 18), lab = paste0("Q", 1:4))
dg <- position_dodge(width = 0.55)
p <- ggplot(bins, aes(bin, obs_rate, colour = site)) +
  geom_vline(xintercept = c(5.5, 10.5, 15.5), colour = GRID, linewidth = 0.6) +
  geom_text(data = lab, aes(x, Inf, label = lab), inherit.aes = FALSE, vjust = 1.4, colour = INK2, size = 3.2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.4, alpha = 0.7, position = dg) +
  geom_line(linewidth = 0.6, position = dg) +
  geom_point(size = 1.6, position = dg) +
  scale_colour_manual(values = SITE_COL) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.1))) +
  labs(x = "Risk bin of the score on MIMIC-IV training cut points (20 equal-count bins)",
       y = "Observed in-hospital mortality") + thm
ragg::agg_png(file.path(OUT, "F21_score_strata.png"), width = 6.5, height = 3.6, units = "in", res = 300); print(p); dev.off()
grDevices::cairo_pdf(file.path(OUT, "F21_score_strata.pdf"), width = 6.5, height = 3.6); print(p); dev.off()
cat("wrote F21_score_strata\n")
