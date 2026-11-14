# paper/make_ladder_figure.R -------------------------------------------------------
# The design ladder (Results 4.3): what each design choice costs and buys.
#   F20_design_ladder (main text, with tab:ladder): one row per model, three panels at eICU,
#     each a dot and the 2.5-97.5 percentile range over the 38 retrained subsamples:
#       AUROC                        dot = full-train model; range over retrained models
#       largest quartile difference  dot = full-train model; range over retrained models
#       Jaccard of the bold channels dot = median over the 703 pairs of retrained models; range over pairs
#     The EBM was not retrained: full-train dot only, no Jaccard.
#   F20_ladder_deciles (appendix): eICU observed minus MIMIC-IV reference mortality per decile of
#     each model's own reference distribution (frozen cut points), full-train models.
# Reads aggregate run tables only.
#
#   Rscript paper/make_ladder_figure.R [strata run] [stabmetrics run] [stabstrata run]
# 2026-10-31: panel B (boxplots) replaced by the dot-and-interval figure; decile panel moved to the appendix.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
a <- commandArgs(trailingOnly = TRUE)
STRATA <- if (length(a) >= 1) a[1] else "out/runs/strata_20261030T151638"
STAB   <- if (length(a) >= 2) a[2] else "out/runs/stabmetrics_20261030T190010"
SSTR   <- if (length(a) >= 3) a[3] else "out/runs/stabstrata_20261030T195855"
OUT <- "paper/figs"

ARMS <- c(llr_sum = "Score", jgamgated = "Joint GAM, gated", jgamungated = "Joint GAM, ungated",
          ebmadd = "EBM, additive", xgb_feat = "XGBoost + SHAP")
STAB_ARM <- c(score = "Score", jgam_gated = "Joint GAM, gated", jgam_ungated = "Joint GAM, ungated", shap = "XGBoost + SHAP")
COL <- c("Score" = "#1d4f9c", "Joint GAM, gated" = "#6f9fdc", "Joint GAM, ungated" = "#eda100",
         "EBM, additive" = "#eb6834", "XGBoost + SHAP" = "#4a3aa7")
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e1"
thm <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2),
        strip.text = element_text(colour = INK, face = "bold", hjust = 0),
        legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(colour = INK, face = "bold", size = 10.5),
        plot.background = element_rect(fill = "white", colour = NA))
rd <- function(run, name) utils::read.csv(file.path(run, "tables", paste0(name, ".csv")), stringsAsFactors = FALSE)
pct <- function(x) unname(stats::quantile(x, c(0.025, 0.975)))
save_fig <- function(p, name, w, h) {
  ragg::agg_png(file.path(OUT, paste0(name, ".png")), width = w, height = h, units = "in", res = 300); print(p); dev.off()
  grDevices::cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h); print(p); dev.off()
  cat("wrote", name, "\n")
}

# --- F20: dot and retraining range per model, at eICU -----------------------------------
PANEL <- c(auroc = "AUROC", gap = "Largest quartile difference (points)", jac = "Jaccard index of the bold channels")
disc <- rd(STRATA, "strata_discrimination")
disc <- disc[disc$site == "eicu" & disc$assignment == "raw score" & disc$arm %in% names(ARMS), ]
cal <- rd(STRATA, "strata_calibration_summary")
cal <- cal[cal$site == "eicu" & cal$partition == "quartile" & cal$arm %in% names(ARMS), ]
if (nrow(disc) != length(ARMS) || nrow(cal) != length(ARMS)) stop("expected one eICU row per model", call. = FALSE)
pts <- rbind(data.frame(model = ARMS[disc$arm], panel = "auroc", value = disc$auroc),
             data.frame(model = ARMS[cal$arm], panel = "gap", value = 100 * cal$max_abs_difference))
pts$lo <- NA_real_; pts$hi <- NA_real_

ab <- rd(SSTR, "auroc_by_bag"); ab <- ab[ab$site == "eicu", ]
ss <- rd(SSTR, "strata_stability")
ss <- ss[ss$site == "eicu" & ss$partition == "quartile" & ss$metric == "max_abs_gap", ]
pr <- rd(STAB, "agreement_pairs"); pr <- pr[pr$contrast == "bag" & pr$site == "eicu", ]
n_bag <- table(ab$arm); n_pair <- table(pr$arm)
if (length(unique(as.vector(n_bag))) != 1L || length(unique(as.vector(n_pair))) != 1L)
  stop("the retrained models do not share one set of subsamples", call. = FALSE)
for (k in names(STAB_ARM)) {
  m <- STAB_ARM[[k]]
  r <- pct(ab$auroc[ab$arm == k]); pts[pts$model == m & pts$panel == "auroc", c("lo", "hi")] <- r
  s <- ss[ss$arm == k, ]; if (nrow(s) != 1L) stop("strata_stability: one row expected for ", k, call. = FALSE)
  pts[pts$model == m & pts$panel == "gap", c("lo", "hi")] <- c(s$p025, s$p975)
  j <- pr$jaccard_ch[pr$arm == k]
  pts <- rbind(pts, data.frame(model = m, panel = "jac", value = stats::median(j), lo = pct(j)[1], hi = pct(j)[2]))
}
pts$model <- factor(pts$model, rev(unname(ARMS)))
pts$panel <- factor(PANEL[pts$panel], PANEL)
print(pts[order(pts$panel, pts$model), ], row.names = FALSE, digits = 3)

p20 <- ggplot(pts, aes(value, model, colour = model)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.9, alpha = 0.55, na.rm = TRUE) +
  geom_point(size = 2.4) +
  facet_wrap(~ panel, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = COL, guide = "none") +
  labs(x = NULL, y = NULL) +
  thm + theme(panel.spacing.x = unit(1.2, "lines"))
cat(sprintf("ranges over %d retrained subsamples, Jaccard over %d pairs
", as.integer(n_bag[1]), as.integer(n_pair[1])))
save_fig(p20, "F20_design_ladder", 7.5, 2.3)

# --- appendix: transport of the deciles -----------------------------------------------------
sc <- rd(STRATA, "strata_calibration")
sc <- sc[sc$site == "eicu" & sc$partition == "decile" & sc$arm %in% names(ARMS), ]
if (nrow(sc) != 10L * length(ARMS)) stop("expected ten eICU deciles for each of ", length(ARMS), " arms", call. = FALSE)
sc$arm <- factor(ARMS[sc$arm], ARMS)
sc$gap <- 100 * (sc$observed - sc$reference)
pA <- ggplot(sc, aes(stratum, gap, colour = arm)) +
  geom_hline(yintercept = 0, colour = INK2, linewidth = 0.4) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
  scale_x_continuous(breaks = 1:10) +
  scale_colour_manual(values = COL, drop = FALSE) +
  labs(x = "Decile of the model's MIMIC-IV reference distribution",
       y = "eICU observed minus MIMIC-IV\nreference mortality (points)") + thm
save_fig(pA, "F20_ladder_deciles", 7.5, 3.6)
