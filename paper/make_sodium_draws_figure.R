# paper/make_sodium_draws_figure.R --------------------------------------------------
# Results, laboratory propensities and the draw count (author, 2026-11-04). The stored final curve of the sodium
# low-side propensity, with every (draws, low draws) group of MIMIC-IV training stays drawn as a point at its
# propensity: observed log-odds of death relative to the reference group (never low over two draws, the group at
# the curve's anchor), sized by stays, shaded by draw count. Shows that the curve's troughs sit on large groups
# drawn once or twice and its peaks where only groups drawn four or more times can land.
#
# Hard rule 1: training rows are read locally and only group aggregates leave R. Groups are shown only with
# >= 50 stays, >= 5 deaths and >= 5 survivors; the group table is written to paper/figs/coefs/sodium_draw_groups.csv.
# Points are unadjusted; the curve is the model's term with the channel's other terms held fixed.
#
#   Rscript paper/make_sodium_draws_figure.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(ggplot2) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
cfg <- tar_read(cfg); tabs <- tar_read(tabs); priors <- tar_read(priors); tr <- tar_read(train_ids)
if (is.data.frame(tr)) tr <- tr$stay_id
OUT <- "paper/figs"; SG <- "sodium"; TERM <- "s(pi_minus)"; KEY <- "sodium/meas"

d  <- signal_frame(SG, "meas", tabs, cfg, priors_for(priors, SG, "final"), stay_ids = tr)
yv <- all.vars(attr(d, "formula"))[1]
sf <- tabs$signal_features[tabs$signal_features$signal == SG, ]; sf <- sf[match(d$stay_id, sf$stay_id), ]
g  <- data.frame(y = d[[yv]], pi = d$pi_minus, n = sf$n_obs, k = sf$k_low)
grp <- aggregate(cbind(stays = 1, deaths = y) ~ n + k, data = g, FUN = sum)
grp$pi <- vapply(seq_len(nrow(grp)), function(i) { v <- unique(round(g$pi[g$n == grp$n[i] & g$k == grp$k[i]], 12))
  if (length(v) != 1L) stop("group n=", grp$n[i], " k=", grp$k[i], " has ", length(v), " propensity values", call. = FALSE); v }, 0)
grp <- grp[grp$stays >= 50 & grp$deaths >= 5 & grp$stays - grp$deaths >= 5, ]
lo <- function(p) log(p / (1 - p))
ref <- grp[grp$n == 2 & grp$k == 0, ]; stopifnot(nrow(ref) == 1L)
grp$rel <- lo(grp$deaths / grp$stays) - lo(ref$deaths / ref$stays)
grp$se  <- sqrt(1 / grp$deaths + 1 / (grp$stays - grp$deaths))
grp$mortality <- grp$deaths / grp$stays
cat(sprintf("%d groups shown, covering %.3f of %d measured training stays\n", nrow(grp), sum(grp$stays) / nrow(g), nrow(g)))

cv <- readRDS("paper/figs/coefs/anchored_curves.rds")
z  <- cv[cv$key == KEY & cv$term == TERM, ]; z <- z[order(z$grid_x), ]
ANCH <- read.csv("paper/figs/coefs/anchors.csv", stringsAsFactors = FALSE)
x0 <- ANCH$anchor[ANCH$key == KEY & ANCH$term == TERM][1]
if (abs(x0 - ref$pi) > 1e-3) stop("curve anchor ", x0, " is not the reference group's propensity ", ref$pi, call. = FALSE)
grp$multiplier <- exp(approx(z$grid_x, z$fit, grp$pi)$y)
write.csv(grp[order(grp$pi), c("pi", "n", "k", "stays", "mortality", "multiplier")], "paper/figs/coefs/sodium_draw_groups.csv", row.names = FALSE)

# Draw count is ordered magnitude: one sequential hue, light to dark (dataviz reference ramp, steps 250, 450, 650),
# in three bands so the troughs (drawn once or twice) and the peaks (four or more) separate at a glance.
SEQ <- c("1 or 2" = "#86b6ef", "3" = "#2a78d6", "4 or more" = "#104281")
grp$draws <- factor(ifelse(grp$n <= 2, "1 or 2", ifelse(grp$n == 3, "3", "4 or more")), names(SEQ))
grp$lab <- ifelse((grp$n == 1 & grp$k %in% 0:1) | (grp$n == 2 & grp$k == 1) | (grp$n == 4 & grp$k == 1) | (grp$n == 5 & grp$k == 0),
                  sprintf("%d of %d low", grp$k, grp$n), NA)
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e1e0d9"
zs <- z[z$supported, ]
p <- ggplot() +
  geom_hline(yintercept = 0, colour = "#c3c2b7", linewidth = 0.4) +
  geom_ribbon(data = zs, aes(grid_x, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), fill = INK2, alpha = 0.12) +
  geom_line(data = zs, aes(grid_x, fit), colour = INK2, linewidth = 0.7) +
  geom_linerange(data = grp, aes(pi, ymin = rel - 1.96 * se, ymax = rel + 1.96 * se, colour = draws), linewidth = 0.4, alpha = 0.4) +
  geom_point(data = grp, aes(pi, rel, size = stays, fill = draws), shape = 21, colour = "#fcfcfb", stroke = 0.6) +
  ggrepel::geom_label_repel(data = grp[!is.na(grp$lab), ], aes(pi, rel, label = lab), size = 2.9, colour = INK2, seed = 3,
                            fill = "#ffffff", label.size = 0, label.padding = 0.12, min.segment.length = 0,
                            segment.colour = "#898781", box.padding = 0.6, point.padding = 0.5, force = 3) +
  scale_fill_manual(values = SEQ, name = "Sodium draws") + scale_colour_manual(values = SEQ, guide = "none") +
  scale_size(range = c(2, 10), breaks = c(100, 1000, 10000), labels = c("100", "1,000", "10,000"), name = "Stays") +
  guides(fill = guide_legend(override.aes = list(size = 4), order = 1),
         size = guide_legend(override.aes = list(fill = "#898781", colour = "#fcfcfb"), order = 2)) +
  labs(x = "Propensity of sodium to fall below the reference range",
       y = "Log-odds of death relative to the reference (nats)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2),
        legend.title = element_text(colour = INK2, size = 9), legend.text = element_text(colour = INK2),
        plot.background = element_rect(fill = "white", colour = NA))
ragg::agg_png(file.path(OUT, "F29_sodium_draws.png"), width = 7.2, height = 4.4, units = "in", res = 300); print(p); dev.off()
grDevices::cairo_pdf(file.path(OUT, "F29_sodium_draws.pdf"), width = 7.2, height = 4.4); print(p); dev.off()
cat("wrote F29_sodium_draws\n")
