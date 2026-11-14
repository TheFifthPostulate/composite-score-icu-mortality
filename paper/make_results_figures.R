# paper/make_results_figures.R -----------------------------------------------------
# Results figures, 2026-10-31. Dropped by the author: the stratum-resolution and bold-count figures, and the
# difference version of F23.
#   F22_discrimination          AUROC and AUPRC lift per model, MIMIC-IV OOF / test / eICU, patient-bootstrap 95% CIs;
#                               severity block on the stays with severity scores (test, eICU)
#   F23_strata_transport        quartile mortality per model: observed (test, eICU) against the reported MIMIC-IV reference
#   F26_design_ladder           largest quartile difference (test, eICU), same-quartile share and bold-channel Jaccard
#                               between retrained models; dots full-train (gap) or median over pairs, bars 2.5-97.5 percentiles
#   F27_channel_calibration     per-channel calibration slope, MIMIC-IV OOF and eICU, profile-likelihood 95% CIs
# Reads aggregate run tables only (hard rule 1).
#
#   Rscript paper/make_results_figures.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(scales) })
R <- "out/runs"; OUT <- "paper/figs"
# Input runs (options added 2026-11-13): --internal --test --external --jgamgated --jgamgatedtest --jgamungated
# --jgamungatedtest --ebmadd --ebmaddtest --ebm --ebmtest --strata --stabstrata --stabmetrics --chansem, each a
# run directory under out/runs. The defaults are the runs of record.
args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
.run <- function(flag, default) {
  v <- .opt(flag, NULL); if (is.null(v)) return(default)
  if (!identical(normalizePath(dirname(v), winslash = "/", mustWork = TRUE), normalizePath(R, winslash = "/", mustWork = TRUE)))
    stop(flag, " must name a run directory under ", R, call. = FALSE)
  basename(v)
}
RN <- list(internal = .run("--internal", "internal_20260909T112643"), test = .run("--test", "test_20260912T215546"),
           external = .run("--external", "external_20260909T112824"),
           jgamgated = .run("--jgamgated", "jgamgated_20261030T151327"), jgamgatedtest = .run("--jgamgatedtest", "jgamgatedtest_20261030T151338"),
           jgamungated = .run("--jgamungated", "jgamungated_20261030T150052"), jgamungatedtest = .run("--jgamungatedtest", "jgamungatedtest_20261030T150101"),
           ebmadd = .run("--ebmadd", "ebmadd_20261028T234749"), ebmaddtest = .run("--ebmaddtest", "ebmaddtest_20261028T234803"),
           ebm = .run("--ebm", "ebm_20261028T224038"), ebmtest = .run("--ebmtest", "ebmtest_20261028T224126"),
           strata = .run("--strata", "strata_20261030T151638"), stabstrata = .run("--stabstrata", "stabstrata_20261030T195855"),
           stabmetrics = .run("--stabmetrics", "stabmetrics_20261030T190010"), chansem = .run("--chansem", "chansem_20260914T084935"))
rd <- function(...) utils::read.csv(file.path(R, ...), stringsAsFactors = FALSE)
one <- function(d, what) { if (nrow(d) != 1L) stop(what, ": expected one row, found ", nrow(d), call. = FALSE); d }

INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e1"
SITE_COL <- c("MIMIC-IV (out of fold)" = "#2a78d6", "MIMIC-IV test" = "#1baf7a", "eICU" = "#eb6834")
MOD_COL <- c("Score" = "#1d4f9c", "Joint GAM, gated" = "#6f9fdc", "Joint GAM, ungated" = "#eda100",
             "EBM, additive" = "#eb6834", "XGBoost + SHAP" = "#4a3aa7")
thm <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2),
        strip.text = element_text(colour = INK, face = "bold", hjust = 0),
        legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(colour = INK, face = "bold", size = 10.5),
        plot.background = element_rect(fill = "white", colour = NA))
save_fig <- function(p, name, w, h) {
  ragg::agg_png(file.path(OUT, paste0(name, ".png")), width = w, height = h, units = "in", res = 300); print(p); dev.off()
  grDevices::cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h); print(p); dev.off()
  cat("wrote", name, "\n")
}
wilson <- function(k, n, z = 1.96) {
  p <- k / n; d <- 1 + z^2 / n; c0 <- (p + z^2 / (2 * n)) / d
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  data.frame(rate = p, lo = c0 - h, hi = c0 + h)
}
pct <- function(x) unname(stats::quantile(x, c(0.025, 0.975)))

# ==== F22: discrimination at the three sites ====================================
I <- file.path(RN$internal, "tables"); T <- file.path(RN$test, "tables"); E <- file.path(RN$external, "tables")
SITES <- names(SITE_COL)
from_summary <- function(file, label, site, model) {
  d <- one(rd(file)[rd(file)$label == label, ], label)
  data.frame(model = model, site = site, auroc = d$auroc, auroc_lo = d$auroc_lo, auroc_hi = d$auroc_hi,
             lift = d$auprc / d$event_rate, lift_lo = d$auprc_lo / d$event_rate, lift_hi = d$auprc_hi / d$event_rate)
}
from_metrics <- function(run, site_key, site, model) {
  d <- rd(run, "tables", "ebm_metrics.csv"); d <- one(d[d$site == site_key, ], paste(run, site_key))
  br <- d$deaths / d$n
  data.frame(model = model, site = site, auroc = d$auroc, auroc_lo = d$auroc_lo, auroc_hi = d$auroc_hi,
             lift = d$auprc / br, lift_lo = d$auprc_lo / br, lift_hi = d$auprc_hi / br)
}
summ3 <- function(stem, model) rbind(
  from_summary(file.path(I, "score_summary_oof.csv"), paste0(stem, "_oof"), SITES[1], model),
  from_summary(file.path(T, "score_summary.csv"), paste0(stem, "_frozen"), SITES[2], model),
  from_summary(file.path(E, "score_summary.csv"), paste0(stem, "_frozen"), SITES[3], model))
met3 <- function(run, run_test, model) rbind(
  from_metrics(run, "mimic_oof", SITES[1], model), from_metrics(run_test, "mimic_test", SITES[2], model),
  from_metrics(run, "eicu", SITES[3], model))
DISC_LEV <- c("Score", "Score, measurement terms only", "Score, within-channel interactions",
              "Joint GAM, gated", "Joint GAM, ungated", "EBM, additive", "EBM, pairwise interactions",
              "XGBoost, constructed variables", "XGBoost, raw summaries")
disc <- rbind(summ3("llr_sum", DISC_LEV[1]), summ3("llr_meas", DISC_LEV[2]), summ3("llr_sum_ti_all", DISC_LEV[3]),
              met3(RN$jgamgated, RN$jgamgatedtest, DISC_LEV[4]),
              met3(RN$jgamungated, RN$jgamungatedtest, DISC_LEV[5]),
              met3(RN$ebmadd, RN$ebmaddtest, DISC_LEV[6]),
              met3(RN$ebm, RN$ebmtest, DISC_LEV[7]),
              summ3("xgb_feat", DISC_LEV[8]), summ3("xgb_raw", DISC_LEV[9]))
SEV_LEV <- c("Score", "APACHE II, acute physiology score", "APACHE II, total", "SOFA")
sev <- do.call(rbind, lapply(list(c("llr_sum_sev", SEV_LEV[1]), c("apache2_aps_sev", SEV_LEV[2]),
                                  c("apache2_total_sev", SEV_LEV[3]), c("sofa_sev", SEV_LEV[4])), function(x) rbind(
  from_summary(file.path(T, "severity_score_summary.csv"), x[1], SITES[2], x[2]),
  from_summary(file.path(E, "severity_score_summary.csv"), x[1], SITES[3], x[2]))))
long <- function(d, lev) {
  out <- rbind(data.frame(d[, c("model", "site")], metric = "AUROC", v = d$auroc, lo = d$auroc_lo, hi = d$auroc_hi),
               data.frame(d[, c("model", "site")], metric = "AUPRC lift (AUPRC / mortality rate)", v = d$lift, lo = d$lift_lo, hi = d$lift_hi))
  out$model <- factor(out$model, rev(lev)); out$site <- factor(out$site, SITES)
  out$metric <- factor(out$metric, c("AUROC", "AUPRC lift (AUPRC / mortality rate)")); out
}
dl <- long(disc, DISC_LEV); sl <- long(sev, SEV_LEV)
print(disc, digits = 3, row.names = FALSE)
disc_plot <- function(d, title) ggplot(d, aes(v, model)) +
  geom_line(aes(group = model), colour = GRID, linewidth = 1.6) +
  geom_errorbar(aes(xmin = lo, xmax = hi, colour = site), width = 0, linewidth = 0.5, orientation = "y") +
  geom_point(aes(colour = site), size = 2) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = SITE_COL, drop = FALSE) +
  labs(title = title, x = NULL, y = NULL) + thm + theme(panel.spacing.x = unit(1.5, "lines"))
p22 <- disc_plot(dl, "A. All stays") / (disc_plot(sl, "B. Stays with severity scores") + theme(strip.text = element_blank(), legend.position = "none")) +
  plot_layout(heights = c(9, 4))
save_fig(p22, "F22_discrimination", 7.5, 5.4)

# ==== F23: quartile transport, every model ======================================
ARM5 <- c(llr_sum = "Score", jgamgated = "Joint GAM, gated", jgamungated = "Joint GAM, ungated",
          ebmadd = "EBM, additive", xgb_feat = "XGBoost + SHAP")
sc <- rd(RN$strata, "tables", "strata_calibration.csv")
sc <- sc[sc$partition == "quartile" & sc$arm %in% names(ARM5), ]
if (nrow(sc) != 2L * 4L * length(ARM5)) stop("expected four quartiles per model at two sites", call. = FALSE)
sc$model <- factor(ARM5[sc$arm], ARM5)
sc$site <- factor(c(mimic_test = "MIMIC-IV test", eicu = "eICU")[sc$site], c("MIMIC-IV test", "eICU"))
L_COL2 <- SITE_COL[c("MIMIC-IV test", "eICU")]
sc$x <- sc$stratum + ifelse(sc$site == "eICU", 0.15, -0.15)
ref <- unique(sc[, c("model", "stratum", "reference")])
# mortality per quartile (author's choice, 2026-10-31), the reported (reference) rate as a black bar, observed rates as dots
p23 <- ggplot(sc) +
  geom_segment(data = ref, aes(x = stratum - 0.38, xend = stratum + 0.38, y = 100 * reference, yend = 100 * reference),
               colour = INK, linewidth = 0.9) +
  geom_errorbar(aes(x = x, ymin = 100 * observed_lo, ymax = 100 * observed_hi, colour = site), width = 0, linewidth = 0.5) +
  geom_point(aes(x = x, y = 100 * observed, colour = site), size = 1.9) +
  facet_wrap(~ model, nrow = 1) +
  scale_x_continuous(breaks = 1:4, labels = paste0("Q", 1:4)) +
  scale_colour_manual(values = L_COL2) +
  labs(x = NULL, y = "Mortality (%)", caption = "Black bar: reported (MIMIC-IV reference) mortality of the quartile. Dots: observed, Wilson 95% CI.") +
  thm + theme(panel.grid.major.x = element_blank(), plot.caption = element_text(colour = INK2, size = 8, hjust = 0))
save_fig(p23, "F23_strata_transport", 7.5, 3.4)
# ==== F26: design ladder, revised ================================================
STAB <- c(score = "Score", jgam_gated = "Joint GAM, gated", jgam_ungated = "Joint GAM, ungated", shap = "XGBoost + SHAP")
L_SITE <- c(test = "MIMIC-IV test", eicu = "eICU"); L_COL <- SITE_COL[c("MIMIC-IV test", "eICU")]
PAN <- c("A. Largest quartile difference (points)", "B. Stays in the same quartile, two retrained models", "C. Jaccard index of the bold channels")
cs <- rd(RN$strata, "tables", "strata_calibration_summary.csv")
cs <- cs[cs$partition == "quartile" & cs$arm %in% names(ARM5) & cs$site %in% c("mimic_test", "eicu"), ]
gap <- data.frame(model = ARM5[cs$arm], site = L_SITE[c(mimic_test = "test", eicu = "eicu")[cs$site]],
                  panel = PAN[1], v = 100 * cs$max_abs_difference, lo = NA_real_, hi = NA_real_)
ss <- rd(RN$stabstrata, "tables", "strata_stability.csv")
ss <- ss[ss$partition == "quartile" & ss$metric == "max_abs_gap", ]
for (i in seq_len(nrow(ss))) {
  j <- gap$model == STAB[[ss$arm[i]]] & gap$site == L_SITE[[ss$site[i]]]
  if (sum(j) != 1L) stop("ladder: no full-train row for ", ss$arm[i], " ", ss$site[i], call. = FALSE)
  gap[j, c("lo", "hi")] <- c(ss$p025[i], ss$p975[i])
}
aa <- rd(RN$stabstrata, "tables", "assignment_agreement.csv"); aa <- aa[aa$metric == "same_quartile", ]
same <- data.frame(model = STAB[aa$arm], site = L_SITE[aa$site], panel = PAN[2], v = 100 * aa$median, lo = 100 * aa$p025, hi = 100 * aa$p975)
pr <- rd(RN$stabmetrics, "tables", "agreement_pairs.csv"); pr <- pr[pr$contrast == "bag", ]
jac <- do.call(rbind, lapply(split(pr, list(pr$arm, pr$site)), function(d) {
  q <- pct(d$jaccard_ch); data.frame(model = STAB[[d$arm[1]]], site = L_SITE[[d$site[1]]], panel = PAN[3],
                                     v = stats::median(d$jaccard_ch), lo = q[1], hi = q[2]) }))
lad <- rbind(gap, same, jac)
lad$model <- factor(lad$model, rev(unname(ARM5))); lad$site <- factor(lad$site, names(L_COL)); lad$panel <- factor(lad$panel, PAN)
print(lad[order(lad$panel, lad$model, lad$site), ], digits = 3, row.names = FALSE)
dg <- position_dodge(width = 0.5)
p26 <- ggplot(lad, aes(v, model, colour = site)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.8, alpha = 0.55, orientation = "y", position = dg, na.rm = TRUE) +
  geom_point(size = 2.1, position = dg) +
  facet_wrap(~ panel, nrow = 1, scales = "free_x", labeller = label_wrap_gen(30)) +
  scale_colour_manual(values = L_COL) +
  labs(x = NULL, y = NULL) + thm + theme(panel.spacing.x = unit(1.3, "lines"))
save_fig(p26, "F26_design_ladder", 7.5, 3.2)

# ==== F27: channel calibration ==================================================
CH <- c(mbp = "Mean blood pressure", heart_rate = "Heart rate", spo2 = "Oxygen saturation", resp_rate = "Respiratory rate",
        glucose = "Glucose", gcs_motor = "Glasgow motor", gcs_eyes = "Glasgow eyes", gcs_verbal = "Glasgow verbal",
        urine_output_rate = "Urine output rate", creatinine = "Creatinine", platelet = "Platelets", hemoglobin = "Hemoglobin",
        temperature = "Temperature", sodium = "Sodium", bicarbonate = "Bicarbonate", bun = "Urea nitrogen",
        wbc = "White cell count", lactate = "Lactate", bilirubin_total = "Bilirubin")
cc <- rd(RN$chansem, "diagnostics", "channel_calibration.csv")
if (nrow(cc) != 38L) stop("expected 19 channels at two sites", call. = FALSE)
cc$channel <- factor(CH[cc$signal], rev(CH))
cc$site <- factor(c(mimic_oof = "MIMIC-IV (out of fold)", eicu = "eICU")[cc$site], c("MIMIC-IV (out of fold)", "eICU"))
p27 <- ggplot(cc, aes(slope, channel, colour = site)) +
  geom_vline(xintercept = 1, colour = INK2, linewidth = 0.4) +
  geom_errorbar(aes(xmin = slope_lo, xmax = slope_hi), width = 0, linewidth = 0.5, orientation = "y", position = position_dodge(width = 0.5)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = SITE_COL[c("MIMIC-IV (out of fold)", "eICU")]) +
  labs(x = "Calibration slope of the channel's weight of evidence (95% CI)", y = NULL) + thm
save_fig(p27, "F27_channel_calibration", 6, 5.2)
