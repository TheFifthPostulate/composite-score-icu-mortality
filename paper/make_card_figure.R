# paper/make_card_figure.R -----------------------------------------------------
# Draw the patient-card exhibits from a card run's exported tables. The stay
# identifier is never read; the figures show derived evidence values only.
#   Rscript paper/make_card_figure.R out/runs/card_<id>
#
# A single-stay run (tables/per_signal.csv) gives F14_patient_card. A profile
# run (run/patient_card.R --profiles; tables/per_signal_<profile>.csv) gives
# one card per reference-risk profile, F14_patient_card_<profile>, and the term
# figures F17_card_terms_<signal>: the fitted terms of each decomposed signal's final model
# (the QC curve table, as F15) with each profile's stay marked on every term.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(scales) })
args <- commandArgs(trailingOnly = TRUE)
# --gamqc <run>: the GAM QC run (option added 2026-11-13; default is the run of record). The option and its
# value are taken out before the positional arguments are read.
.gi <- match("--gamqc", args)
GAMQC <- if (is.na(.gi) || .gi == length(args)) "out/runs/gamqc_20260909T180107" else args[.gi + 1L]
pos <- if (is.na(.gi)) args else args[-c(.gi, .gi + 1L)]
RUN <- if (length(pos)) pos[1] else stop("pass the card run directory")
# Optional second argument: the output directory (default paper/figs), e.g. the
# run's own directory for a review card that is not a paper figure.
OUT_DIR <- if (length(pos) >= 2L) pos[2] else "paper/figs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
G <- file.path(GAMQC, "tables")
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e1"; BLUE <- "#2a78d6"
thm <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2),
        plot.title = element_text(colour = INK, face = "bold", size = 11),
        plot.subtitle = element_text(colour = INK2, size = 9), legend.position = "none",
        plot.background = element_rect(fill = "white", colour = NA))
PROFILE_LAB <- c(low = "Low reference risk", intermediate = "Ambiguous reference risk", high = "High reference risk")
PROFILE_COL <- c("Low reference risk" = "#1baf7a", "Ambiguous reference risk" = "#eda100",
                 "High reference risk" = "#eb6834")
save_both <- function(p, name, w, h) {
  for (ext in c("png", "pdf")) {
    f <- file.path(OUT_DIR, paste0(name, ".", ext))
    if (ext == "png") ragg::agg_png(f, width = w, height = h, units = "in", res = 300) else grDevices::cairo_pdf(f, width = w, height = h)
    print(p); dev.off()
  }
  cat("wrote", name, "\n")
}

# Card rows follow the frozen reporting partition (config/domains.csv), grouped
# by physiological system in the display order of paper/domain_order.R.
source("paper/domain_order.R")

ref_path <- file.path(RUN, "tables", "signal_reference.csv")
has_ref <- file.exists(ref_path)
rq_all <- if (has_ref) read.csv(ref_path, stringsAsFactors = FALSE) else NULL

# --- one card -------------------------------------------------------------------
draw_card <- function(per, sm, title, ct = NULL) {
  # Rows in the domain display order, first domain at the top; each domain is a
  # facet row, so its channels sit together under one label.
  # Paper names; a dagger marks a channel with a term outside training support.
  stopifnot(setequal(per$signal, DOM$signal))
  per <- per[match(DOM$signal, per$signal), ]
  per$sig <- as.character(per$signal)
  per$domain <- factor(DOM_LAB[DOM$domain], unname(DOM_LAB[DOMAIN_ORDER]))
  lab <- unname(ifelse(per$signal %in% names(CH_LAB), CH_LAB[per$signal], per$signal))
  # Units in the channel label when the card carries the variable table, whose
  # tail and trend values are in the channel's units (Table tab:bounds).
  UNIT <- c(mbp = "mmHg", heart_rate = "/min", spo2 = "%", resp_rate = "/min", glucose = "mg/dL",
            urine_output_rate = "mL/kg/h", creatinine = "mg/dL", platelet = "10^9/L", hemoglobin = "g/dL",
            temperature = "\u00b0C", sodium = "mmol/L", bicarbonate = "mmol/L", bun = "mg/dL", wbc = "10^9/L",
            lactate = "mmol/L", bilirubin_total = "mg/dL")
  if (!is.null(ct)) lab <- ifelse(per$sig %in% names(UNIT), paste0(lab, " (", UNIT[per$sig], ")"), lab)
  flagged <- if ("n_terms_outside" %in% names(per)) per$n_terms_outside > 0 else rep(FALSE, nrow(per))
  lab[flagged] <- paste0(lab[flagged], " †")
  lab_of <- stats::setNames(lab, per$signal)
  dom_of <- stats::setNames(per$domain, per$signal)
  per$signal <- factor(lab, rev(lab))
  per$pct_train[!per$measured] <- NA
  per$kind <- ifelse(per$measured, "measured", "unmeasured (L = 0 by assignment)")
  # The channel's training distribution, drawn behind each point on the same
  # nats axis: a light strip for the 1st to 99th percentile and a darker one
  # for the 10th to 90th, with a tick at the median.
  STRIP <- "#d9d6cf"; STRIP2 <- "#bfbbb1"
  # Domains as facet rows: labeled at the left of the evidence panel (outside
  # the channel names), unlabeled in the percentile panel.
  dom_facet <- function(switch = NULL) list(
    facet_grid(rows = vars(domain), scales = "free_y", space = "free_y", switch = switch),
    theme(panel.spacing.y = unit(3, "pt")))
  # The nats axis is the log of the likelihood ratio, so a likelihood-ratio
  # axis is the same axis with other labels: nats below, likelihood ratio above.
  LR_BREAKS <- c(0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 50)
  lr_lab <- function(x) paste0(format(x, drop0trailing = TRUE, trim = TRUE), "\u00d7")
  pA <- ggplot(per, aes(L_full, signal))
  if (has_ref) {
    rq <- rq_all; rq$domain <- dom_of[rq$signal]
    rq$signal <- factor(unname(lab_of[rq$signal]), levels(per$signal))
    pA <- pA +
      geom_tile(data = rq, aes(x = (p01 + p99) / 2, y = signal, width = p99 - p01, height = 0.64),
                inherit.aes = FALSE, fill = STRIP, colour = NA) +
      geom_tile(data = rq, aes(x = (p10 + p90) / 2, y = signal, width = p90 - p10, height = 0.64),
                inherit.aes = FALSE, fill = STRIP2, colour = NA) +
      geom_tile(data = rq, aes(x = p50, y = signal, width = 0.012, height = 0.64),
                inherit.aes = FALSE, fill = "white", colour = NA)
  }
  pA <- pA +
    geom_vline(xintercept = 0, colour = INK2, linewidth = 0.3) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 1.2, colour = BLUE, alpha = 0.6, orientation = "y") +
    geom_point(aes(shape = kind), colour = BLUE, size = 2.4, fill = "white") +
    scale_shape_manual(values = c("measured" = 16, "unmeasured (L = 0 by assignment)" = 21)) +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Likelihood ratio of the channel's observations",
                                           breaks = log(LR_BREAKS), labels = lr_lab(LR_BREAKS))) +
    labs(x = "Evidence L (nats) with 95% posterior band", y = NULL,
         title = "Evidence per channel, by system",
         subtitle = if (has_ref) paste0("Grey strips: the channel's MIMIC-IV training distribution\n(1st-99th, 10th-90th percentile; white tick = median).\nOpen circles: unmeasured, L = 0.",
                                        if (any(flagged)) "\n† A term of this channel lies outside its training support;\n   the value is extrapolated." else "")
                    else "Open circles: unmeasured channels, L = 0 by assignment") + thm + dom_facet("y") +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1, colour = INK, face = "bold", size = 8.5))
  # Where the value and its band end points sit within the channel's MIMIC-IV
  # training distribution (measured stays only; unmeasured rows stay empty).
  pP <- ggplot(per, aes(pct_train, signal)) +
    geom_vline(xintercept = 50, colour = INK2, linewidth = 0.3, linetype = "22") +
    geom_errorbar(aes(xmin = pct_lo, xmax = pct_hi), width = 0, linewidth = 1.2, colour = "#c9722f", alpha = 0.55, orientation = "y", na.rm = TRUE) +
    geom_point(colour = "#c9722f", size = 2.4, na.rm = TRUE) +
    scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    labs(x = "Percentile within channel (MIMIC-IV training)", y = NULL,
         title = "Rarity",
         subtitle = "Percentile within\nits channel; band\nends mapped alike") +
    thm + dom_facet() + theme(axis.text.y = element_blank(), strip.text.y = element_blank())
  # THE VARIABLE TABLE (author, 2026-10-18; revised 2026-10-28,
  # docs/patient_card_changes_queue_20261028.md). The stay's random variables in
  # columns, one channel per row, with two lines where both sides of the
  # reference range enter the model. Columns: percentile, side, propensity,
  # deviated / observed hours, tail, trend, interventions.
  # A value's contribution to the channel's evidence is the anchored term of
  # R/15_anchors.R, relative to a quiet day on the channel (run/patient_card.R,
  # card_terms); an intervention's contribution sums all of its terms
  # (exposure, intensity, admission indicator). Every value carries its
  # likelihood multiplier ("×2.4"); a contribution of at least log 2 toward
  # death is set in bold, every other multiplier in grey (author, 2026-10-29).
  # The bold rule stays one-sided: strong death-ward evidence should surface in
  # any stratum. A fixed landmark, not a ranking. The percentile of the
  # channel's evidence within its channel is bold at 85 or above; the channel's
  # likelihood multiplier exp(L) follows it, grey below 85, and a grey second
  # line gives the quiet-day multiplier exp(reference evidence). The quiet-day
  # multiplier times the multipliers of the row's values equals the channel's,
  # which is asserted below. The masking state (run/patient_card.R, state_of()) is in
  # the header, not the table. Column borders are faint; rows
  # have none. The posterior intervals of the propensities stay in the
  # runner's table but are not drawn: with one or two observed hours they span
  # most of the unit interval, and the hour count already states how few hours
  # a value rests on.
  pT <- NULL
  if ("pi_minus" %in% names(per) && !is.null(ct)) {
    STRONG <- log(2)
    PCT_BOLD <- 85
    # A likelihood multiplier, as the evidence panel's top axis writes it, in
    # parentheses (author, 2026-10-28).
    # The format is chosen on the value as printed, so 9.96 reads "×10", not "×10.0".
    # Below 1 two decimals, so a survival-ward multiplier keeps two digits (×0.30).
    mf <- function(v) {
      m <- exp(v)
      paste0("(×", ifelse(round(m, 1) >= 10, sprintf("%.0f", m), ifelse(m < 0.995, sprintf("%.2f", m), sprintf("%.1f", m))), ")")
    }
    # The table is drawn as rich text, so only the multiplier of a value below
    # the bold rule is grey; a bold cell is wrapped in <b> by the layer.
    MULT_GREY <- "#8f8c85"
    gm <- function(v) paste0("<span style='color:", MULT_GREY, "'>", mf(v), "</span>")
    # A value followed by its multiplier: plain inside a bold cell, grey otherwise.
    sv <- function(txt, v) ifelse(is.na(v), txt, ifelse(v >= STRONG, paste(txt, mf(v)), paste(txt, gm(v))))
    # The contributions shown on each channel's row, for the reconciliation check.
    shown <- list()
    IV_LAB <- c(vasopressor = "vasopressor", inotrope = "inotrope", invasive_vent = "ventilation", fio2 = "FiO2",
                rrt = "RRT", sedation_benzo = "benzo", sedation_propofol = "propofol",
                sedation_dexmed = "dexmed", insulin = "insulin", diuretic = "diuretic",
                transfusion_prbc = "red cells", transfusion_platelet = "platelets")
    DIG <- c(temperature = 1, creatinine = 1, lactate = 1, bilirubin_total = 1, hemoglobin = 1,
             urine_output_rate = 2, wbc = 1)
    dig <- function(sg) if (sg %in% names(DIG)) DIG[[sg]] else 0
    pf  <- function(x) ifelse(100 * x < 10, sprintf("%.1f", 100 * x), sprintf("%.0f", 100 * x))
    contrib <- function(sg, pat) {
      z <- ct[ct$signal == sg & grepl(pat, ct$term), , drop = FALSE]
      if (nrow(z)) sum(z$contribution) else NA_real_
    }
    rows <- list()
    add <- function(sg, col, line, text, value = NA_real_, grey = FALSE, bold = isTRUE(value >= STRONG)) {
      rows[[length(rows) + 1L]] <<- data.frame(sig = sg, col = col, line = line, text = text,
                                               bold = bold, grey = grey, stringsAsFactors = FALSE)
      if (!is.na(value)) shown[[sg]] <<- c(shown[[sg]], value)
    }
    # The masking state is written once in the header (Masking State), not as a
    # column (author, 2026-10-28). A channel whose every charted hour was removed
    # reads "masked"; "not measured" only when nothing was charted.
    for (i in seq_len(nrow(per))) {
      r <- per[i, ]; sg <- r$sig
      if (!r$measured) { add(sg, "pct", 1, if (r$n_removed > 0) "masked" else "not measured", grey = TRUE); next }
      pr <- round(r$pct_train)
      add(sg, "pct", 1, paste0(sprintf("%.0f (%.0f-%.0f) ", r$pct_train, r$pct_lo, r$pct_hi),
                               if (pr >= PCT_BOLD) mf(r$L_full) else gm(r$L_full)), bold = pr >= PCT_BOLD)
      ref_ev <- unique(ct$reference_evidence[ct$signal == sg])
      if (length(ref_ev) != 1L) stop("expected one reference evidence for ", sg)
      add(sg, "pct", 2, paste("quiet day", gm(ref_ev)), bold = FALSE)
      shown[[sg]] <- c(shown[[sg]], ref_ev)
      sides <- c("low", "high")[c(!is.na(r$pi_minus), !is.na(r$pi_plus))]
      for (ln in seq_along(sides)) {
        sd <- sides[ln]
        p  <- if (sd == "low") r$pi_minus else r$pi_plus
        k  <- if (sd == "low") r$k_low else r$k_high
        cp <- contrib(sg, if (sd == "low") "^s\\(pi_minus\\)$" else "^s\\(pi_plus\\)$")
        add(sg, "side", ln, sd)
        add(sg, "prop", ln, paste0(sv(paste0(pf(p), "%"), cp), if (ln == 1 && isTRUE(r$heavy_shrinkage)) " ‡" else ""), cp)
        add(sg, "hours", ln, sprintf("%d/%d=%s%%", k, r$n_obs, pf(k / r$n_obs)))
        # The tail on this side, labeled by its statistic. A tail that the
        # excursion-side rule leaves out of the model is shown in grey with a
        # dash in place of a contribution (run/patient_card.R, stay_parameters()).
        STAT <- c(q05 = "p5", q95 = "p95", value_min = "min", value_max = "max")
        tvar <- if (sd == "low") r$tail_low_var else r$tail_high_var
        tst  <- if (sd == "low") r$tail_low_stat else r$tail_high_stat
        tv   <- if (sd == "low") r$tail_low else r$tail_high
        tmod <- isTRUE(if (sd == "low") r$tail_low_modeled else r$tail_high_modeled)
        # ONE DRAW: the minimum and the maximum are the same reading, so the
        # tail is shown once, as "singular" (author, 2026-10-28), with the summed contribution of the
        # modeled tail terms; how a single reading's evidence divides between
        # the min and max terms carries no meaning of its own.
        single <- isTRUE(r$n_obs == 1) && isTRUE(r$tail_low == r$tail_high) && length(sides) == 2L
        if (single) {
          if (ln == 1L) {
            tvars <- c(if (isTRUE(r$tail_low_modeled)) r$tail_low_var, if (isTRUE(r$tail_high_modeled)) r$tail_high_var)
            cs <- sum(vapply(tvars, function(v) contrib(sg, paste0("^s\\(", v, "\\)$")), numeric(1)))
            add(sg, "tail", ln, sv(paste("singular", formatC(r$tail_low, format = "f", digits = dig(sg))), cs), cs)
          }
        } else if (!is.na(tv)) {
          tnum <- paste(STAT[[tst]], formatC(tv, format = "f", digits = dig(sg)))
          if (tmod) {
            ctl <- contrib(sg, paste0("^s\\(", tvar, "\\)$"))
            add(sg, "tail", ln, sv(tnum, ctl), ctl)
          } else {
            add(sg, "tail", ln, paste(tnum, "—"), grey = TRUE)
          }
        }
      }
      # One line per shown intervention, so each carries its own bold; three
      # lines (the Glasgow components' three sedatives) take the three-line row
      # layout.
      n_iv <- 0L
      if (!is.na(r$exposures) && nzchar(r$exposures)) {
        it   <- lapply(strsplit(r$exposures, ";", fixed = TRUE)[[1]], function(e) strsplit(e, "|", fixed = TRUE)[[1]])
        ivs  <- vapply(it, `[`, "", 1); kind <- vapply(it, `[`, "", 2); val <- as.numeric(vapply(it, `[`, "", 3))
        adm  <- vapply(it, function(e) identical(e[4], "1"), logical(1))
        civ  <- vapply(ivs, function(iv) contrib(sg, paste0("(^|\\()", iv, "__")), numeric(1))
        nm   <- ifelse(ivs %in% names(IV_LAB), IV_LAB[ivs], ivs)
        vtxt <- ifelse(kind == "exposure_frac", sprintf("%.0f%%", 100 * val), sprintf("%.0f h", val))
        # An unexposed intervention sits at its anchor, so its contribution is
        # zero; it is left off (author, 2026-10-28).
        keep <- !(val == 0 & !adm & (is.na(civ) | abs(civ) < 0.005))
        txt  <- sv(paste0(nm, " ", vtxt, ifelse(adm, ", from admission", "")), civ)[keep]
        val2 <- civ[keep]
        if (length(txt) > 3L) stop("more than three interventions shown for ", sg)
        for (j in seq_along(txt)) add(sg, "iv", j, txt[j], val2[j])
        n_iv <- length(txt)
      }
      # "N/A" for an unpaired channel, which has no intervention (author, 2026-10-28);
      # a dash for a paired channel with none shown.
      if (n_iv == 0L) add(sg, "iv", 1, if (is.na(r$exposures) || !nzchar(r$exposures)) "N/A" else "—", grey = TRUE)
      if (is.na(r$trend_defined)) {
        add(sg, "trend", 1, "—", grey = TRUE)
      } else if (!isTRUE(r$trend_defined)) {
        add(sg, "trend", 1, "n/d", grey = TRUE)
      } else {
        ctr <- contrib(sg, "^s\\(trend\\)$")
        add(sg, "trend", 1, sv(formatC(r$trend, format = "f", digits = if (startsWith(sg, "gcs_")) 1 else dig(sg), flag = "+"), ctr), ctr)
      }
    }
    tb <- do.call(rbind, rows)
    # Reconciliation: on every measured channel the quiet-day evidence plus the
    # contributions shown on its row must equal the channel's evidence, so the
    # printed multipliers multiply to the channel's (up to rounding). An
    # unexposed intervention left off the row contributes under 0.005 nats.
    for (sg in per$sig[per$measured]) {
      gap <- sum(shown[[sg]]) - per$L_full[per$sig == sg]
      if (abs(gap) > 0.02) stop(sprintf("row of %s does not reconcile with its channel evidence (gap %.3f nats)", sg, gap))
    }
    # Offset of the lines of a row from its centre, in row units: the low and high
    # sides, or the first and second intervention line, of one channel; three
    # intervention lines spread their row over three offsets in smaller type.
    LINE_DY <- 0.22; LINE_DY3 <- 0.3
    nl <- tapply(tb$line, tb$sig, max)
    tb$nudge <- ifelse(nl[tb$sig] == 1, 0,
                ifelse(nl[tb$sig] == 2, ifelse(tb$line == 1, LINE_DY, -LINE_DY),
                       c(LINE_DY3, 0, -LINE_DY3)[tb$line]))
    COLS <- c(pct = 0, side = 19, prop = 25, hours = 38, tail = 51, trend = 72, iv = 84)
    HDR  <- c("Percentile of\nChannel score", "Side", "Deviation\nPropensity", "Deviated h /\nObserved h",
              "Tail", "Trend /24 h", "Interventions\n(% of hours exposed, or # of hours dosed)")
    XMAX <- 127
    tb$x <- COLS[tb$col]
    tb$signal <- factor(lab_of[tb$sig], levels(per$signal))
    tb$domain <- dom_of[tb$sig]
    tb$text <- ifelse(tb$bold, paste0("<b>", tb$text, "</b>"), tb$text)
    TXT <- 2.6; TXT3 <- 2.15
    tb$size <- ifelse(nl[tb$sig] == 3, TXT3, TXT)
    lay <- function(d, dy) ggtext::geom_richtext(data = d, aes(x = x, y = signal, label = text, colour = grey, size = size),
                                     hjust = 0, family = "mono", fill = NA, label.colour = NA,
                                     label.padding = grid::unit(0, "pt"), position = position_nudge(y = dy))
    pT <- ggplot(per, aes(y = signal)) +
      geom_blank() +
      geom_vline(xintercept = COLS[-1] - 0.8, colour = "#dedbd4", linewidth = 0.3) +
      lapply(sort(unique(tb$nudge)), function(dy) lay(tb[tb$nudge == dy, ], dy)) +
      scale_size_identity() +
      scale_colour_manual(values = c("FALSE" = INK, "TRUE" = "#5e5c57")) +
      scale_x_continuous(limits = c(0, XMAX), expand = c(0, 0), position = "top", breaks = COLS, labels = HDR) +
      coord_cartesian(clip = "off") +
      # The key, rewritten by the author 2026-10-28: the two bold rules, each
      # with the reference its multiplier is against; the ‡ line only where the
      # symbol appears (author, 2026-10-28). The other columns explain themselves or are keyed in
      # the figure caption.
      labs(x = NULL, y = NULL, title = "Random variables of the stay",
           subtitle = paste0(
             "Bold value: its likelihood multiplier (×) toward death is 2 or more, against a quiet day on the channel.\n",
             "Bold percentile: 85th or above among MIMIC-IV training stays, with the channel's multiplier against the training stays measured on that channel.\n",
             "Grey: every other multiplier. The quiet-day multiplier times the multipliers of a row's values gives the channel's multiplier.",
             if (any(per$measured & per$heavy_shrinkage %in% TRUE)) "\n‡ The prior outweighs the stay's own observed hours." else "")) +
      thm + dom_facet() +
      theme(axis.text.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            strip.text.y = element_blank(),
            axis.text.x.top = element_text(hjust = 0, vjust = 0, colour = INK, face = "bold", size = 8.5))
  }
  # Header text, one line per annotation so the quartile line can be bold.
  # Order: quartile first (bold), then decile, then the 20-bin stratum.
  hdr <- data.frame(
    # The evidence-toward-death/survival line was removed by the author
    # (2026-10-28) together with claim B: distracting, and it showed nothing concrete.
    txt = c(sprintf("Composite score (sum of channel evidence)  %.2f nats  [95%% band %.2f to %.2f]", sm$llr_sum, sm$sum_lo, sm$sum_hi),
            sprintf("  percentile among MIMIC-IV training stays: %.0f (band: %.0f to %.0f)", sm$llr_sum_pct, sm$sum_pct_lo, sm$sum_pct_hi),
            "Reference stratum (MIMIC-IV training, frozen cut points):",
            sprintf("  quartile %d of 4, %s: mortality %.1f%% (%.1f-%.1f%%)", sm$quartile, sm$quartile_label,
                    100 * sm$quartile_rate, 100 * sm$quartile_rate_lo, 100 * sm$quartile_rate_hi),
            sprintf("  decile %d of 10: %.1f%% (%.1f-%.1f%%)", sm$decile, 100 * sm$decile_rate, 100 * sm$decile_rate_lo, 100 * sm$decile_rate_hi),
            sprintf("  bin %d of 20: %.1f%% (%.1f-%.1f%%)   band spans bins %d to %d", sm$bin20, 100 * sm$bin20_rate, 100 * sm$bin20_rate_lo, 100 * sm$bin20_rate_hi, sm$bin20_lo, sm$bin20_hi),
            sprintf("Cohort rate %.1f%%.  %d of 19 channels measured.  %d draws.", 100 * sm$cohort_rate, sm$n_measured, sm$n_draws)),
    face = c("plain", "plain", "plain", "bold", "plain", "plain", "plain"))
  # The Masking State lines follow the score lines; one header line is 0.1
  # layout units, so the table keeps its size as the header grows.
  hdr <- rbind(hdr, data.frame(txt = mask_lines(per), face = "plain"))
  hdr$y <- 1 - (seq_len(nrow(hdr)) - 1) / nrow(hdr)
  pB <- ggplot(hdr) + geom_text(aes(x = 0, y = y, label = txt, fontface = face), hjust = 0, vjust = 1,
                                size = 3.2, colour = INK, family = "mono") +
    xlim(0, 1) + ylim(0, 1) + theme_void() + labs(title = title) +
    theme(plot.title = element_text(colour = INK, face = "bold", size = 11))
  if (is.null(pT)) return(pB / (pA | pP) + plot_layout(heights = c(1, 3.2), widths = c(2.2, 1)))
  pB / ((pA | pT) + plot_layout(widths = c(1.45, 3.6))) + plot_layout(heights = c(0.1 * nrow(hdr), 4.2))
}
# MASKING STATE (author, 2026-10-28): one line per distinct state of
# run/patient_card.R state_of(), naming the channels it applies to, in the
# card's domain order; "N/A" when no channel has one. Wrapped at MASK_WRAP
# characters with a hanging indent.
MASK_WRAP <- 130
mask_lines <- function(per) {
  # Channel codes: draw_card() keeps them in `sig` after relabeling `signal`.
  sg <- as.character(if ("sig" %in% names(per)) per$sig else per$signal)
  st <- if ("state" %in% names(per)) per$state else rep(NA_character_, nrow(per))
  st[is.na(st)] <- ""          # read.csv reads an all-empty column as NA
  o <- match(DOM$signal, sg); sg <- sg[o]; st <- st[o]
  names_of <- function(s) {
    and_join <- function(x) if (length(x) == 1L) x else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
    if (all(startsWith(s, "gcs_"))) return(paste("GCS", and_join(sub("^gcs_", "", s))))
    and_join(tolower(unname(CH_LAB[s])))
  }
  if (!any(nzchar(st))) return("Masking State: N/A")
  c("Masking State:", unlist(lapply(unique(st[nzchar(st)]), function(s) {
    strwrap(paste0(s, if (grepl("removed$", s)) " from " else " for ", names_of(sg[st == s])),
            width = MASK_WRAP, indent = 2, exdent = 4)
  })))
}
# Size of a card: with the random-variable column the card is set landscape. A
# header line is 0.1 of the 5.0 layout units the 9.8-inch card had with its
# eight score lines, so each header line is 0.196 inch; the card now has seven
# score lines (the evidence-balance line removed), plus its Masking State lines.
CARD_W <- function(per) if ("pi_minus" %in% names(per)) 15 else 9
CARD_H <- function(per) (if ("pi_minus" %in% names(per)) 9.8 - 0.196 else 7.5) + 0.196 * length(mask_lines(per))

tabs_dir <- file.path(RUN, "tables")
if (file.exists(file.path(tabs_dir, "per_signal.csv"))) {
  per <- read.csv(file.path(tabs_dir, "per_signal.csv"), stringsAsFactors = FALSE)
  sm  <- read.csv(file.path(tabs_dir, "summary.csv"), stringsAsFactors = FALSE)
  ctf <- file.path(tabs_dir, "card_terms.csv")
  ct  <- if (file.exists(ctf)) read.csv(ctf, stringsAsFactors = FALSE) else NULL
  sc  <- yaml::read_yaml(file.path(RUN, "manifest.yml"))$state_case
  ttl <- sprintf("Patient card, %s%s", sm$site, if (length(sc) && nzchar(sc)) paste0(": state case ", gsub("_", " ", sc)) else "")
  save_both(draw_card(per, sm, ttl, ct), "F14_patient_card", CARD_W(per), CARD_H(per))
  quit(save = "no")
}

profiles <- names(PROFILE_LAB)[file.exists(file.path(tabs_dir, paste0("per_signal_", names(PROFILE_LAB), ".csv")))]
if (length(profiles) != 3L) stop("expected three profile cards under ", tabs_dir)
for (tag in profiles) {
  per <- read.csv(file.path(tabs_dir, paste0("per_signal_", tag, ".csv")), stringsAsFactors = FALSE)
  sm  <- read.csv(file.path(tabs_dir, paste0("summary_", tag, ".csv")), stringsAsFactors = FALSE)
  ctf <- file.path(tabs_dir, paste0("card_terms_", tag, ".csv"))
  ct  <- if (file.exists(ctf)) read.csv(ctf, stringsAsFactors = FALSE) else NULL
  save_both(draw_card(per, sm, sprintf("Patient card, %s: %s", sm$site, PROFILE_LAB[[tag]]), ct),
            paste0("F14_patient_card_", tag), CARD_W(per), CARD_H(per))
}

# --- F17: the three stays on the fitted terms, one figure per decomposed signal --
tm_all <- do.call(rbind, lapply(profiles, function(tag) {
  z <- read.csv(file.path(tabs_dir, paste0("terms_", tag, ".csv")), stringsAsFactors = FALSE)
  z$profile <- factor(PROFILE_LAB[[tag]], PROFILE_LAB); z
}))
# Anchored evidence functions (paper/make_anchored_curves.R): each term relative
# to its reference value, matching the anchored term values of the card.
cv_all <- readRDS("paper/figs/coefs/anchored_curves.rds")
ANCH   <- read.csv("paper/figs/coefs/anchors.csv", stringsAsFactors = FALSE)
pc_all <- read.csv("paper/figs/coefs/parametric_coefs.csv", stringsAsFactors = FALSE)
TERM_LAB <- c("s(pi_minus)" = "Deviation propensity, low side", "s(pi_plus)" = "Deviation propensity, high side",
              "s(q05_delta)" = "Frequency-adjusted tail extremity, low side", "s(q95_delta)" = "Frequency-adjusted tail extremity, high side",
              "s(trend)" = "Trend", "s(invasive_vent__exposure_frac)" = "Invasive ventilation, exposure fraction",
              "invasive_vent__present_at_admission" = "Invasive ventilation, admission indicator",
              "s(vasopressor__exposure_frac)" = "Vasopressor, exposure fraction",
              "s(vasopressor__lambda)" = "Vasopressor, exposure-adjusted intensity",
              "vasopressor__present_at_admission" = "Vasopressor, admission indicator")
Y_LAB <- "Relative to reference (nats)"
UNIT <- c(resp_rate = "breaths/min", heart_rate = "beats/min")
x_label <- function(sg, t) {
  u <- if (sg %in% names(UNIT)) UNIT[[sg]] else "units"
  if (t == "s(pi_minus)") return("Propensity to fall below the reference range")
  if (t == "s(pi_plus)")  return("Propensity to rise above the reference range")
  if (endsWith(t, "_delta)")) return(sprintf("Extremity relative to expected (%s)", u))
  if (t == "s(trend)") return(sprintf("Slope of hourly medians (%s per 24 h)", u))
  if (endsWith(t, "__exposure_frac)")) return("Fraction of the 24-hour window exposed")
  if (endsWith(t, "__lambda)")) return("Intensity relative to what exposure predicts")
  NULL
}
# A term shrunk to (near) zero must look flat: every panel spans at least
# this many nats on its y axis, centred on its own data.
MIN_SPAN <- 0.5
# Propensities first, then tail extremities, then trend, then the intervention terms:
# a tail extremity is read after the propensity on its side.
rank_of <- function(t) ifelse(startsWith(t, "s(pi_"), 1, ifelse(endsWith(t, "_delta)"), 2, ifelse(t == "s(trend)", 3, 4)))

for (sg in unique(tm_all$signal)) {
  tm  <- tm_all[tm_all$signal == sg, ]
  key <- paste(sg, tm$model[1], sep = "/")
  cv  <- cv_all[cv_all$key == key, ]
  pc  <- pc_all[pc_all$key == key, ]
  # Consistency check between the two sources: the QC curve evaluated at a
  # stay's value must equal the stay's term contribution from the bundle. One
  # aggregate number is printed (the largest disagreement), no stay values.
  sm_rows <- tm[tm$term %in% cv$term, ]
  dev <- vapply(seq_len(nrow(sm_rows)), function(i) {
    cz <- cv[cv$term == sm_rows$term[i], ]; cz <- cz[order(cz$grid_x), ]
    if (sm_rows$x[i] < min(cz$grid_x) || sm_rows$x[i] > max(cz$grid_x)) return(NA_real_)
    abs(stats::approx(cz$grid_x, cz$fit, xout = sm_rows$x[i])$y - sm_rows$contribution[i])
  }, numeric(1))
  cat(sprintf("%s: curve vs card term contribution: max |difference| %.4f nats over %d in-grid term values (%d off grid)\n",
              sg, max(dev, na.rm = TRUE), sum(!is.na(dev)), sum(is.na(dev))))

  terms_order <- unique(tm$term); terms_order <- terms_order[order(rank_of(terms_order), seq_along(terms_order))]
  panels <- lapply(terms_order, function(t) {
    pts <- tm[tm$term == t, ]
    if (t %in% cv$term) {
      cz <- cv[cv$term == t, ]; cz <- cz[order(cz$grid_x), ]
      cz$lo <- cz$fit - 1.96 * cz$se; cz$hi <- cz$fit + 1.96 * cz$se
      x0 <- ANCH$anchor[ANCH$key == key & ANCH$term == t][1]
      g <- ggplot() + geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
        geom_vline(xintercept = x0, colour = INK2, linewidth = 0.3, linetype = "22") +
        geom_ribbon(data = cz, aes(grid_x, ymin = lo, ymax = hi), fill = BLUE, alpha = 0.12)
      g <- if (cz$kind[1] == "atoms") g + geom_point(data = cz, aes(grid_x, fit), colour = BLUE, size = 1.1, alpha = 0.6)
           else g + geom_line(data = cz, aes(grid_x, fit), colour = BLUE, linewidth = 0.6, alpha = 0.8)
    } else {
      cf <- pc[pc$term == t, ]
      if (nrow(cf) != 1L) stop("no parametric coefficient for ", key, " ", t)
      base <- data.frame(grid_x = c(0, 1), fit = c(0, cf$estimate), lo = c(0, cf$estimate - 1.96 * cf$se),
                         hi = c(0, cf$estimate + 1.96 * cf$se))
      g <- ggplot() + geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
        geom_errorbar(data = base, aes(grid_x, ymin = lo, ymax = hi), width = 0.06, colour = BLUE, alpha = 0.5) +
        scale_x_continuous(breaks = c(0, 1), labels = c("not at admission", "at admission"), limits = c(-0.3, 1.3))
    }
    g + geom_errorbar(data = pts, aes(x, ymin = lo, ymax = hi, colour = profile), width = 0, linewidth = 1) +
      geom_point(data = pts, aes(x, contribution, colour = profile), size = 2.4) +
      scale_colour_manual(values = PROFILE_COL, drop = FALSE) +
      expand_limits(y = c(-MIN_SPAN / 2, MIN_SPAN / 2)) +
      labs(title = if (t %in% names(TERM_LAB)) TERM_LAB[[t]] else t, x = x_label(sg, t), y = Y_LAB) + thm +
      theme(plot.title = element_text(size = 9))
  })
  refe <- tm$reference_evidence[1]
  p17 <- wrap_plots(panels, ncol = 3, guides = "collect") +
    plot_annotation(subtitle = sprintf(paste0("Each term relative to its reference value (dotted line). Reference evidence of the ",
                                              "channel: %+.2f nats; a stay's evidence for the channel is this plus its terms."), refe),
                    theme = theme(plot.subtitle = element_text(size = 9, colour = INK2))) &
    theme(legend.position = "bottom", legend.title = element_blank())
  save_both(p17, paste0("F17_card_terms_", sg), 9.5, 3.2 * ceiling(length(panels) / 3) + 0.6)
}
