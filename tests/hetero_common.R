# tests/hetero_common.R --------------------------------------------------------
# Shared declarations and helpers for the channel heterogeneity stages
# (docs/heterogeneity_plan_20261027.md). Sourced by tests/hetero_*.R; defines
# nothing that fits.
# ------------------------------------------------------------------------------

# The percentile levels of a channel's own out-of-fold weights of evidence among
# measured MIMIC-IV training stays (author, 2026-10-27). High levels mark a
# channel extreme toward death, low levels extreme toward survival.
HETERO_HIGH <- c(80, 90)
HETERO_LOW  <- c(5, 10)
# The four settings; the first is primary for the confirmation rule.
HETERO_SETTINGS <- data.frame(high = c(80, 80, 90, 90), low = c(5, 10, 10, 5),
                              primary = c(TRUE, FALSE, FALSE, FALSE))
HETERO_MIN_CELL <- 10L          # no rate or quantile from fewer stays (DUA floor)
HETERO_SITES <- c("mimic_oof", "eicu")

.hopt <- function(args, flag, default) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L]
}

#' The stage 1 objects for one site.
hetero_read <- function(stage1, site) {
  if (!site %in% HETERO_SITES) stop("unknown site: ", site, call. = FALSE)
  chk <- read.csv(file.path(stage1, "diagnostics", "hetero_assemble_checks.csv"), stringsAsFactors = FALSE)
  if (any(chk$status != "PASS")) stop("stage 1 run ", stage1, " has failed checks", call. = FALSE)
  rd <- function(nm) readRDS(file.path(stage1, "tables", paste0(nm, ".rds")))
  # The measured matrix comes from the stage A run that stage 1 read, by its
  # manifest, so the two can never be from different runs.
  staytab <- yaml::read_yaml(file.path(stage1, "manifest.yml"))$staytab
  if (is.null(staytab)) stop("stage 1 manifest names no staytab run", call. = FALSE)
  out <- list(tab = rd(paste0("hetero_stay_", site)), ext = rd(paste0("hetero_extreme_", site)),
              quartile_breaks = rd("quartile_breaks"),
              measured = readRDS(file.path(staytab, "tables", paste0("measured_", site, ".rds"))))
  if (!identical(out$tab$stay_id, rownames(out$measured)) && !is.null(rownames(out$measured))) {
    stop("stage 1 and stage A rows are not aligned at ", site, call. = FALSE)
  }
  if (!identical(colnames(out$measured), colnames(out$ext[[1]]))) stop("channel order differs", call. = FALSE)
  out
}

#' Deaths-minus-survivors difference in the mean of each column of X (0/1 or
#' numeric), with a patient-bootstrap percentile interval. The draws resample
#' patients within the rows given, so the interval is conditional on the
#' stratum, as the per-bin intervals elsewhere are.
#' @return data frame, one row per column of X
boot_diff <- function(X, y, group, n_boot, seed) {
  X <- as.matrix(X); y <- as.integer(y)
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  est_of <- function(j) colMeans(X[j[y[j] == 1L], , drop = FALSE]) - colMeans(X[j[y[j] == 0L], , drop = FALSE])
  out <- data.frame(measure = colnames(X), deaths = n1, survivors = n0,
                    mean_deaths = NA_real_, mean_survivors = NA_real_, difference = NA_real_,
                    lo = NA_real_, hi = NA_real_, n_boot = n_boot, boot_unit = NA_character_,
                    stringsAsFactors = FALSE)
  if (n1 < HETERO_MIN_CELL || n0 < HETERO_MIN_CELL) return(out)
  out$mean_deaths <- colMeans(X[y == 1L, , drop = FALSE])
  out$mean_survivors <- colMeans(X[y == 0L, , drop = FALSE])
  out$difference <- out$mean_deaths - out$mean_survivors
  bu <- .boot_units(length(y), group)
  B <- with_seed(seed, vapply(seq_len(n_boot), function(b) {
    j <- bu$draw()
    if (sum(y[j]) == 0L || sum(y[j]) == length(j)) return(rep(NA_real_, ncol(X)))
    est_of(j)
  }, numeric(ncol(X))))
  B <- matrix(B, nrow = ncol(X))
  out$lo <- apply(B, 1, stats::quantile, 0.025, na.rm = TRUE)
  out$hi <- apply(B, 1, stats::quantile, 0.975, na.rm = TRUE)
  out$boot_unit <- bu$unit
  out
}

#' Quantiles of x with the floor applied.
q_floor <- function(x, probs = c(0.25, 0.5, 0.75)) {
  if (length(x) < HETERO_MIN_CELL) return(stats::setNames(rep(NA_real_, length(probs)), paste0("q", probs * 100)))
  stats::setNames(stats::quantile(x, probs, names = FALSE), paste0("q", probs * 100))
}

#' Strata of one partition: list of row-index vectors named by stratum.
hetero_strata <- function(tab, partition) {
  v <- tab[[partition]]
  split(seq_len(nrow(tab)), factor(v, levels = sort(unique(v))))
}

#' Displayed items of one channel on the patient card and their contributions,
#' counted as paper/make_card_figure.R draws them (moved here from
#' tests/hetero_card_flags.R on 2026-10-28, unchanged, so the card-flag analysis
#' and tests/bold_values_replicates.R count with one definition). `tt` holds the
#' anchored term columns and the display columns ("x:<iv>__..."), one row per
#' stay measured on the channel. Items: each modeled propensity side; each
#' modeled tail, or one combined tail for a single observed hour on two sides;
#' the trend when established; each paired intervention shown (terms summed).
#' @param strong the bold threshold in nats (log 2 on the card), passed by the caller
#' @param tails_eq per row: whether the low and high tail values are one reading
#'   (tails_equal()). The card shows one "singular" tail only then
#'   (paper/make_card_figure.R: n_obs == 1 AND tail_low == tail_high AND two sides);
#'   with one observed hour, a raw minimum and maximum can still differ. Required
#'   since 2026-10-30 (Codex review, finding 1): the first version merged the
#'   tails on n_obs == 1 alone.
#' @return matrix (n_items, n_bold) per row; attribute "items" holds the item
#'   contributions, NA where an item is not shown, in a fixed column order
#' Whether a stay's low and high tail values are one reading, as the card tests
#' it: the statistic family of run/patient_card.R (q05/q95 when the channel's
#' modeled tail is a quantile, value_min/value_max otherwise), read from the
#' signal_features rows `z` (columns q05, q95, value_min, value_max).
tails_equal <- function(z, sg, cfg) {
  dv  <- grep("_delta$", all.vars(build_formula(sg, "meas", cfg)), value = TRUE)
  fam <- if (any(delta_base_of(dv) %in% c("q05", "q95"))) c("q05", "q95") else c("value_min", "value_max")
  z[[fam[1]]] == z[[fam[2]]]
}

channel_items <- function(tt, sg, n_obs, trend_def, cfg, strong, tails_eq) {
  cols <- names(tt); it <- list()
  pis <- intersect(c("s(pi_minus)", "s(pi_plus)"), cols)
  for (p in pis) it[[p]] <- tt[[p]]
  tl <- grep("_delta\\)$", cols, value = TRUE)
  if (length(tails_eq) != length(n_obs)) stop("channel_items: tails_eq must have one value per row", call. = FALSE)
  single <- n_obs == 1 & (tails_eq %in% TRUE) & length(pis) == 2L
  if (length(tl)) {
    it[["tail_singular"]] <- ifelse(single, rowSums(tt[, tl, drop = FALSE]), NA_real_)
    for (v in tl) it[[v]] <- ifelse(single, NA_real_, tt[[v]])
  }
  if ("s(trend)" %in% cols) it[["trend"]] <- ifelse(trend_def %in% TRUE, tt[["s(trend)"]], NA_real_)
  for (iv in interventions_of(sg, cfg)) {
    ic <- grep(paste0("(^|\\()", iv, "__"), cols, value = TRUE)
    if (!length(ic)) next
    civ <- rowSums(tt[, ic, drop = FALSE])
    val <- tt[[paste0("x:", iv, "__exposure_frac")]]
    if (is.null(val)) val <- tt[[paste0("x:", iv, "__n_hours")]]
    adm <- tt[[paste0("x:", iv, "__present_at_admission")]]
    if (is.null(adm)) adm <- 0
    keep <- !(val == 0 & !(adm %in% 1) & abs(civ) < 0.005)
    it[[iv]] <- ifelse(keep, civ, NA_real_)
  }
  M <- do.call(cbind, it)
  out <- cbind(n_items = rowSums(!is.na(M)), n_bold = rowSums(!is.na(M) & M >= strong))
  attr(out, "items") <- M
  out
}
